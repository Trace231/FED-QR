args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
} else {
  file.path(getwd(), "scripts", "run_nyc_taxi_large_scale_experiment.R")
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

source(file.path(root, "scripts", "load_rfedqr.R"))
make_experiment_dirs(root)

data_path <- file.path(root, "data", "processed", "nyc_taxi_yellow_qr_2024q1_sample.csv")
if (!file.exists(data_path)) {
  stop("NYC Taxi processed data not found. Run scripts/prepare_nyc_taxi_data.sh first.")
}

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) as.integer(value) else default
}

env_num <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) as.numeric(value) else default
}

max_rows <- env_int("NYC_MODEL_MAX_ROWS", NA_integer_)
rounds <- env_int("NYC_ROUNDS", 350)
trace_every <- env_int("NYC_TRACE_EVERY", 25)
clients_per_round_cap <- env_int("NYC_CLIENTS_PER_ROUND", 12)
batch_size <- env_int("NYC_BATCH_SIZE", 200)
seed <- env_int("NYC_SEED", 20260526)
tau <- env_num("NYC_TAU", 0.9)

message("Reading NYC Taxi modeling CSV...")
taxi <- read.csv(data_path)
if (!is.na(max_rows) && nrow(taxi) > max_rows) {
  set.seed(seed)
  taxi <- taxi[sample.int(nrow(taxi), max_rows), ]
}

taxi <- taxi[complete.cases(taxi), ]
taxi$client_id <- as.integer(taxi$client_id)

numeric_features <- c(
  "log_trip_distance",
  "duration_min",
  "passenger_count",
  "congestion_surcharge",
  "airport_fee"
)
for (v in numeric_features) {
  taxi[[paste0(v, "_z")]] <- as.numeric(scale(taxi[[v]]))
}

taxi$duration_log_z <- as.numeric(scale(log1p(taxi$duration_min)))
taxi$hour_sin <- sin(2 * pi * taxi$pickup_hour / 24)
taxi$hour_cos <- cos(2 * pi * taxi$pickup_hour / 24)
taxi$dow_sin <- sin(2 * pi * taxi$pickup_dow / 7)
taxi$dow_cos <- cos(2 * pi * taxi$pickup_dow / 7)
taxi$pickup_borough <- factor(taxi$pickup_borough)
taxi$vendor_id <- factor(taxi$vendor_id)
taxi$rate_code <- factor(taxi$rate_code)
taxi$payment_type <- factor(taxi$payment_type)

formula_nyc <- log_total_amount ~
  log_trip_distance_z + duration_log_z + passenger_count_z +
  congestion_surcharge_z + airport_fee_z +
  hour_sin + hour_cos + dow_sin + dow_cos +
  pickup_borough + vendor_id + rate_code + payment_type

X <- model.matrix(formula_nyc, data = taxi)
y <- taxi$log_total_amount
client_raw <- split(seq_along(y), taxi$client_id)
client_indices <- client_raw[lengths(client_raw) >= 100]
client_indices <- client_indices[order(as.integer(names(client_indices)))]

keep_idx <- sort(unlist(client_indices, use.names = FALSE))
X <- X[keep_idx, , drop = FALSE]
y <- y[keep_idx]
client_id_keep <- taxi$client_id[keep_idx]
client_indices <- split(seq_along(y), client_id_keep)
client_indices <- client_indices[order(as.integer(names(client_indices)))]

n_clients <- length(client_indices)
clients_per_round <- min(clients_per_round_cap, n_clients)

design_summary <- data.frame(
  n = nrow(X),
  p = ncol(X),
  tau = tau,
  n_clients = n_clients,
  min_client_n = min(lengths(client_indices)),
  median_client_n = median(lengths(client_indices)),
  max_client_n = max(lengths(client_indices)),
  clients_per_round = clients_per_round,
  batch_size = batch_size,
  rounds = rounds
)
write.csv(design_summary,
          file.path(root, "results", "nyc_taxi_large_scale_design.csv"),
          row.names = FALSE)

client_summary <- do.call(rbind, lapply(names(client_indices), function(cid) {
  idx <- client_indices[[cid]]
  data.frame(
    client_id = as.integer(cid),
    n = length(idx),
    y_mean = mean(y[idx]),
    y_q50 = as.numeric(quantile(y[idx], 0.5)),
    y_q90 = as.numeric(quantile(y[idx], 0.9)),
    distance_mean = mean(taxi$trip_distance[keep_idx][idx]),
    duration_mean = mean(taxi$duration_min[keep_idx][idx])
  )
}))
write.csv(client_summary,
          file.path(root, "results", "nyc_taxi_large_scale_clients.csv"),
          row.names = FALSE)

methods <- c("QR box-dual", "QR box-dual long", "FSPG-smooth", "FedSubGrad", "FedSPD-check")
result <- run_fedqr_methods(
  methods,
  X,
  y,
  client_indices = client_indices,
  tau = tau,
  rounds = rounds,
  clients_per_round = clients_per_round,
  batch_size = batch_size,
  seed = seed,
  trace_every = trace_every,
  term_names = colnames(X)
)
summary_tbl <- result$summary
trace_tbl <- result$trace
save_experiment_outputs("nyc_taxi_large_scale", result, root)

method_cols <- c(
  "QR box-dual" = "#2A6FBB",
  "QR box-dual long" = "#1B9E77",
  "FSPG-smooth" = "#984EA3",
  "FedSubGrad" = "#4DAF4A",
  "FedSPD-check" = "#D95F02"
)

plot_convergence(
  trace_tbl,
  file.path(root, "figures", "nyc_taxi_large_scale_convergence.png"),
  title = sprintf("NYC Taxi large-scale federated QR, tau=%.2f", tau),
  method_cols = method_cols
)
plot_final_gap(
  summary_tbl,
  file.path(root, "figures", "nyc_taxi_large_scale_final_gap.png"),
  title = sprintf("NYC Taxi final objective, n=%s, clients=%d",
                  format(nrow(X), big.mark = ","), n_clients),
  method_cols = method_cols
)

cat("\nNYC Taxi design summary:\n")
print(design_summary)
cat("\nNYC Taxi method summary:\n")
print(summary_tbl)
cat("\nFigures written to:\n")
cat(file.path(root, "figures", "nyc_taxi_large_scale_convergence.png"), "\n")
cat(file.path(root, "figures", "nyc_taxi_large_scale_final_gap.png"), "\n")
