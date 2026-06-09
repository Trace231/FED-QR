args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
} else {
  file.path(getwd(), "scripts", "run_household_power_large_scale_experiment.R")
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

source(file.path(root, "scripts", "load_rfedqr.R"))
make_experiment_dirs(root)

data_path <- file.path(root, "data", "processed", "household_power_consumption.csv")
if (!file.exists(data_path)) {
  stop("Household power data not found. Run scripts/prepare_household_power_data.sh first.")
}

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) as.integer(value) else default
}

env_num <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) as.numeric(value) else default
}

max_rows <- env_int("POWER_MODEL_MAX_ROWS", NA_integer_)
rounds <- env_int("POWER_ROUNDS", 400)
trace_every <- env_int("POWER_TRACE_EVERY", 25)
clients_per_round_cap <- env_int("POWER_CLIENTS_PER_ROUND", 10)
batch_size <- env_int("POWER_BATCH_SIZE", 5000)
seed <- env_int("POWER_SEED", 20260602)
tau <- env_num("POWER_TAU", 0.9)
lambda <- env_num("POWER_LAMBDA", 0)
fedspd_rounds <- env_int("POWER_FEDSPD_ROUNDS", rounds)

message("Reading UCI household power CSV...")
power <- read.csv(data_path, na.strings = c("?", ""), stringsAsFactors = FALSE)
names(power) <- make.names(names(power))

numeric_cols <- c(
  "Global_active_power",
  "Global_reactive_power",
  "Voltage",
  "Global_intensity",
  "Sub_metering_1",
  "Sub_metering_2",
  "Sub_metering_3"
)
for (v in numeric_cols) {
  power[[v]] <- as.numeric(power[[v]])
}

power <- power[complete.cases(power[, c("Date", "Time", numeric_cols)]), ]
if (!is.na(max_rows) && nrow(power) > max_rows) {
  set.seed(seed)
  power <- power[sample.int(nrow(power), max_rows), ]
}

date <- as.Date(power$Date, format = "%d/%m/%Y")
hour <- as.integer(substr(power$Time, 1, 2))
minute <- as.integer(substr(power$Time, 4, 5))
dow <- as.POSIXlt(date)$wday
month <- as.integer(format(date, "%m"))
year_month <- format(date, "%Y-%m")

power$log_active_power <- log1p(power$Global_active_power)
power$reactive_z <- as.numeric(scale(power$Global_reactive_power))
power$voltage_z <- as.numeric(scale(power$Voltage))
power$sub1_z <- as.numeric(scale(power$Sub_metering_1))
power$sub2_z <- as.numeric(scale(power$Sub_metering_2))
power$sub3_z <- as.numeric(scale(power$Sub_metering_3))
power$hour_sin <- sin(2 * pi * (hour + minute / 60) / 24)
power$hour_cos <- cos(2 * pi * (hour + minute / 60) / 24)
power$dow_sin <- sin(2 * pi * dow / 7)
power$dow_cos <- cos(2 * pi * dow / 7)
power$month <- factor(month)
power$year_month <- factor(year_month)

formula_power <- log_active_power ~
  reactive_z + voltage_z + sub1_z + sub2_z + sub3_z +
  hour_sin + hour_cos + dow_sin + dow_cos + month

X <- model.matrix(formula_power, data = power)
y <- power$log_active_power

client_raw <- split(seq_along(y), droplevels(power$year_month), drop = TRUE)
client_indices <- client_raw[lengths(client_raw) >= 1000]
client_indices <- client_indices[order(names(client_indices))]

keep_idx <- sort(unlist(client_indices, use.names = FALSE))
X <- X[keep_idx, , drop = FALSE]
y <- y[keep_idx]
power_keep <- power[keep_idx, , drop = FALSE]
client_indices <- split(seq_along(y), droplevels(power_keep$year_month), drop = TRUE)
client_indices <- client_indices[order(names(client_indices))]

n_clients <- length(client_indices)
clients_per_round <- min(clients_per_round_cap, n_clients)

design_summary <- data.frame(
  dataset = "UCI Individual Household Electric Power Consumption",
  n = nrow(X),
  p = ncol(X),
  tau = tau,
  lambda = lambda,
  n_clients = n_clients,
  client_definition = "calendar month",
  min_client_n = min(lengths(client_indices)),
  median_client_n = median(lengths(client_indices)),
  max_client_n = max(lengths(client_indices)),
  clients_per_round = clients_per_round,
  batch_size = batch_size,
  rounds = rounds,
  fedspd_rounds = fedspd_rounds
)
write.csv(design_summary,
          file.path(root, "results", "household_power_large_scale_design.csv"),
          row.names = FALSE)

client_summary <- do.call(rbind, lapply(names(client_indices), function(cid) {
  idx <- client_indices[[cid]]
  data.frame(
    client_id = cid,
    n = length(idx),
    y_mean = mean(y[idx]),
    y_q50 = as.numeric(quantile(y[idx], 0.5)),
    y_q90 = as.numeric(quantile(y[idx], 0.9)),
    active_power_mean = mean(power_keep$Global_active_power[idx]),
    reactive_power_mean = mean(power_keep$Global_reactive_power[idx]),
    voltage_mean = mean(power_keep$Voltage[idx])
  )
}))
write.csv(client_summary,
          file.path(root, "results", "household_power_large_scale_clients.csv"),
          row.names = FALSE)

methods <- c("QR box-dual", "QR box-dual long", "FSPG-smooth", "FedSubGrad", "FedSPD-check")
result <- run_fedqr_methods(
  methods,
  X,
  y,
  client_indices = client_indices,
  tau = tau,
  lambda = lambda,
  rounds = rounds,
  clients_per_round = clients_per_round,
  batch_size = batch_size,
  seed = seed,
  trace_every = trace_every,
  term_names = colnames(X),
  method_controls = list("FedSPD-check" = list(rounds = fedspd_rounds))
)
summary_tbl <- result$summary
trace_tbl <- result$trace
save_experiment_outputs("household_power_large_scale", result, root)

method_cols <- c(
  "QR box-dual" = "#2A6FBB",
  "QR box-dual long" = "#1B9E77",
  "FSPG-smooth" = "#984EA3",
  "FedSubGrad" = "#4DAF4A",
  "FedSPD-check" = "#D95F02"
)

png(file.path(root, "figures", "household_power_client_heterogeneity.png"),
    width = 1500, height = 900, res = 180)
op <- par(mar = c(8, 5, 4, 1))
plot(
  seq_len(nrow(client_summary)),
  client_summary$y_q90,
  type = "b",
  pch = 19,
  col = "#2A6FBB",
  ylim = range(c(client_summary$y_q90, client_summary$y_q50, client_summary$y_mean)),
  xaxt = "n",
  xlab = "Calendar-month client",
  ylab = "Client response quantile / mean",
  main = "Household power monthly-client heterogeneity"
)
lines(seq_len(nrow(client_summary)), client_summary$y_q50, type = "b", pch = 17, col = "#1B9E77")
lines(seq_len(nrow(client_summary)), client_summary$y_mean, type = "b", pch = 15, col = "#984EA3")
axis(1, at = seq_len(nrow(client_summary)), labels = client_summary$client_id, las = 2, cex.axis = 0.6)
legend("topright", legend = c("q90", "q50", "mean"), col = c("#2A6FBB", "#1B9E77", "#984EA3"),
       pch = c(19, 17, 15), lwd = 2, bty = "n")
par(op)
dev.off()

plot_convergence(
  trace_tbl,
  file.path(root, "figures", "household_power_large_scale_convergence.png"),
  title = sprintf("Household power large-scale federated QR, tau=%.2f", tau),
  method_cols = method_cols
)
plot_final_gap(
  summary_tbl,
  file.path(root, "figures", "household_power_large_scale_final_gap.png"),
  title = sprintf("Household power final objective, n=%s, clients=%d",
                  format(nrow(X), big.mark = ","), n_clients),
  method_cols = method_cols
)

cat("\nHousehold power design summary:\n")
print(design_summary)
cat("\nHousehold power method summary:\n")
print(summary_tbl)
cat("\nFigures written to:\n")
cat(file.path(root, "figures", "household_power_client_heterogeneity.png"), "\n")
cat(file.path(root, "figures", "household_power_large_scale_convergence.png"), "\n")
cat(file.path(root, "figures", "household_power_large_scale_final_gap.png"), "\n")
