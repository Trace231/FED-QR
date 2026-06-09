args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
} else {
  file.path(getwd(), "scripts", "run_federated_simulation.R")
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

source(file.path(root, "R", "prox.R"))
source(file.path(root, "R", "objective.R"))
source(file.path(root, "R", "simulate.R"))
source(file.path(root, "R", "qr_pdhg.R"))
source(file.path(root, "R", "fed_qr_spd.R"))

dir.create(file.path(root, "results"), recursive = TRUE, showWarnings = FALSE)

summarize_partition <- function(dat, client_indices, partition_name) {
  do.call(rbind, lapply(seq_along(client_indices), function(j) {
    idx <- client_indices[[j]]
    data.frame(
      partition = partition_name,
      client = j,
      n = length(idx),
      y_mean = mean(dat$y[idx]),
      y_sd = sd(dat$y[idx]),
      y_q10 = as.numeric(quantile(dat$y[idx], 0.10)),
      y_q50 = as.numeric(quantile(dat$y[idx], 0.50)),
      y_q90 = as.numeric(quantile(dat$y[idx], 0.90))
    )
  }))
}

run_one <- function(dat, client_indices, partition_name, config_name,
                    clients_per_round, batch_size, rounds = 800) {
  fit <- fed_qr_spd(
    dat$X, dat$y,
    client_indices = client_indices,
    tau = dat$tau,
    lambda = 0,
    penalty = "none",
    rounds = rounds,
    clients_per_round = clients_per_round,
    batch_size = batch_size,
    step_rule = "box",
    client_weighting = "renormalized",
    beta_true = dat$beta,
    trace_every = 10,
    seed = 20260526
  )

  trace <- fit$trace
  trace$partition <- partition_name
  trace$config <- config_name
  trace$clients_per_round <- clients_per_round
  trace$batch_size <- batch_size

  summary <- data.frame(
    partition = partition_name,
    config = config_name,
    rounds = rounds,
    clients_per_round = clients_per_round,
    batch_size = batch_size,
    final_objective = fit$objective,
    final_beta_l2_error = sqrt(sum((fit$beta - dat$beta)^2)),
    eta = fit$eta,
    sigma = fit$sigma
  )

  list(summary = summary, trace = trace, fit = fit)
}

dat <- make_qr_sim(n = 1200, p = 12, tau = 0.5, noise = "asymmetric", seed = 20260526)
n_clients <- 8

partitions <- list(
  iid = iid_partition(nrow(dat$X), n_clients = n_clients, seed = 20260526),
  noniid = dirichlet_partition(dat$y, n_clients = n_clients, alpha = 0.25, seed = 20260526)
)

partition_summary <- do.call(rbind, Map(
  function(idx, nm) summarize_partition(dat, idx, nm),
  partitions,
  names(partitions)
))
write.csv(partition_summary,
          file.path(root, "results", "federated_partition_summary.csv"),
          row.names = FALSE)

max_client_n <- max(lengths(partitions$iid), lengths(partitions$noniid))
configs <- data.frame(
  config = c("full_clients_full_batch", "r1_only", "r2_only", "r1_r2"),
  clients_per_round = c(n_clients, ceiling(n_clients / 2), n_clients, ceiling(n_clients / 2)),
  batch_size = c(max_client_n, max_client_n, 40, 40)
)

runs <- list()
k <- 1
for (partition_name in names(partitions)) {
  for (i in seq_len(nrow(configs))) {
    message(sprintf(
      "Running partition=%s config=%s",
      partition_name,
      configs$config[i]
    ))
    runs[[k]] <- run_one(
      dat,
      partitions[[partition_name]],
      partition_name = partition_name,
      config_name = configs$config[i],
      clients_per_round = configs$clients_per_round[i],
      batch_size = configs$batch_size[i]
    )
    k <- k + 1
  }
}

summary_tbl <- do.call(rbind, lapply(runs, `[[`, "summary"))
trace_tbl <- do.call(rbind, lapply(runs, `[[`, "trace"))

central <- qr_pdhg(
  dat$X, dat$y,
  tau = dat$tau,
  max_iter = 2000,
  step_rule = "box",
  trace_every = 25,
  seed = 20260526
)

central_summary <- data.frame(
  partition = "central",
  config = "central_qr_pdhg",
  rounds = 2000,
  clients_per_round = NA_integer_,
  batch_size = NA_integer_,
  final_objective = central$objective,
  final_beta_l2_error = sqrt(sum((central$beta - dat$beta)^2)),
  eta = central$eta,
  sigma = central$sigma
)

summary_tbl <- rbind(summary_tbl, central_summary)

write.csv(summary_tbl,
          file.path(root, "results", "federated_simulation_summary.csv"),
          row.names = FALSE)
write.csv(trace_tbl,
          file.path(root, "results", "federated_simulation_trace.csv"),
          row.names = FALSE)

cat("\nPartition summary:\n")
print(partition_summary)

cat("\nFederated simulation summary:\n")
print(summary_tbl)

