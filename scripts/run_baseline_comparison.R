args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
} else {
  file.path(getwd(), "scripts", "run_baseline_comparison.R")
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

source(file.path(root, "R", "prox.R"))
source(file.path(root, "R", "objective.R"))
source(file.path(root, "R", "simulate.R"))
source(file.path(root, "R", "qr_pdhg.R"))
source(file.path(root, "R", "fed_qr_spd.R"))
source(file.path(root, "R", "baselines.R"))

dir.create(file.path(root, "results"), recursive = TRUE, showWarnings = FALSE)

run_method <- function(method, dat, client_indices, tau, regime, rounds,
                       clients_per_round, batch_size) {
  if (method == "fedspd_generic") {
    fit <- fed_qr_spd(
      dat$X, dat$y,
      client_indices = client_indices,
      tau = tau,
      rounds = rounds,
      clients_per_round = clients_per_round,
      batch_size = batch_size,
      step_rule = "generic",
      beta_true = dat$beta,
      trace_every = 10,
      seed = 20260526
    )
  } else if (method == "qrfedspd_box") {
    fit <- fed_qr_spd(
      dat$X, dat$y,
      client_indices = client_indices,
      tau = tau,
      rounds = rounds,
      clients_per_round = clients_per_round,
      batch_size = batch_size,
      step_rule = "box",
      beta_true = dat$beta,
      trace_every = 10,
      seed = 20260526
    )
  } else if (method == "fed_subgrad") {
    fit <- fed_subgrad_qr(
      dat$X, dat$y,
      client_indices = client_indices,
      tau = tau,
      rounds = rounds,
      clients_per_round = clients_per_round,
      batch_size = batch_size,
      beta_true = dat$beta,
      trace_every = 10,
      seed = 20260526
    )
  } else if (method == "fspg_smooth") {
    fit <- fed_smooth_qr(
      dat$X, dat$y,
      client_indices = client_indices,
      tau = tau,
      rounds = rounds,
      clients_per_round = clients_per_round,
      batch_size = batch_size,
      smooth_mu = 0.1,
      beta_true = dat$beta,
      trace_every = 10,
      seed = 20260526
    )
  } else {
    stop("Unknown method: ", method)
  }

  trace <- fit$trace
  trace$method <- method
  trace$tau <- tau
  trace$regime <- regime

  summary <- data.frame(
    tau = tau,
    regime = regime,
    method = method,
    rounds = rounds,
    clients_per_round = clients_per_round,
    batch_size = batch_size,
    final_objective = fit$objective,
    final_beta_l2_error = sqrt(sum((fit$beta - dat$beta)^2))
  )

  list(summary = summary, trace = trace)
}

taus <- c(0.5, 0.9)
methods <- c("fed_subgrad", "fspg_smooth", "fedspd_generic", "qrfedspd_box")
rounds <- 1000
n_clients <- 8

runs <- list()
k <- 1
for (tau in taus) {
  dat <- make_qr_sim(n = 1200, p = 12, tau = tau, noise = "asymmetric", seed = 20260526)
  client_indices <- dirichlet_partition(dat$y, n_clients = n_clients, alpha = 0.25, seed = 20260526)
  regimes <- data.frame(
    regime = c("full_clients_full_batch", "r1_r2_stochastic"),
    clients_per_round = c(n_clients, ceiling(n_clients / 2)),
    batch_size = c(max(lengths(client_indices)), 40)
  )

  for (r in seq_len(nrow(regimes))) {
    for (method in methods) {
      message(sprintf(
        "Running tau=%.2f regime=%s method=%s",
        tau,
        regimes$regime[r],
        method
      ))
      runs[[k]] <- run_method(
        method, dat, client_indices, tau,
        regime = regimes$regime[r],
        rounds = rounds,
        clients_per_round = regimes$clients_per_round[r],
        batch_size = regimes$batch_size[r]
      )
      k <- k + 1
    }
  }
}

summary_tbl <- do.call(rbind, lapply(runs, `[[`, "summary"))
trace_tbl <- do.call(rbind, lapply(runs, `[[`, "trace"))

write.csv(summary_tbl,
          file.path(root, "results", "baseline_comparison_summary.csv"),
          row.names = FALSE)
write.csv(trace_tbl,
          file.path(root, "results", "baseline_comparison_trace.csv"),
          row.names = FALSE)

cat("\nBaseline comparison summary:\n")
print(summary_tbl)
