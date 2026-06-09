args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
} else {
  file.path(getwd(), "scripts", "run_qr_adaptation_comparison.R")
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

source(file.path(root, "R", "prox.R"))
source(file.path(root, "R", "objective.R"))
source(file.path(root, "R", "simulate.R"))
source(file.path(root, "R", "qr_pdhg.R"))
source(file.path(root, "R", "fedspd_dp_qr.R"))

dir.create(file.path(root, "results"), recursive = TRUE, showWarnings = FALSE)

run_case <- function(dat, client_indices, tau, loss, K, label,
                     rounds = 250, Q = 5) {
  message(sprintf("Running tau=%.2f loss=%s K=%d", tau, loss, K))
  fit <- fedspd_dp_qr(
    dat$X, dat$y,
    client_indices = client_indices,
    tau = tau,
    smooth_mu = 0.05,
    loss = loss,
    rounds = rounds,
    Q = Q,
    clients_per_round = K,
    batch_size = 30,
    rho = 20,
    gamma_rule = "sqrt",
    gamma0 = 4,
    dp = FALSE,
    beta_ref = dat$beta,
    trace_every = 5,
    seed = 20260526
  )

  trace <- fit$trace
  trace$tau <- tau
  trace$loss <- loss
  trace$K <- K
  trace$label <- label

  summary <- data.frame(
    tau = tau,
    loss = loss,
    label = label,
    K = K,
    Q = Q,
    rounds = rounds,
    final_qr_objective = fit$qr_objective,
    final_smooth_objective = fit$smooth_objective,
    final_consensus_gap = fit$consensus_gap,
    final_beta_l2_error = sqrt(sum((fit$beta - dat$beta)^2))
  )

  list(summary = summary, trace = trace, fit = fit)
}

taus <- c(0.5, 0.75, 0.9)
losses <- c("smooth", "check")
n_clients <- 20

runs <- list()
k <- 1
central_rows <- list()

for (tau in taus) {
  dat <- make_qr_sim(n = 800, p = 10, tau = tau, noise = "asymmetric", seed = 20260526)
  client_indices <- iid_partition(nrow(dat$X), n_clients = n_clients, seed = 20260526)

  central <- qr_pdhg(
    dat$X, dat$y,
    tau = tau,
    max_iter = 2500,
    step_rule = "box",
    trace_every = 50,
    seed = 20260526
  )

  central_rows[[length(central_rows) + 1]] <- data.frame(
    tau = tau,
    loss = "central",
    label = "central_qr_pdhg",
    K = NA_integer_,
    Q = NA_integer_,
    rounds = 2500,
    final_qr_objective = central$objective,
    final_smooth_objective = smooth_qr_objective(
      dat$X, dat$y, central$beta,
      tau = tau,
      mu = 0.05
    ),
    final_consensus_gap = NA_real_,
    final_beta_l2_error = sqrt(sum((central$beta - dat$beta)^2))
  )

  for (loss in losses) {
    for (case in list(
      list(K = n_clients, label = "full_clients"),
      list(K = ceiling(n_clients / 2), label = "half_clients")
    )) {
      runs[[k]] <- run_case(
        dat, client_indices,
        tau = tau,
        loss = loss,
        K = case$K,
        label = case$label
      )
      k <- k + 1
    }
  }
}

summary_tbl <- do.call(rbind, lapply(runs, `[[`, "summary"))
trace_tbl <- do.call(rbind, lapply(runs, `[[`, "trace"))
summary_tbl <- rbind(summary_tbl, do.call(rbind, central_rows))

write.csv(summary_tbl,
          file.path(root, "results", "qr_adaptation_comparison_summary.csv"),
          row.names = FALSE)
write.csv(trace_tbl,
          file.path(root, "results", "qr_adaptation_comparison_trace.csv"),
          row.names = FALSE)

cat("\nQR adaptation comparison summary:\n")
print(summary_tbl)

