args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
} else {
  file.path(getwd(), "scripts", "run_fedspd_dp_reproduction.R")
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

source(file.path(root, "R", "prox.R"))
source(file.path(root, "R", "objective.R"))
source(file.path(root, "R", "simulate.R"))
source(file.path(root, "R", "qr_pdhg.R"))
source(file.path(root, "R", "fedspd_dp_qr.R"))

dir.create(file.path(root, "results"), recursive = TRUE, showWarnings = FALSE)

run_case <- function(dat, client_indices, label, K, Q, dp,
                     rounds = 200, rho = 20, gamma_rule = "sqrt",
                     epsilon = 5) {
  message(sprintf(
    "Running %s: K=%d Q=%d dp=%s gamma=%s",
    label, K, Q, dp, gamma_rule
  ))
  fit <- fedspd_dp_qr(
    dat$X, dat$y,
    client_indices = client_indices,
    tau = dat$tau,
    smooth_mu = 0.05,
    rounds = rounds,
    Q = Q,
    clients_per_round = K,
    batch_size = 30,
    rho = rho,
    gamma_rule = gamma_rule,
    gamma0 = 4,
    dp = dp,
    epsilon = epsilon,
    delta = 1e-4,
    G = 1,
    phi = 1,
    d_lambda = 1,
    d_x = 1,
    beta_ref = dat$beta,
    trace_every = 5,
    seed = 20260526
  )

  trace <- fit$trace
  trace$label <- label
  trace$K <- K
  trace$Q <- Q
  trace$dp <- dp

  summary <- data.frame(
    label = label,
    K = K,
    Q = Q,
    dp = dp,
    gamma_rule = gamma_rule,
    epsilon = epsilon,
    rounds = rounds,
    final_qr_objective = fit$qr_objective,
    final_smooth_objective = fit$smooth_objective,
    final_consensus_gap = fit$consensus_gap,
    final_beta_l2_error = sqrt(sum((fit$beta - dat$beta)^2))
  )

  list(summary = summary, trace = trace, fit = fit)
}

dat <- make_qr_sim(n = 800, p = 10, tau = 0.5, noise = "asymmetric", seed = 20260526)
n_clients <- 20
client_indices <- iid_partition(nrow(dat$X), n_clients = n_clients, seed = 20260526)

central <- qr_pdhg(
  dat$X, dat$y,
  tau = dat$tau,
  max_iter = 2500,
  step_rule = "box",
  trace_every = 50,
  seed = 20260526
)

cases <- list(
  list(label = "K_full_Q5_noDP", K = 20, Q = 5, dp = FALSE),
  list(label = "K_half_Q5_noDP", K = 10, Q = 5, dp = FALSE),
  list(label = "K_quarter_Q5_noDP", K = 5, Q = 5, dp = FALSE),
  list(label = "K_half_Q1_noDP", K = 10, Q = 1, dp = FALSE),
  list(label = "K_half_Q10_noDP", K = 10, Q = 10, dp = FALSE),
  list(label = "K_half_Q5_DP", K = 10, Q = 5, dp = TRUE),
  list(label = "K_full_Q5_paper_gamma", K = 20, Q = 5, dp = FALSE,
       gamma_rule = "paper", epsilon = 100)
)

runs <- lapply(cases, function(case) {
  run_case(
    dat,
    client_indices,
    label = case$label,
    K = case$K,
    Q = case$Q,
    dp = case$dp,
    gamma_rule = if (is.null(case$gamma_rule)) "sqrt" else case$gamma_rule,
    epsilon = if (is.null(case$epsilon)) 5 else case$epsilon
  )
})

summary_tbl <- do.call(rbind, lapply(runs, `[[`, "summary"))
trace_tbl <- do.call(rbind, lapply(runs, `[[`, "trace"))

central_summary <- data.frame(
  label = "central_qr_pdhg",
  K = NA_integer_,
  Q = NA_integer_,
  dp = FALSE,
  gamma_rule = NA_character_,
  epsilon = NA_real_,
  rounds = 2500,
  final_qr_objective = central$objective,
  final_smooth_objective = smooth_qr_objective(
    dat$X, dat$y, central$beta,
    tau = dat$tau,
    mu = 0.05
  ),
  final_consensus_gap = NA_real_,
  final_beta_l2_error = sqrt(sum((central$beta - dat$beta)^2))
)
summary_tbl <- rbind(summary_tbl, central_summary)

write.csv(summary_tbl,
          file.path(root, "results", "fedspd_dp_reproduction_summary.csv"),
          row.names = FALSE)
write.csv(trace_tbl,
          file.path(root, "results", "fedspd_dp_reproduction_trace.csv"),
          row.names = FALSE)

cat("\nFedSPD-DP reproduction summary:\n")
print(summary_tbl)
