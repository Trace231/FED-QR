args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
} else {
  file.path(getwd(), "scripts", "run_admm_validation.R")
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

source(file.path(root, "R", "prox.R"))
source(file.path(root, "R", "objective.R"))
source(file.path(root, "R", "simulate.R"))
source(file.path(root, "R", "qr_pdhg.R"))
source(file.path(root, "R", "qr_admm.R"))

dir.create(file.path(root, "results"), recursive = TRUE, showWarnings = FALSE)

taus <- c(0.5, 0.75, 0.9)
runs <- lapply(taus, function(tau) {
  dat <- make_qr_sim(n = 500, p = 8, tau = tau, noise = "asymmetric", seed = 20260526)

  pdhg <- qr_pdhg(
    dat$X, dat$y,
    tau = tau,
    max_iter = 5000,
    step_rule = "generic",
    trace_every = 100,
    seed = 20260526
  )

  admm <- qr_admm(
    dat$X, dat$y,
    tau = tau,
    rho = 0.3,
    max_iter = 10000,
    trace_every = 100,
    check_convergence = FALSE
  )

  fed_admm <- fed_qr_admm(
    dat$X, dat$y,
    client_indices = iid_partition(nrow(dat$X), n_clients = 10, seed = 20260526),
    tau = tau,
    rounds = 500,
    clients_per_round = 10,
    batch_size = 50,
    rho_consensus = 1,
    rho_residual = 0.3,
    inner_iter = 30,
    beta_ref = dat$beta,
    trace_every = 25,
    seed = 20260526
  )

  data.frame(
    tau = tau,
    method = c("central_qr_pdhg", "central_qr_admm", "fed_qr_admm"),
    objective = c(pdhg$objective, admm$objective, fed_admm$objective),
    objective_gap_vs_pdhg = c(0, admm$objective - pdhg$objective, fed_admm$objective - pdhg$objective),
    beta_l2_error = c(
      sqrt(sum((pdhg$beta - dat$beta)^2)),
      sqrt(sum((admm$beta - dat$beta)^2)),
      sqrt(sum((fed_admm$beta - dat$beta)^2))
    )
  )
})

summary_tbl <- do.call(rbind, runs)
write.csv(summary_tbl,
          file.path(root, "results", "admm_validation_summary.csv"),
          row.names = FALSE)

cat("\nADMM validation summary:\n")
print(summary_tbl)

