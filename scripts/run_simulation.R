args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
} else {
  file.path(getwd(), "scripts", "run_simulation.R")
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
source(file.path(root, "R", "prox.R"))
source(file.path(root, "R", "objective.R"))
source(file.path(root, "R", "simulate.R"))
source(file.path(root, "R", "qr_pdhg.R"))

dir.create(file.path(root, "results"), recursive = TRUE, showWarnings = FALSE)

fit_one <- function(tau, step_rule, lambda = 0, penalty = "none") {
  dat <- make_qr_sim(n = 600, p = 12, tau = tau, noise = "asymmetric", seed = 20260526)
  fit <- qr_pdhg(
    dat$X, dat$y,
    tau = tau,
    lambda = lambda,
    penalty = penalty,
    max_iter = 1500,
    step_rule = step_rule,
    trace_every = 25,
    seed = 20260526
  )

  beta_err <- sqrt(sum((fit$beta - dat$beta)^2))
  out <- data.frame(
    tau = tau,
    step_rule = step_rule,
    lambda = lambda,
    penalty = penalty,
    objective = fit$objective,
    beta_l2_error = beta_err,
    eta = fit$eta,
    sigma = fit$sigma
  )

  trace <- fit$trace
  trace$tau <- tau
  trace$step_rule <- step_rule
  trace

  list(summary = out, trace = trace, fit = fit, data = dat)
}

taus <- c(0.5, 0.75, 0.9)
rules <- c("generic", "box", "tau_adaptive")
runs <- list()
k <- 1
for (tau in taus) {
  for (rule in rules) {
    runs[[k]] <- fit_one(tau, rule)
    k <- k + 1
  }
}

summary_tbl <- do.call(rbind, lapply(runs, `[[`, "summary"))
trace_tbl <- do.call(rbind, lapply(runs, `[[`, "trace"))

write.csv(summary_tbl, file.path(root, "results", "simulation_summary.csv"), row.names = FALSE)
write.csv(trace_tbl, file.path(root, "results", "simulation_trace.csv"), row.names = FALSE)

print(summary_tbl)

if (requireNamespace("quantreg", quietly = TRUE)) {
  dat <- make_qr_sim(n = 600, p = 12, tau = 0.5, noise = "asymmetric", seed = 20260526)
  df <- data.frame(y = dat$y, dat$X[, -1, drop = FALSE])
  rq_fit <- quantreg::rq(y ~ ., data = df, tau = 0.5)
  pdhg_fit <- qr_pdhg(dat$X, dat$y, tau = 0.5, max_iter = 2500, step_rule = "generic")
  cmp <- data.frame(
    term = names(coef(rq_fit)),
    quantreg = as.numeric(coef(rq_fit)),
    qr_pdhg = as.numeric(pdhg_fit$beta)
  )
  cmp$abs_diff <- abs(cmp$quantreg - cmp$qr_pdhg)
  write.csv(cmp, file.path(root, "results", "quantreg_comparison.csv"), row.names = FALSE)
  print(cmp)
} else {
  message("Package 'quantreg' is not installed; skipped quantreg comparison.")
}
