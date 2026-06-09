args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
} else {
  file.path(getwd(), "scripts", "run_box_dual_comparison.R")
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

source(file.path(root, "R", "prox.R"))
source(file.path(root, "R", "objective.R"))
source(file.path(root, "R", "simulate.R"))
source(file.path(root, "R", "qr_pdhg.R"))
source(file.path(root, "R", "fedspd_dp_qr.R"))
source(file.path(root, "R", "qr_box_fed_pdhg.R"))

dir.create(file.path(root, "results"), recursive = TRUE, showWarnings = FALSE)

run_box <- function(dat, client_indices, tau, label, K, batch_size,
                    rounds = 1500) {
  fit <- qr_box_fed_pdhg(
    dat$X, dat$y,
    client_indices = client_indices,
    tau = tau,
    rounds = rounds,
    clients_per_round = K,
    batch_size = batch_size,
    step_rule = "operator",
    aggregation = "cached",
    beta_ref = dat$beta,
    trace_every = 10,
    seed = 20260526
  )

  trace <- fit$trace
  trace$tau <- tau
  trace$method <- "qr_box_fed_pdhg"
  trace$label <- label
  trace$K <- K
  trace$batch_size <- batch_size

  summary <- data.frame(
    tau = tau,
    method = "qr_box_fed_pdhg",
    label = label,
    K = K,
    batch_size = batch_size,
    rounds = rounds,
    objective = fit$objective,
    beta_l2_error = sqrt(sum((fit$beta - dat$beta)^2)),
    dual_min = min(unlist(lapply(fit$clients, `[[`, "v"), use.names = FALSE)),
    dual_max = max(unlist(lapply(fit$clients, `[[`, "v"), use.names = FALSE))
  )

  list(summary = summary, trace = trace)
}

run_fedspd_check <- function(dat, client_indices, tau, label, K, batch_size,
                             rounds = 250) {
  fit <- fedspd_dp_qr(
    dat$X, dat$y,
    client_indices = client_indices,
    tau = tau,
    loss = "check",
    rounds = rounds,
    Q = 5,
    clients_per_round = K,
    batch_size = batch_size,
    rho = 20,
    gamma_rule = "sqrt",
    gamma0 = 4,
    dp = FALSE,
    beta_ref = dat$beta,
    trace_every = 10,
    seed = 20260526
  )

  trace <- data.frame(
    round = fit$trace$round,
    objective = fit$trace$qr_objective,
    beta_l2_error = fit$trace$beta_l2_error,
    selected_clients = fit$trace$selected_clients,
    mean_selected_n = NA_real_,
    dual_min = NA_real_,
    dual_max = NA_real_,
    tau = tau,
    method = "fedspd_check",
    label = label,
    K = K,
    batch_size = batch_size
  )

  summary <- data.frame(
    tau = tau,
    method = "fedspd_check",
    label = label,
    K = K,
    batch_size = batch_size,
    rounds = rounds,
    objective = fit$qr_objective,
    beta_l2_error = sqrt(sum((fit$beta - dat$beta)^2)),
    dual_min = NA_real_,
    dual_max = NA_real_
  )

  list(summary = summary, trace = trace)
}

taus <- c(0.5, 0.75, 0.9)
n_clients <- 20
runs <- list()
central_rows <- list()
k <- 1

for (tau in taus) {
  dat <- make_qr_sim(n = 800, p = 10, tau = tau, noise = "asymmetric", seed = 20260526)
  client_indices <- iid_partition(nrow(dat$X), n_clients = n_clients, seed = 20260526)
  full_batch <- max(lengths(client_indices))

  central <- qr_pdhg(
    dat$X, dat$y,
    tau = tau,
    max_iter = 2500,
    step_rule = "generic",
    trace_every = 50,
    seed = 20260526
  )

  central_rows[[length(central_rows) + 1]] <- data.frame(
    tau = tau,
    method = "central_qr_pdhg",
    label = "central",
    K = NA_integer_,
    batch_size = NA_integer_,
    rounds = 2500,
    objective = central$objective,
    beta_l2_error = sqrt(sum((central$beta - dat$beta)^2)),
    dual_min = NA_real_,
    dual_max = NA_real_
  )

  cases <- list(
    list(label = "full_clients_full_batch", K = n_clients, batch_size = full_batch),
    list(label = "half_clients_full_batch", K = ceiling(n_clients / 2), batch_size = full_batch),
    list(label = "full_clients_minibatch", K = n_clients, batch_size = 30),
    list(label = "half_clients_minibatch", K = ceiling(n_clients / 2), batch_size = 30)
  )

  for (case in cases) {
    message(sprintf("Running box tau=%.2f %s", tau, case$label))
    runs[[k]] <- run_box(
      dat, client_indices, tau,
      label = case$label,
      K = case$K,
      batch_size = case$batch_size
    )
    k <- k + 1
  }

  message(sprintf("Running FedSPD-check reference tau=%.2f", tau))
  runs[[k]] <- run_fedspd_check(
    dat, client_indices, tau,
    label = "full_clients_minibatch",
    K = n_clients,
    batch_size = 30
  )
  k <- k + 1
}

summary_tbl <- do.call(rbind, lapply(runs, `[[`, "summary"))
trace_tbl <- do.call(rbind, lapply(runs, `[[`, "trace"))
summary_tbl <- rbind(summary_tbl, do.call(rbind, central_rows))

write.csv(summary_tbl,
          file.path(root, "results", "box_dual_comparison_summary.csv"),
          row.names = FALSE)
write.csv(trace_tbl,
          file.path(root, "results", "box_dual_comparison_trace.csv"),
          row.names = FALSE)

cat("\nBox-dual comparison summary:\n")
print(summary_tbl)

