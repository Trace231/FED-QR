args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
} else {
  file.path(getwd(), "scripts", "run_hard_baseline_comparison_with_plots.R")
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

source(file.path(root, "R", "prox.R"))
source(file.path(root, "R", "objective.R"))
source(file.path(root, "R", "simulate.R"))
source(file.path(root, "R", "qr_pdhg.R"))
source(file.path(root, "R", "fedspd_dp_qr.R"))
source(file.path(root, "R", "qr_box_fed_pdhg.R"))
source(file.path(root, "R", "qr_admm.R"))
source(file.path(root, "R", "baselines.R"))

dir.create(file.path(root, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "figures"), recursive = TRUE, showWarnings = FALSE)

trace_from_fed_gradient <- function(fit, method, scenario, heterogeneity, seed, target_obj) {
  transform(
    fit$trace,
    method = method,
    scenario = scenario,
    heterogeneity = heterogeneity,
    seed = seed,
    gap = objective - target_obj
  )
}

run_method <- function(method, dat, scenario, K, batch_size, target_obj, seed, rounds = 600) {
  if (method == "QR box-dual") {
    fit <- qr_box_fed_pdhg(
      dat$X, dat$y,
      client_indices = dat$client_indices,
      tau = dat$tau,
      rounds = rounds,
      clients_per_round = K,
      batch_size = batch_size,
      step_rule = "operator",
      aggregation = "cached",
      beta_ref = dat$beta,
      trace_every = 5,
      seed = seed
    )
    trace <- transform(
      fit$trace,
      method = method,
      scenario = scenario,
      heterogeneity = dat$heterogeneity,
      seed = seed,
      gap = objective - target_obj
    )
    beta <- fit$beta
    obj <- fit$objective
  } else if (method == "FedSPD-check") {
    fit <- fedspd_dp_qr(
      dat$X, dat$y,
      client_indices = dat$client_indices,
      tau = dat$tau,
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
      trace_every = 5,
      seed = seed
    )
    trace <- data.frame(
      round = fit$trace$round,
      objective = fit$trace$qr_objective,
      beta_l2_error = fit$trace$beta_l2_error,
      selected_clients = fit$trace$selected_clients,
      mean_selected_n = NA_real_,
      dual_min = NA_real_,
      dual_max = NA_real_,
      method = method,
      scenario = scenario,
      heterogeneity = dat$heterogeneity,
      seed = seed,
      gap = fit$trace$qr_objective - target_obj
    )
    beta <- fit$beta
    obj <- fit$qr_objective
  } else if (method == "FedSubGrad") {
    fit <- fed_subgrad_qr(
      dat$X, dat$y,
      client_indices = dat$client_indices,
      tau = dat$tau,
      rounds = rounds,
      clients_per_round = K,
      batch_size = batch_size,
      beta_true = dat$beta,
      trace_every = 5,
      seed = seed
    )
    trace <- trace_from_fed_gradient(fit, method, scenario, dat$heterogeneity, seed, target_obj)
    beta <- fit$beta
    obj <- fit$objective
  } else if (method == "FedQR-ADMM") {
    fit <- fed_qr_admm(
      dat$X, dat$y,
      client_indices = dat$client_indices,
      tau = dat$tau,
      rounds = rounds,
      clients_per_round = K,
      batch_size = batch_size,
      rho_consensus = 1,
      rho_residual = 0.3,
      inner_iter = 25,
      beta_ref = dat$beta,
      trace_every = 5,
      seed = seed
    )
    trace <- data.frame(
      round = fit$trace$round,
      objective = fit$trace$objective,
      method = method,
      scenario = scenario,
      heterogeneity = dat$heterogeneity,
      seed = seed,
      gap = fit$trace$objective - target_obj
    )
    beta <- fit$beta
    obj <- fit$objective
  } else if (method == "FSPG-smooth") {
    fit <- fed_smooth_qr(
      dat$X, dat$y,
      client_indices = dat$client_indices,
      tau = dat$tau,
      rounds = rounds,
      clients_per_round = K,
      batch_size = batch_size,
      smooth_mu = 0.05,
      beta_true = dat$beta,
      trace_every = 5,
      seed = seed
    )
    trace <- trace_from_fed_gradient(fit, method, scenario, dat$heterogeneity, seed, target_obj)
    beta <- fit$beta
    obj <- fit$objective
  } else {
    stop("Unknown method: ", method)
  }

  summary <- data.frame(
    method = method,
    scenario = scenario,
    heterogeneity = dat$heterogeneity,
    seed = seed,
    final_objective = obj,
    final_gap = obj - target_obj,
    beta_l2_error = sqrt(sum((beta - dat$beta)^2))
  )

  list(summary = summary, trace = trace)
}

seeds <- c(20260526, 20260527, 20260528)
heterogeneity_levels <- c("mild", "hard", "extreme")
n_clients <- 20
rounds <- 600
methods <- c("QR box-dual", "FedQR-ADMM", "FedSubGrad", "FSPG-smooth", "FedSPD-check")

scenarios <- data.frame(
  scenario = c("Deterministic", "R1 + R2"),
  K = c(n_clients, 4),
  batch_size = c(100000, 10)
)

all_runs <- list()
central_rows <- list()
k <- 1
cr <- 1

for (seed in seeds) {
  for (heterogeneity in heterogeneity_levels) {
    dat <- make_hard_federated_qr_sim(
      n_clients = n_clients,
      n_per_client = 80,
      p = 15,
      tau = 0.9,
      heterogeneity = heterogeneity,
      seed = seed
    )

    central <- qr_pdhg(
      dat$X, dat$y,
      tau = dat$tau,
      max_iter = 4000,
      step_rule = "generic",
      trace_every = 100,
      seed = seed
    )
    target_obj <- central$objective
    central_rows[[cr]] <- data.frame(
      method = "Central QR-PDHG",
      scenario = "Central",
      heterogeneity = heterogeneity,
      seed = seed,
      final_objective = central$objective,
      final_gap = 0,
      beta_l2_error = sqrt(sum((central$beta - dat$beta)^2))
    )
    cr <- cr + 1

    full_batch <- max(lengths(dat$client_indices))
    for (s in seq_len(nrow(scenarios))) {
      scenario <- scenarios$scenario[s]
      K <- scenarios$K[s]
      batch_size <- if (scenarios$batch_size[s] >= 100000) full_batch else scenarios$batch_size[s]
      for (method in methods) {
        message(sprintf(
          "seed=%d hetero=%s scenario=%s method=%s",
          seed, heterogeneity, scenario, method
        ))
        all_runs[[k]] <- run_method(
          method, dat, scenario, K, batch_size, target_obj, seed, rounds
        )
        k <- k + 1
      }
    }
  }
}

summary_tbl <- do.call(rbind, lapply(all_runs, `[[`, "summary"))
trace_min <- lapply(lapply(all_runs, `[[`, "trace"), function(d) {
  d[, c("round", "objective", "method", "scenario", "heterogeneity", "seed", "gap")]
})
trace_tbl <- do.call(rbind, trace_min)
summary_tbl <- rbind(summary_tbl, do.call(rbind, central_rows))

agg_tbl <- aggregate(
  cbind(final_gap, beta_l2_error) ~ method + scenario + heterogeneity,
  data = subset(summary_tbl, method != "Central QR-PDHG"),
  FUN = mean
)

trace_agg <- aggregate(
  cbind(objective, gap) ~ method + scenario + heterogeneity + round,
  data = trace_tbl,
  FUN = mean
)

write.csv(summary_tbl,
          file.path(root, "results", "hard_baseline_summary.csv"),
          row.names = FALSE)
write.csv(agg_tbl,
          file.path(root, "results", "hard_baseline_aggregate.csv"),
          row.names = FALSE)
write.csv(trace_tbl,
          file.path(root, "results", "hard_baseline_trace.csv"),
          row.names = FALSE)
write.csv(trace_agg,
          file.path(root, "results", "hard_baseline_trace_aggregate.csv"),
          row.names = FALSE)

method_cols <- c(
  "QR box-dual" = "#2A6FBB",
  "FedQR-ADMM" = "#E6AB02",
  "FedSubGrad" = "#4DAF4A",
  "FSPG-smooth" = "#984EA3",
  "FedSPD-check" = "#D95F02"
)

plot_baseline_gap <- function() {
  png(file.path(root, "figures", "hard_baseline_final_gap.png"),
      width = 1800, height = 1100, res = 180)
  op <- par(mar = c(8, 5, 4, 1), mfrow = c(1, 3))
  on.exit(par(op), add = TRUE)
  for (h in heterogeneity_levels) {
    d <- subset(agg_tbl, heterogeneity == h & scenario == "R1 + R2")
    d <- d[match(names(method_cols), d$method), ]
    vals <- d$final_gap
    names(vals) <- d$method
    barplot(
      vals,
      las = 2,
      col = method_cols[d$method],
      ylab = "Final objective gap",
      main = paste("R1 + R2:", h),
      ylim = c(0, max(vals, na.rm = TRUE) * 1.2)
    )
  }
  dev.off()
}

plot_baseline_curves <- function() {
  png(file.path(root, "figures", "hard_baseline_r1r2_convergence.png"),
      width = 1800, height = 1100, res = 180)
  op <- par(mar = c(5, 5, 4, 1), mfrow = c(1, 3))
  on.exit(par(op), add = TRUE)
  for (h in heterogeneity_levels) {
    d <- subset(trace_agg, heterogeneity == h & scenario == "R1 + R2")
    y <- pmax(d$gap, 0)
    plot(
      NA,
      xlim = range(d$round),
      ylim = c(0, quantile(y, 0.98, na.rm = TRUE)),
      xlab = "Communication round",
      ylab = "Objective gap",
      main = paste("R1 + R2:", h)
    )
    for (m in names(method_cols)) {
      dm <- subset(d, method == m)
      lines(dm$round, pmax(dm$gap, 0), col = method_cols[[m]], lwd = 2)
    }
    legend("topright", legend = names(method_cols), col = method_cols, lwd = 2, bty = "n", cex = 0.75)
  }
  dev.off()
}

plot_baseline_log_gap <- function() {
  png(file.path(root, "figures", "hard_baseline_final_gap_log.png"),
      width = 1800, height = 1100, res = 180)
  op <- par(mar = c(8, 5, 4, 1), mfrow = c(1, 3))
  on.exit(par(op), add = TRUE)
  for (h in heterogeneity_levels) {
    d <- subset(agg_tbl, heterogeneity == h & scenario == "R1 + R2")
    d <- d[match(names(method_cols), d$method), ]
    vals <- pmax(d$final_gap, 1e-6)
    names(vals) <- d$method
    barplot(
      log10(vals),
      las = 2,
      col = method_cols[d$method],
      ylab = "log10(final objective gap)",
      main = paste("R1 + R2:", h)
    )
  }
  dev.off()
}

plot_baseline_gap()
plot_baseline_curves()
plot_baseline_log_gap()

cat("\nHard baseline aggregate:\n")
print(agg_tbl)

cat("\nFigures written to:\n")
cat(file.path(root, "figures", "hard_baseline_final_gap.png"), "\n")
cat(file.path(root, "figures", "hard_baseline_r1r2_convergence.png"), "\n")
cat(file.path(root, "figures", "hard_baseline_final_gap_log.png"), "\n")
