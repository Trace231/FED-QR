args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
} else {
  file.path(getwd(), "scripts", "run_hd_baseline_comparison_with_plots.R")
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

source(file.path(root, "scripts", "load_rfedqr.R"))
make_experiment_dirs(root)

data_path <- file.path(root, "data", "processed", "heart_disease_uci_4center.csv")
if (!file.exists(data_path)) {
  stop("Heart Disease data not found. Run scripts/inspect_heart_disease.R first.")
}

hd <- read.csv(data_path)

vars <- c(
  "center", "thalach", "age", "sex", "cp", "trestbps", "chol",
  "fbs", "restecg", "exang", "oldpeak", "num"
)
hd_model <- hd[, vars]
hd_model <- hd_model[complete.cases(hd_model), ]

factor_vars <- c("sex", "cp", "fbs", "restecg", "exang", "num")
for (v in factor_vars) {
  hd_model[[v]] <- factor(hd_model[[v]])
}

numeric_vars <- c("age", "trestbps", "chol", "oldpeak")
for (v in numeric_vars) {
  hd_model[[paste0(v, "_z")]] <- as.numeric(scale(hd_model[[v]]))
}
hd_model$thalach_z <- as.numeric(scale(hd_model$thalach))

formula_hd <- thalach_z ~ age_z + sex + cp + trestbps_z + chol_z +
  fbs + restecg + exang + oldpeak_z + num

X <- model.matrix(formula_hd, data = hd_model)
y <- hd_model$thalach_z
center <- hd_model$center
client_indices <- split(seq_along(y), center)
client_indices <- client_indices[order(names(client_indices))]
n_clients <- length(client_indices)
full_batch <- max(lengths(client_indices))
mini_batch <- 20

response_summary <- do.call(rbind, lapply(names(client_indices), function(center_name) {
  d <- hd_model[client_indices[[center_name]], ]
  data.frame(
    center = center_name,
    n = nrow(d),
    thalach_mean = mean(d$thalach),
    thalach_sd = sd(d$thalach),
    thalach_q10 = as.numeric(quantile(d$thalach, 0.10)),
    thalach_q25 = as.numeric(quantile(d$thalach, 0.25)),
    thalach_q50 = as.numeric(quantile(d$thalach, 0.50)),
    thalach_q75 = as.numeric(quantile(d$thalach, 0.75)),
    thalach_q90 = as.numeric(quantile(d$thalach, 0.90))
  )
}))

write.csv(response_summary,
          file.path(root, "results", "hd_center_response_summary.csv"),
          row.names = FALSE)

target_for_tau <- function(tau, seed) {
  central <- qr_pdhg(
    X, y,
    tau = tau,
    max_iter = 5000,
    step_rule = "box",
    trace_every = 100,
    seed = seed
  )

  target_beta <- central$beta
  target_obj <- central$objective
  target_method <- "Central QR-PDHG"

  if (requireNamespace("quantreg", quietly = TRUE)) {
    rq_fit <- quantreg::rq(formula_hd, data = hd_model, tau = tau)
    rq_beta <- as.numeric(coef(rq_fit))
    rq_obj <- qr_objective(X, y, rq_beta, tau = tau)
    if (rq_obj <= target_obj + 1e-7) {
      target_beta <- rq_beta
      target_obj <- rq_obj
      target_method <- "quantreg::rq"
    }
  }

  list(
    beta = target_beta,
    objective = target_obj,
    method = target_method,
    central_objective = central$objective,
    central_beta = central$beta
  )
}

trace_standard <- function(trace, objective_col, method, scenario, tau, seed, target_obj) {
  data.frame(
    tau = tau,
    method = method,
    scenario = scenario,
    seed = seed,
    round = trace$round,
    objective = trace[[objective_col]],
    gap = trace[[objective_col]] - target_obj
  )
}

run_method <- function(method, tau, scenario, K, batch_size, target, seed, rounds = 600) {
  control <- switch(
    match_fedqr_method(method),
    "QR box-dual" = list(step_rule = "box", aggregation = "cached", beta_ref = target$beta),
    "FedQR-ADMM" = list(rho_consensus = 1, rho_residual = 0.3, inner_iter = 25, beta_ref = target$beta),
    "FedSubGrad" = list(beta_true = target$beta),
    "FSPG-smooth" = list(smooth_mu = 0.05, beta_true = target$beta),
    "FedSPD-check" = list(Q = 5, rho = 20, gamma_rule = "sqrt", gamma0 = 4, dp = FALSE, beta_ref = target$beta)
  )
  fit <- fit_fedqr(
    method = method,
    X = X,
    y = y,
    client_indices = client_indices,
    tau = tau,
    rounds = rounds,
    clients_per_round = K,
    batch_size = batch_size,
    trace_every = 5,
    seed = seed,
    control = control
  )
  trace <- data.frame(
    tau = tau,
    method = fit$method,
    scenario = scenario,
    seed = seed,
    round = fit$trace$round,
    objective = fit$trace$objective,
    gap = fit$trace$objective - target$objective
  )
  beta <- fit$beta
  obj <- fit$objective

  summary <- data.frame(
    tau = tau,
    method = fit$method,
    scenario = scenario,
    seed = seed,
    final_objective = obj,
    final_gap = obj - target$objective,
    beta_l2_to_target = sqrt(sum((beta - target$beta)^2)),
    target_method = target$method,
    target_objective = target$objective
  )

  list(summary = summary, trace = trace)
}

taus <- c(0.5, 0.75, 0.9)
seeds <- c(20260526, 20260527, 20260528)
rounds <- 600
methods <- c("QR box-dual", "FedQR-ADMM", "FedSubGrad", "FSPG-smooth", "FedSPD-check")

scenarios <- data.frame(
  scenario = c("Full", "R1 client", "R2 user", "R1 + R2"),
  K = c(n_clients, 2, n_clients, 2),
  batch_size = c(full_batch, full_batch, mini_batch, mini_batch)
)

targets <- setNames(lapply(taus, target_for_tau, seed = seeds[1]), as.character(taus))
target_tbl <- do.call(rbind, lapply(names(targets), function(tau_name) {
  target <- targets[[tau_name]]
  data.frame(
    tau = as.numeric(tau_name),
    target_method = target$method,
    target_objective = target$objective,
    central_qr_pdhg_objective = target$central_objective
  )
}))
write.csv(target_tbl,
          file.path(root, "results", "hd_baseline_targets.csv"),
          row.names = FALSE)

all_runs <- list()
k <- 1
for (tau in taus) {
  target <- targets[[as.character(tau)]]
  for (s in seq_len(nrow(scenarios))) {
    scenario <- scenarios$scenario[s]
    K <- scenarios$K[s]
    batch_size <- scenarios$batch_size[s]
    scenario_seeds <- if (scenario == "Full") seeds[1] else seeds

    for (seed in scenario_seeds) {
      for (method in methods) {
        message(sprintf(
          "tau=%.2f scenario=%s seed=%d method=%s",
          tau, scenario, seed, method
        ))
        all_runs[[k]] <- run_method(
          method = method,
          tau = tau,
          scenario = scenario,
          K = K,
          batch_size = batch_size,
          target = target,
          seed = seed,
          rounds = rounds
        )
        k <- k + 1
      }
    }
  }
}

summary_tbl <- do.call(rbind, lapply(all_runs, `[[`, "summary"))
trace_tbl <- do.call(rbind, lapply(all_runs, `[[`, "trace"))

agg_tbl <- aggregate(
  cbind(final_gap, beta_l2_to_target) ~ tau + method + scenario,
  data = summary_tbl,
  FUN = mean
)
agg_sd <- aggregate(
  cbind(final_gap, beta_l2_to_target) ~ tau + method + scenario,
  data = summary_tbl,
  FUN = sd
)
names(agg_sd)[names(agg_sd) == "final_gap"] <- "final_gap_sd"
names(agg_sd)[names(agg_sd) == "beta_l2_to_target"] <- "beta_l2_to_target_sd"
agg_tbl <- merge(agg_tbl, agg_sd, by = c("tau", "method", "scenario"), all.x = TRUE)
agg_tbl$final_gap_sd[is.na(agg_tbl$final_gap_sd)] <- 0
agg_tbl$beta_l2_to_target_sd[is.na(agg_tbl$beta_l2_to_target_sd)] <- 0

trace_agg <- aggregate(
  cbind(objective, gap) ~ tau + method + scenario + round,
  data = trace_tbl,
  FUN = mean
)

write.csv(summary_tbl,
          file.path(root, "results", "hd_baseline_summary.csv"),
          row.names = FALSE)
write.csv(agg_tbl,
          file.path(root, "results", "hd_baseline_aggregate.csv"),
          row.names = FALSE)
write.csv(trace_tbl,
          file.path(root, "results", "hd_baseline_trace.csv"),
          row.names = FALSE)
write.csv(trace_agg,
          file.path(root, "results", "hd_baseline_trace_aggregate.csv"),
          row.names = FALSE)

method_cols <- c(
  "QR box-dual" = "#2A6FBB",
  "FedQR-ADMM" = "#E6AB02",
  "FedSubGrad" = "#4DAF4A",
  "FSPG-smooth" = "#984EA3",
  "FedSPD-check" = "#D95F02"
)

scenario_levels <- scenarios$scenario

plot_center_distribution <- function() {
  png(file.path(root, "figures", "hd_center_response_distribution.png"),
      width = 1500, height = 1000, res = 180)
  op <- par(mar = c(6, 5, 4, 1))
  on.exit(par(op), add = TRUE)
  boxplot(
    thalach ~ center,
    data = hd_model,
    col = "#BFD7EA",
    border = "#30475E",
    ylab = "Maximum heart rate (thalach)",
    xlab = "",
    main = "Heart Disease response distribution by center",
    las = 2
  )
  mtext(
    paste0(names(client_indices), " n=", lengths(client_indices), collapse = "    "),
    side = 3,
    line = 0.2,
    cex = 0.75
  )
  dev.off()
}

plot_gap_grid <- function() {
  png(file.path(root, "figures", "hd_baseline_final_gap_log.png"),
      width = 2200, height = 1500, res = 180)
  op <- par(mar = c(6.5, 4.5, 3.2, 0.5), mfrow = c(length(taus), length(scenario_levels)))
  on.exit(par(op), add = TRUE)
  for (tau in taus) {
    for (scenario in scenario_levels) {
      d <- agg_tbl[agg_tbl$tau == tau & agg_tbl$scenario == scenario, ]
      d <- d[match(names(method_cols), d$method), ]
      vals <- pmax(d$final_gap, 1e-8)
      names(vals) <- d$method
      barplot(
        log10(vals),
        las = 2,
        col = method_cols[d$method],
        ylab = "log10(gap)",
        main = sprintf("tau=%.2f | %s", tau, scenario)
      )
      abline(h = -4, col = "#666666", lty = 3)
    }
  }
  dev.off()
}

plot_r1r2_curves <- function() {
  png(file.path(root, "figures", "hd_baseline_r1r2_convergence.png"),
      width = 1800, height = 1050, res = 180)
  op <- par(mar = c(5, 5, 4, 1), mfrow = c(1, length(taus)))
  on.exit(par(op), add = TRUE)
  for (tau in taus) {
    d <- trace_agg[trace_agg$tau == tau & trace_agg$scenario == "R1 + R2", ]
    y_lim <- quantile(pmax(d$gap, 0), 0.98, na.rm = TRUE)
    y_lim <- max(y_lim, 1e-4)
    plot(
      NA,
      xlim = range(d$round),
      ylim = c(0, y_lim),
      xlab = "Communication round",
      ylab = "Objective gap",
      main = sprintf("R1 + R2, tau=%.2f", tau)
    )
    for (m in names(method_cols)) {
      dm <- subset(d, method == m)
      lines(dm$round, pmax(dm$gap, 0), col = method_cols[[m]], lwd = 2)
    }
    legend("topright", legend = names(method_cols), col = method_cols, lwd = 2, bty = "n", cex = 0.75)
  }
  dev.off()
}

plot_r1r2_final_gap <- function() {
  png(file.path(root, "figures", "hd_baseline_r1r2_final_gap.png"),
      width = 1800, height = 1050, res = 180)
  op <- par(mar = c(7, 5, 4, 1), mfrow = c(1, length(taus)))
  on.exit(par(op), add = TRUE)
  for (tau in taus) {
    d <- agg_tbl[agg_tbl$tau == tau & agg_tbl$scenario == "R1 + R2", ]
    d <- d[match(names(method_cols), d$method), ]
    vals <- pmax(d$final_gap, 0)
    names(vals) <- d$method
    barplot(
      vals,
      las = 2,
      col = method_cols[d$method],
      ylab = "Final objective gap",
      main = sprintf("R1 + R2, tau=%.2f", tau),
      ylim = c(0, max(vals, na.rm = TRUE) * 1.2)
    )
  }
  dev.off()
}

plot_center_distribution()
plot_gap_grid()
plot_r1r2_curves()
plot_r1r2_final_gap()

cat("\nHeart Disease client sizes:\n")
print(data.frame(center = names(client_indices), n = lengths(client_indices)))

cat("\nTargets:\n")
print(target_tbl)

cat("\nHD baseline aggregate:\n")
print(agg_tbl[order(agg_tbl$tau, agg_tbl$scenario, agg_tbl$final_gap), ])

cat("\nFigures written to:\n")
cat(file.path(root, "figures", "hd_center_response_distribution.png"), "\n")
cat(file.path(root, "figures", "hd_baseline_final_gap_log.png"), "\n")
cat(file.path(root, "figures", "hd_baseline_r1r2_convergence.png"), "\n")
cat(file.path(root, "figures", "hd_baseline_r1r2_final_gap.png"), "\n")
