args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
} else {
  file.path(getwd(), "scripts", "run_advanced_innovation_experiment_with_plots.R")
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

source(file.path(root, "scripts", "load_rfedqr.R"))
make_experiment_dirs(root)

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) as.integer(value) else default
}

rounds <- env_int("ADVANCED_INNOVATION_ROUNDS", 300)
trace_every <- env_int("ADVANCED_INNOVATION_TRACE_EVERY", 25)
sim_seeds <- seq(20260609, length.out = env_int("ADVANCED_INNOVATION_SIM_SEEDS", 2))
hd_seeds <- seq(20260619, length.out = env_int("ADVANCED_INNOVATION_HD_SEEDS", 1))

methods <- c(
  "QR box-dual",
  "QR box-dual stale",
  "QR box-dual robust",
  "QR box-dual stale+robust",
  "QR box-dual adaptive",
  "QR box-dual adaptive+VR",
  "FSPG-smooth",
  "FedSPD-check"
)

method_controls <- list(
  "QR box-dual" = list(step_rule = "box"),
  "QR box-dual stale" = list(step_rule = "box", staleness_rate = 0.05, staleness_floor = 0.35),
  "QR box-dual robust" = list(step_rule = "box", client_weighting = "uniform"),
  "QR box-dual stale+robust" = list(
    step_rule = "box",
    staleness_rate = 0.05,
    staleness_floor = 0.35,
    client_weighting = "uniform"
  ),
  "QR box-dual adaptive" = list(
    step_rule = "box",
    staleness_rate = 0.05,
    staleness_floor = 0.35,
    adaptive_client_power = 1,
    adaptive_client_blend = 0.85,
    adaptive_client_floor = 0.05,
    adaptive_client_smooth = 0.5
  ),
  "QR box-dual adaptive+VR" = list(
    step_rule = "box",
    staleness_rate = 0.05,
    staleness_floor = 0.35,
    adaptive_client_power = 1,
    adaptive_client_blend = 0.85,
    adaptive_client_floor = 0.05,
    adaptive_client_smooth = 0.5,
    variance_reduction = "control_variate",
    vr_alpha = 0.05,
    vr_blend = 0.1,
    vr_max_correction_ratio = 0.5
  ),
  "FedSPD-check" = list(Q = 3, gamma0 = 4)
)

target_for <- function(X, y, tau, seed) {
  fit <- qr_pdhg(
    X, y,
    tau = tau,
    max_iter = 5000,
    step_rule = "box",
    trace_every = 100,
    seed = seed
  )
  list(beta = fit$beta, objective = fit$objective)
}

collect_run <- function(dataset, setting, X, y, client_indices, tau, seed,
                        target, clients_per_round, batch_size, rounds) {
  result <- run_fedqr_methods(
    methods,
    X,
    y,
    client_indices = client_indices,
    tau = tau,
    rounds = rounds,
    clients_per_round = clients_per_round,
    batch_size = batch_size,
    seed = seed,
    trace_every = trace_every,
    method_controls = method_controls,
    term_names = colnames(X),
    verbose = TRUE
  )

  summary <- result$summary
  summary$dataset <- dataset
  summary$setting <- setting
  summary$tau <- tau
  summary$seed <- seed
  summary$target_objective <- target$objective
  summary$target_gap <- summary$objective - target$objective
  summary <- summary[, c(
    "dataset", "setting", "tau", "seed", "method",
    "objective", "target_objective", "target_gap",
    "gap_to_best_observed", "beta_norm"
  )]

  trace <- result$trace
  trace$dataset <- dataset
  trace$setting <- setting
  trace$tau <- tau
  trace$seed <- seed
  trace$target_gap <- trace$objective - target$objective

  client_loss <- do.call(rbind, lapply(result$fits, function(fit) {
    out <- client_loss_summary(X, y, fit$beta, client_indices, tau = tau)
    out$dataset <- dataset
    out$setting <- setting
    out$tau <- tau
    out$seed <- seed
    out$method <- fit$method
    out
  }))
  client_loss <- client_loss[, c(
    "dataset", "setting", "tau", "seed", "method",
    "global_mean_loss", "client_mean_loss", "worst_client_loss",
    "client_loss_sd", "client_q90_loss", "min_client_loss"
  )]

  calibration <- do.call(rbind, lapply(result$fits, function(fit) {
    raw <- calibration_summary(X, y, fit$beta, tau = tau, client_indices = client_indices)
    raw$calibration <- "raw"

    global_cal <- calibrate_quantile_intercept(X, y, fit$beta, tau = tau, mode = "global")
    global <- calibration_summary(X, y, global_cal$beta, tau = tau, client_indices = client_indices)
    global$calibration <- "global_intercept"

    client_cal <- calibrate_quantile_intercept(
      X, y, fit$beta,
      tau = tau,
      client_indices = client_indices,
      mode = "client_offset"
    )
    client <- calibration_summary(
      X, y, client_cal$beta,
      tau = tau,
      client_indices = client_indices,
      offsets = client_cal$offsets
    )
    client$calibration <- "client_offset"

    out <- rbind(raw, global, client)
    out$dataset <- dataset
    out$setting <- setting
    out$tau <- tau
    out$seed <- seed
    out$method <- fit$method
    out
  }))
  calibration <- calibration[, c(
    "dataset", "setting", "tau", "seed", "method", "calibration",
    "global_coverage", "global_coverage_error",
    "mean_client_coverage_error", "worst_client_coverage_error",
    "coverage_sd", "min_client_coverage", "max_client_coverage"
  )]

  list(summary = summary, trace = trace, client_loss = client_loss, calibration = calibration)
}

runs <- list()
k <- 1

for (heterogeneity in c("hard", "extreme")) {
  for (tau in c(0.9, 0.95)) {
    for (seed in sim_seeds) {
      set.seed(seed + 1000)
      n_per_client <- sample(seq(40, 140, by = 10), 20, replace = TRUE)
      dat <- make_hard_federated_qr_sim(
        n_clients = 20,
        n_per_client = n_per_client,
        p = 15,
        tau = tau,
        heterogeneity = heterogeneity,
        seed = seed
      )
      target <- target_for(dat$X, dat$y, tau, seed)
      runs[[k]] <- collect_run(
        dataset = "simulation",
        setting = heterogeneity,
        X = dat$X,
        y = dat$y,
        client_indices = dat$client_indices,
        tau = tau,
        seed = seed,
        target = target,
        clients_per_round = 4,
        batch_size = 20,
        rounds = rounds
      )
      k <- k + 1
    }
  }
}

hd_path <- file.path(root, "data", "processed", "heart_disease_uci_4center.csv")
if (file.exists(hd_path)) {
  hd <- read.csv(hd_path)
  vars <- c(
    "center", "thalach", "age", "sex", "cp", "trestbps", "chol",
    "fbs", "restecg", "exang", "oldpeak", "num"
  )
  hd_model <- hd[, vars]
  hd_model <- hd_model[complete.cases(hd_model), ]
  for (v in c("sex", "cp", "fbs", "restecg", "exang", "num")) {
    hd_model[[v]] <- factor(hd_model[[v]])
  }
  for (v in c("age", "trestbps", "chol", "oldpeak")) {
    hd_model[[paste0(v, "_z")]] <- as.numeric(scale(hd_model[[v]]))
  }
  hd_model$thalach_z <- as.numeric(scale(hd_model$thalach))
  formula_hd <- thalach_z ~ age_z + sex + cp + trestbps_z + chol_z +
    fbs + restecg + exang + oldpeak_z + num
  X_hd <- model.matrix(formula_hd, data = hd_model)
  y_hd <- hd_model$thalach_z
  clients_hd <- split(seq_along(y_hd), hd_model$center)
  clients_hd <- clients_hd[order(names(clients_hd))]
  full_batch <- max(lengths(clients_hd))

  for (tau in c(0.5, 0.75, 0.9)) {
    target <- target_for(X_hd, y_hd, tau, hd_seeds[1])
    for (seed in hd_seeds) {
      runs[[k]] <- collect_run(
        dataset = "heart_disease",
        setting = "four_center_r1r2",
        X = X_hd,
        y = y_hd,
        client_indices = clients_hd,
        tau = tau,
        seed = seed,
        target = target,
        clients_per_round = 2,
        batch_size = min(20, full_batch),
        rounds = rounds
      )
      k <- k + 1
    }
  }
}

summary_tbl <- do.call(rbind, lapply(runs, `[[`, "summary"))
trace_tbl <- do.call(rbind, lapply(runs, `[[`, "trace"))
client_loss_tbl <- do.call(rbind, lapply(runs, `[[`, "client_loss"))
calibration_tbl <- do.call(rbind, lapply(runs, `[[`, "calibration"))

write.csv(summary_tbl,
          file.path(root, "results", "advanced_innovation_summary.csv"),
          row.names = FALSE)
write.csv(trace_tbl,
          file.path(root, "results", "advanced_innovation_trace.csv"),
          row.names = FALSE)
write.csv(client_loss_tbl,
          file.path(root, "results", "advanced_innovation_client_loss.csv"),
          row.names = FALSE)
write.csv(calibration_tbl,
          file.path(root, "results", "advanced_innovation_calibration.csv"),
          row.names = FALSE)

method_cols <- fedqr_default_method_cols(methods)

agg_gap <- aggregate(target_gap ~ dataset + setting + tau + method,
                     data = summary_tbl, FUN = mean)
png(file.path(root, "figures", "advanced_innovation_final_gap.png"),
    width = 2200, height = 1300, res = 180)
op <- par(mar = c(7, 5, 4, 1), mfrow = c(2, 4))
for (key in unique(paste(agg_gap$dataset, agg_gap$setting, agg_gap$tau, sep = " | "))) {
  d <- agg_gap[paste(agg_gap$dataset, agg_gap$setting, agg_gap$tau, sep = " | ") == key, ]
  d <- d[match(methods, d$method), ]
  vals <- pmax(d$target_gap, 1e-8)
  names(vals) <- d$method
  barplot(log10(vals), las = 2, cex.names = 0.65, col = method_cols[d$method],
          ylab = "log10 gap to central QR", main = key)
}
par(op)
dev.off()

agg_fair <- aggregate(worst_client_loss ~ dataset + setting + tau + method,
                      data = client_loss_tbl, FUN = mean)
png(file.path(root, "figures", "advanced_innovation_client_fairness.png"),
    width = 2200, height = 1300, res = 180)
op <- par(mar = c(7, 5, 4, 1), mfrow = c(2, 4))
for (key in unique(paste(agg_fair$dataset, agg_fair$setting, agg_fair$tau, sep = " | "))) {
  d <- agg_fair[paste(agg_fair$dataset, agg_fair$setting, agg_fair$tau, sep = " | ") == key, ]
  d <- d[match(methods, d$method), ]
  vals <- d$worst_client_loss
  names(vals) <- d$method
  barplot(vals, las = 2, cex.names = 0.65, col = method_cols[d$method],
          ylab = "Worst-client check loss", main = key)
}
par(op)
dev.off()

cal_plot <- calibration_tbl[calibration_tbl$calibration %in% c("raw", "global_intercept"), ]
agg_cal <- aggregate(global_coverage_error ~ dataset + setting + tau + method + calibration,
                     data = cal_plot, FUN = mean)
png(file.path(root, "figures", "advanced_innovation_calibration.png"),
    width = 2200, height = 1300, res = 180)
op <- par(mar = c(7, 5, 4, 1), mfrow = c(2, 4))
for (key in unique(paste(agg_cal$dataset, agg_cal$setting, agg_cal$tau, sep = " | "))) {
  d <- agg_cal[paste(agg_cal$dataset, agg_cal$setting, agg_cal$tau, sep = " | ") == key, ]
  raw <- d[d$calibration == "raw", ]
  global <- d[d$calibration == "global_intercept", ]
  raw <- raw[match(methods, raw$method), ]
  global <- global[match(methods, global$method), ]
  vals <- rbind(raw = raw$global_coverage_error, calibrated = global$global_coverage_error)
  colnames(vals) <- methods
  barplot(vals, beside = TRUE, las = 2, cex.names = 0.65,
          col = c("#999999", "#E69F00"), ylab = "Global coverage error", main = key)
  legend("topright", legend = rownames(vals), fill = c("#999999", "#E69F00"), bty = "n", cex = 0.7)
}
par(op)
dev.off()

stale_trace <- trace_tbl[trace_tbl$method %in% c("QR box-dual stale", "QR box-dual stale+robust"), ]
png(file.path(root, "figures", "advanced_innovation_staleness_trace.png"),
    width = 1600, height = 1000, res = 180)
op <- par(mar = c(5, 5, 4, 1))
plot(NA, xlim = range(stale_trace$round), ylim = range(stale_trace$mean_stale_weight, na.rm = TRUE),
     xlab = "Communication round", ylab = "Mean raw stale weight",
     main = "Staleness-aware cached direction weights")
for (method in unique(stale_trace$method)) {
  d <- aggregate(mean_stale_weight ~ round + method,
                 data = stale_trace[stale_trace$method == method, ], FUN = mean)
  lines(d$round, d$mean_stale_weight, col = method_cols[[method]], lwd = 2)
}
legend("bottomleft", legend = unique(stale_trace$method),
       col = method_cols[unique(stale_trace$method)], lwd = 2, bty = "n")
par(op)
dev.off()

cat("\nAdvanced innovation summary aggregate:\n")
print(agg_gap[order(agg_gap$dataset, agg_gap$setting, agg_gap$tau, agg_gap$target_gap), ])
cat("\nClient fairness aggregate:\n")
print(agg_fair[order(agg_fair$dataset, agg_fair$setting, agg_fair$tau, agg_fair$worst_client_loss), ])
cat("\nCalibration aggregate:\n")
print(agg_cal[order(agg_cal$dataset, agg_cal$setting, agg_cal$tau, agg_cal$global_coverage_error), ])
