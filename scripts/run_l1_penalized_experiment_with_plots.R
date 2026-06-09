args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
} else {
  file.path(getwd(), "scripts", "run_l1_penalized_experiment_with_plots.R")
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

support_metrics <- function(beta_hat, beta_true, tol = 1e-3) {
  keep <- if (!is.null(names(beta_hat)) && names(beta_hat)[1] == "X.Intercept.") {
    seq_along(beta_hat)[-1]
  } else if (!is.null(names(beta_true)) && names(beta_true)[1] == "X.Intercept.") {
    seq_along(beta_true)[-1]
  } else {
    seq_along(beta_true)[-1]
  }

  true_support <- abs(beta_true[keep]) > 1e-10
  selected <- abs(beta_hat[keep]) > tol
  tp <- sum(selected & true_support)
  fp <- sum(selected & !true_support)
  fn <- sum(!selected & true_support)
  data.frame(
    selected_size = sum(selected),
    true_size = sum(true_support),
    tp = tp,
    fp = fp,
    fn = fn,
    tpr = if (sum(true_support) == 0) NA_real_ else tp / sum(true_support),
    fdr = if (sum(selected) == 0) 0 else fp / sum(selected),
    l2_error = sqrt(sum((beta_hat - beta_true)^2))
  )
}

run_method <- function(method, dat, lambda, target_obj, seed, rounds = 800) {
  n_clients <- length(dat$client_indices)
  K <- 4
  batch_size <- 10

  if (method == "QR box-dual") {
    fit <- qr_box_fed_pdhg(
      dat$X, dat$y,
      client_indices = dat$client_indices,
      tau = dat$tau,
      lambda = lambda,
      penalty = "l1",
      rounds = rounds,
      clients_per_round = K,
      batch_size = batch_size,
      step_rule = "operator",
      aggregation = "cached",
      beta_ref = dat$beta,
      trace_every = 5,
      seed = seed
    )
    beta <- fit$beta
    obj <- fit$objective
    trace <- transform(
      fit$trace,
      method = method,
      lambda = lambda,
      seed = seed,
      gap = objective - target_obj
    )
  } else if (method == "FSPG-smooth") {
    fit <- fed_smooth_qr(
      dat$X, dat$y,
      client_indices = dat$client_indices,
      tau = dat$tau,
      lambda = lambda,
      penalty = "l1",
      rounds = rounds,
      clients_per_round = K,
      batch_size = batch_size,
      smooth_mu = 0.05,
      beta_true = dat$beta,
      trace_every = 5,
      seed = seed
    )
    beta <- fit$beta
    obj <- fit$objective
    trace <- transform(
      fit$trace,
      method = method,
      lambda = lambda,
      seed = seed,
      gap = objective - target_obj
    )
  } else if (method == "FedSubGrad") {
    fit <- fed_subgrad_qr(
      dat$X, dat$y,
      client_indices = dat$client_indices,
      tau = dat$tau,
      lambda = lambda,
      penalty = "l1",
      rounds = rounds,
      clients_per_round = K,
      batch_size = batch_size,
      beta_true = dat$beta,
      trace_every = 5,
      seed = seed
    )
    beta <- fit$beta
    obj <- fit$objective
    trace <- transform(
      fit$trace,
      method = method,
      lambda = lambda,
      seed = seed,
      gap = objective - target_obj
    )
  } else if (method == "FedQR-ADMM") {
    fit <- fed_qr_admm(
      dat$X, dat$y,
      client_indices = dat$client_indices,
      tau = dat$tau,
      lambda = lambda,
      penalty = "l1",
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
    beta <- fit$beta
    obj <- fit$objective
    trace <- data.frame(
      round = fit$trace$round,
      objective = fit$trace$objective,
      beta_l2_error = fit$trace$beta_l2_error,
      method = method,
      lambda = lambda,
      seed = seed,
      gap = fit$trace$objective - target_obj
    )
  } else {
    stop("Unknown method: ", method)
  }

  supp <- support_metrics(beta, dat$beta, tol = 1e-3)
  summary <- cbind(
    data.frame(
      method = method,
      lambda = lambda,
      seed = seed,
      final_objective = obj,
      final_gap = obj - target_obj
    ),
    supp
  )

  list(summary = summary, trace = trace)
}

seeds <- c(20260526, 20260527, 20260528)
methods <- c("QR box-dual", "FSPG-smooth", "FedSubGrad", "FedQR-ADMM")
lambdas <- c(0.001, 0.005, 0.01, 0.02)

all_runs <- list()
central_rows <- list()
k <- 1
cr <- 1

for (seed in seeds) {
  dat <- make_hard_federated_qr_sim(
    n_clients = 20,
    n_per_client = 70,
    p = 60,
    sparsity = 8,
    tau = 0.9,
    heterogeneity = "hard",
    seed = seed
  )

  for (lambda in lambdas) {
    central <- qr_box_fed_pdhg(
      dat$X, dat$y,
      client_indices = list(seq_len(nrow(dat$X))),
      tau = dat$tau,
      lambda = lambda,
      penalty = "l1",
      rounds = 3000,
      clients_per_round = 1,
      batch_size = nrow(dat$X),
      step_rule = "operator",
      aggregation = "cached",
      beta_ref = dat$beta,
      trace_every = 100,
      seed = seed
    )
    target_obj <- central$objective
    central_rows[[cr]] <- cbind(
      data.frame(
        method = "Central L1 QR box-dual",
        lambda = lambda,
        seed = seed,
        final_objective = central$objective,
        final_gap = 0
      ),
      support_metrics(central$beta, dat$beta, tol = 1e-3)
    )
    cr <- cr + 1

    for (method in methods) {
      message(sprintf("seed=%d lambda=%.3f method=%s", seed, lambda, method))
      all_runs[[k]] <- run_method(method, dat, lambda, target_obj, seed, rounds = 800)
      k <- k + 1
    }
  }
}

summary_tbl <- do.call(rbind, lapply(all_runs, `[[`, "summary"))
trace_tbl <- do.call(rbind, lapply(lapply(all_runs, `[[`, "trace"), function(d) {
  keep <- intersect(c("round", "objective", "method", "lambda", "seed", "gap"), names(d))
  d[, keep]
}))
central_tbl <- do.call(rbind, central_rows)
summary_all <- rbind(summary_tbl, central_tbl)

agg_tbl <- aggregate(
  cbind(final_gap, selected_size, tp, fp, fn, tpr, fdr, l2_error) ~ method + lambda,
  data = summary_all,
  FUN = mean
)

trace_agg <- aggregate(
  cbind(objective, gap) ~ method + lambda + round,
  data = trace_tbl,
  FUN = mean
)

write.csv(summary_all,
          file.path(root, "results", "l1_penalized_summary.csv"),
          row.names = FALSE)
write.csv(agg_tbl,
          file.path(root, "results", "l1_penalized_aggregate.csv"),
          row.names = FALSE)
write.csv(trace_tbl,
          file.path(root, "results", "l1_penalized_trace.csv"),
          row.names = FALSE)
write.csv(trace_agg,
          file.path(root, "results", "l1_penalized_trace_aggregate.csv"),
          row.names = FALSE)

method_cols <- c(
  "QR box-dual" = "#2A6FBB",
  "FSPG-smooth" = "#984EA3",
  "FedSubGrad" = "#4DAF4A",
  "FedQR-ADMM" = "#E6AB02",
  "Central L1 QR box-dual" = "#555555"
)

plot_l1_gap <- function() {
  png(file.path(root, "figures", "l1_penalized_final_gap.png"),
      width = 1600, height = 1100, res = 180)
  op <- par(mar = c(5, 5, 3, 1))
  on.exit(par(op), add = TRUE)
  plot(
    NA,
    xlim = range(lambdas),
    ylim = c(0, max(subset(agg_tbl, method != "Central L1 QR box-dual")$final_gap) * 1.05),
    xlab = "lambda",
    ylab = "Final objective gap",
    main = "L1 penalized hard non-IID R1+R2"
  )
  for (m in setdiff(names(method_cols), "Central L1 QR box-dual")) {
    d <- subset(agg_tbl, method == m)
    lines(d$lambda, d$final_gap, type = "b", pch = 19, col = method_cols[[m]], lwd = 2)
  }
  legend("topleft", legend = setdiff(names(method_cols), "Central L1 QR box-dual"),
         col = method_cols[setdiff(names(method_cols), "Central L1 QR box-dual")],
         lwd = 2, pch = 19, bty = "n", cex = 0.8)
  dev.off()
}

plot_l1_selection <- function() {
  png(file.path(root, "figures", "l1_penalized_selection.png"),
      width = 1800, height = 900, res = 180)
  op <- par(mar = c(5, 5, 3, 1), mfrow = c(1, 2))
  on.exit(par(op), add = TRUE)
  for (metric in c("tpr", "fdr")) {
    plot(
      NA,
      xlim = range(lambdas),
      ylim = c(0, 1),
      xlab = "lambda",
      ylab = toupper(metric),
      main = paste("Support", toupper(metric))
    )
    for (m in names(method_cols)) {
      d <- subset(agg_tbl, method == m)
      lines(d$lambda, d[[metric]], type = "b", pch = 19, col = method_cols[[m]], lwd = 2)
    }
    legend("right", legend = names(method_cols), col = method_cols, lwd = 2, pch = 19,
           bty = "n", cex = 0.72)
  }
  dev.off()
}

plot_l1_support_size <- function() {
  png(file.path(root, "figures", "l1_penalized_support_size.png"),
      width = 1600, height = 1100, res = 180)
  op <- par(mar = c(5, 5, 3, 1))
  on.exit(par(op), add = TRUE)
  plot(
    NA,
    xlim = range(lambdas),
    ylim = c(0, max(agg_tbl$selected_size) * 1.1),
    xlab = "lambda",
    ylab = "Selected support size",
    main = "L1 selected variables"
  )
  for (m in names(method_cols)) {
    d <- subset(agg_tbl, method == m)
    lines(d$lambda, d$selected_size, type = "b", pch = 19, col = method_cols[[m]], lwd = 2)
  }
  abline(h = 8, lty = 2, col = "#333333")
  legend("topright", legend = names(method_cols), col = method_cols, lwd = 2, pch = 19,
         bty = "n", cex = 0.72)
  dev.off()
}

plot_l1_gap()
plot_l1_selection()
plot_l1_support_size()

cat("\nL1 penalized aggregate:\n")
print(agg_tbl)

cat("\nFigures written to:\n")
cat(file.path(root, "figures", "l1_penalized_final_gap.png"), "\n")
cat(file.path(root, "figures", "l1_penalized_selection.png"), "\n")
cat(file.path(root, "figures", "l1_penalized_support_size.png"), "\n")

