args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
} else {
  file.path(getwd(), "scripts", "run_nonconvex_penalty_experiment_with_plots.R")
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

source(file.path(root, "R", "prox.R"))
source(file.path(root, "R", "objective.R"))
source(file.path(root, "R", "simulate.R"))
source(file.path(root, "R", "qr_box_fed_pdhg.R"))

dir.create(file.path(root, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "figures"), recursive = TRUE, showWarnings = FALSE)

support_metrics <- function(beta_hat, beta_true, tol = 1e-3) {
  keep <- seq_along(beta_true)[-1]
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
    tpr = tp / sum(true_support),
    fdr = if (sum(selected) == 0) 0 else fp / sum(selected),
    l2_error = sqrt(sum((beta_hat - beta_true)^2))
  )
}

run_box_penalty <- function(dat, penalty, lambda, seed, target_obj, rounds = 1000) {
  fit <- qr_box_fed_pdhg(
    dat$X, dat$y,
    client_indices = dat$client_indices,
    tau = dat$tau,
    lambda = lambda,
    penalty = penalty,
    rounds = rounds,
    clients_per_round = 4,
    batch_size = 10,
    step_rule = "operator",
    aggregation = "cached",
    beta_ref = dat$beta,
    trace_every = 5,
    seed = seed
  )

  trace <- transform(
    fit$trace,
    penalty = penalty,
    lambda = lambda,
    seed = seed,
    gap = objective - target_obj
  )

  summary <- cbind(
    data.frame(
      penalty = penalty,
      lambda = lambda,
      seed = seed,
      final_objective = fit$objective,
      final_gap = fit$objective - target_obj
    ),
    support_metrics(fit$beta, dat$beta)
  )

  list(summary = summary, trace = trace, beta = fit$beta)
}

seeds <- c(20260526, 20260527, 20260528)
penalties <- c("l1", "mcp", "scad")
lambdas <- c(0.005, 0.01, 0.02, 0.04, 0.08)

runs <- list()
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

  for (penalty in penalties) {
    for (lambda in lambdas) {
      message(sprintf("seed=%d penalty=%s lambda=%.3f", seed, penalty, lambda))
      central <- qr_box_fed_pdhg(
        dat$X, dat$y,
        client_indices = list(seq_len(nrow(dat$X))),
        tau = dat$tau,
        lambda = lambda,
        penalty = penalty,
        rounds = 3500,
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
          penalty = paste0("central_", penalty),
          lambda = lambda,
          seed = seed,
          final_objective = central$objective,
          final_gap = 0
        ),
        support_metrics(central$beta, dat$beta)
      )
      cr <- cr + 1

      runs[[k]] <- run_box_penalty(dat, penalty, lambda, seed, target_obj)
      k <- k + 1
    }
  }
}

summary_tbl <- do.call(rbind, lapply(runs, `[[`, "summary"))
trace_tbl <- do.call(rbind, lapply(lapply(runs, `[[`, "trace"), function(d) {
  d[, c("round", "objective", "penalty", "lambda", "seed", "gap")]
}))
central_tbl <- do.call(rbind, central_rows)
summary_all <- rbind(summary_tbl, central_tbl)

agg_tbl <- aggregate(
  cbind(final_gap, selected_size, tp, fp, fn, tpr, fdr, l2_error) ~ penalty + lambda,
  data = summary_all,
  FUN = mean
)

trace_agg <- aggregate(
  cbind(objective, gap) ~ penalty + lambda + round,
  data = trace_tbl,
  FUN = mean
)

write.csv(summary_all,
          file.path(root, "results", "nonconvex_penalty_summary.csv"),
          row.names = FALSE)
write.csv(agg_tbl,
          file.path(root, "results", "nonconvex_penalty_aggregate.csv"),
          row.names = FALSE)
write.csv(trace_tbl,
          file.path(root, "results", "nonconvex_penalty_trace.csv"),
          row.names = FALSE)
write.csv(trace_agg,
          file.path(root, "results", "nonconvex_penalty_trace_aggregate.csv"),
          row.names = FALSE)

cols <- c("l1" = "#2A6FBB", "mcp" = "#D95F02", "scad" = "#984EA3")

plot_penalty_gap <- function() {
  png(file.path(root, "figures", "nonconvex_penalty_final_gap.png"),
      width = 1600, height = 1100, res = 180)
  op <- par(mar = c(5, 5, 3, 1))
  on.exit(par(op), add = TRUE)
  d0 <- subset(agg_tbl, penalty %in% penalties)
  plot(
    NA,
    xlim = range(lambdas),
    ylim = c(0, max(d0$final_gap) * 1.05),
    xlab = "lambda",
    ylab = "Final objective gap vs central same penalty",
    main = "QR box-dual with L1/MCP/SCAD"
  )
  for (pen in penalties) {
    d <- subset(agg_tbl, penalty == pen)
    lines(d$lambda, d$final_gap, type = "b", pch = 19, col = cols[[pen]], lwd = 2)
  }
  legend("topright", legend = penalties, col = cols, lwd = 2, pch = 19, bty = "n")
  dev.off()
}

plot_penalty_selection <- function() {
  png(file.path(root, "figures", "nonconvex_penalty_selection.png"),
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
    for (pen in penalties) {
      d <- subset(agg_tbl, penalty == pen)
      lines(d$lambda, d[[metric]], type = "b", pch = 19, col = cols[[pen]], lwd = 2)
    }
    legend("right", legend = penalties, col = cols, lwd = 2, pch = 19, bty = "n")
  }
  dev.off()
}

plot_penalty_support <- function() {
  png(file.path(root, "figures", "nonconvex_penalty_support_size.png"),
      width = 1600, height = 1100, res = 180)
  op <- par(mar = c(5, 5, 3, 1))
  on.exit(par(op), add = TRUE)
  d0 <- subset(agg_tbl, penalty %in% penalties)
  plot(
    NA,
    xlim = range(lambdas),
    ylim = c(0, max(d0$selected_size) * 1.1),
    xlab = "lambda",
    ylab = "Selected support size",
    main = "Selected variables by penalty"
  )
  for (pen in penalties) {
    d <- subset(agg_tbl, penalty == pen)
    lines(d$lambda, d$selected_size, type = "b", pch = 19, col = cols[[pen]], lwd = 2)
  }
  abline(h = 8, lty = 2, col = "#333333")
  legend("topright", legend = penalties, col = cols, lwd = 2, pch = 19, bty = "n")
  dev.off()
}

plot_penalty_gap()
plot_penalty_selection()
plot_penalty_support()

cat("\nNonconvex penalty aggregate:\n")
print(agg_tbl)

cat("\nFigures written to:\n")
cat(file.path(root, "figures", "nonconvex_penalty_final_gap.png"), "\n")
cat(file.path(root, "figures", "nonconvex_penalty_selection.png"), "\n")
cat(file.path(root, "figures", "nonconvex_penalty_support_size.png"), "\n")
