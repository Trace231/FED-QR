args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
} else {
  file.path(getwd(), "scripts", "run_expanded_scenario_suite_with_plots.R")
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

make_unbalanced_sizes <- function() {
  c(30, 35, 40, 45, 50, 60, 70, 80, 90, 100,
    110, 120, 130, 140, 150, 160, 170, 180, 190, 200)
}

make_scenario_data <- function(scenario, seed) {
  if (scenario == "low_participation") {
    dat <- make_hard_federated_qr_sim(
      n_clients = 20, n_per_client = 80, p = 15,
      tau = 0.9, heterogeneity = "hard", seed = seed
    )
    return(list(dat = dat, K = 2, batch_size = 10, label = "K=2, batch=10"))
  }

  if (scenario == "tiny_user_batch") {
    dat <- make_hard_federated_qr_sim(
      n_clients = 20, n_per_client = 80, p = 15,
      tau = 0.9, heterogeneity = "hard", seed = seed
    )
    return(list(dat = dat, K = 4, batch_size = 5, label = "K=4, batch=5"))
  }

  if (scenario == "extreme_tau") {
    dat <- make_hard_federated_qr_sim(
      n_clients = 20, n_per_client = 80, p = 15,
      tau = 0.95, heterogeneity = "hard", seed = seed
    )
    return(list(dat = dat, K = 4, batch_size = 10, label = "tau=0.95"))
  }

  if (scenario == "extreme_heterogeneity") {
    dat <- make_hard_federated_qr_sim(
      n_clients = 20, n_per_client = 80, p = 15,
      tau = 0.9, heterogeneity = "extreme", seed = seed
    )
    return(list(dat = dat, K = 4, batch_size = 10, label = "extreme heterogeneity"))
  }

  if (scenario == "unbalanced_clients") {
    dat <- make_hard_federated_qr_sim(
      n_clients = 20, n_per_client = make_unbalanced_sizes(), p = 15,
      tau = 0.9, heterogeneity = "hard", seed = seed
    )
    return(list(dat = dat, K = 4, batch_size = 10, label = "unbalanced client sizes"))
  }

  if (scenario == "high_dim_sparse") {
    dat <- make_hard_federated_qr_sim(
      n_clients = 20, n_per_client = 80, p = 60, sparsity = 8,
      tau = 0.9, heterogeneity = "hard", seed = seed
    )
    return(list(dat = dat, K = 4, batch_size = 10, label = "p=60 sparse"))
  }

  stop("Unknown scenario: ", scenario)
}

standard_trace <- function(trace, method, scenario, seed, target_obj) {
  data.frame(
    round = trace$round,
    objective = trace$objective,
    method = method,
    scenario = scenario,
    seed = seed,
    gap = trace$objective - target_obj
  )
}

run_method <- function(method, dat, scenario, K, batch_size, target, seed, rounds = 400) {
  tau <- dat$tau

  if (method == "QR box-dual") {
    fit <- qr_box_fed_pdhg(
      dat$X, dat$y,
      client_indices = dat$client_indices,
      tau = tau,
      rounds = rounds,
      clients_per_round = K,
      batch_size = batch_size,
      step_rule = "operator",
      aggregation = "cached",
      beta_ref = target$beta,
      trace_every = 5,
      seed = seed
    )
    beta <- fit$beta
    obj <- fit$objective
    trace <- standard_trace(fit$trace, method, scenario, seed, target$objective)
  } else if (method == "FSPG-smooth") {
    fit <- fed_smooth_qr(
      dat$X, dat$y,
      client_indices = dat$client_indices,
      tau = tau,
      rounds = rounds,
      clients_per_round = K,
      batch_size = batch_size,
      smooth_mu = 0.05,
      beta_true = target$beta,
      trace_every = 5,
      seed = seed
    )
    beta <- fit$beta
    obj <- fit$objective
    trace <- standard_trace(fit$trace, method, scenario, seed, target$objective)
  } else if (method == "FedSubGrad") {
    fit <- fed_subgrad_qr(
      dat$X, dat$y,
      client_indices = dat$client_indices,
      tau = tau,
      rounds = rounds,
      clients_per_round = K,
      batch_size = batch_size,
      beta_true = target$beta,
      trace_every = 5,
      seed = seed
    )
    beta <- fit$beta
    obj <- fit$objective
    trace <- standard_trace(fit$trace, method, scenario, seed, target$objective)
  } else if (method == "FedSPD-check") {
    fit <- fedspd_dp_qr(
      dat$X, dat$y,
      client_indices = dat$client_indices,
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
      beta_ref = target$beta,
      trace_every = 5,
      seed = seed
    )
    beta <- fit$beta
    obj <- fit$qr_objective
    trace <- data.frame(
      round = fit$trace$round,
      objective = fit$trace$qr_objective,
      method = method,
      scenario = scenario,
      seed = seed,
      gap = fit$trace$qr_objective - target$objective
    )
  } else if (method == "FedQR-ADMM") {
    fit <- fed_qr_admm(
      dat$X, dat$y,
      client_indices = dat$client_indices,
      tau = tau,
      rounds = min(rounds, 300),
      clients_per_round = K,
      batch_size = batch_size,
      rho_consensus = 1,
      rho_residual = 0.3,
      inner_iter = 20,
      beta_ref = target$beta,
      trace_every = 5,
      seed = seed
    )
    beta <- fit$beta
    obj <- fit$objective
    trace <- standard_trace(fit$trace, method, scenario, seed, target$objective)
  } else {
    stop("Unknown method: ", method)
  }

  summary <- data.frame(
    scenario = scenario,
    seed = seed,
    method = method,
    tau = tau,
    heterogeneity = dat$heterogeneity,
    n = nrow(dat$X),
    p = ncol(dat$X),
    n_clients = length(dat$client_indices),
    K = K,
    batch_size = batch_size,
    final_objective = obj,
    final_gap = obj - target$objective,
    beta_l2_to_target = sqrt(sum((beta - target$beta)^2))
  )

  list(summary = summary, trace = trace)
}

scenario_ids <- c(
  "low_participation",
  "tiny_user_batch",
  "extreme_tau",
  "extreme_heterogeneity",
  "unbalanced_clients",
  "high_dim_sparse"
)
seeds <- c(20260526, 20260527)
methods <- c("QR box-dual", "FSPG-smooth", "FedSubGrad", "FedSPD-check", "FedQR-ADMM")
rounds <- 400

all_runs <- list()
target_rows <- list()
scenario_rows <- list()
k <- 1
tr <- 1
sr <- 1

for (seed in seeds) {
  for (scenario in scenario_ids) {
    setup <- make_scenario_data(scenario, seed)
    dat <- setup$dat

    central <- qr_pdhg(
      dat$X, dat$y,
      tau = dat$tau,
      max_iter = 5000,
      step_rule = "box",
      trace_every = 100,
      seed = seed
    )
    target <- list(beta = central$beta, objective = central$objective)
    target_rows[[tr]] <- data.frame(
      scenario = scenario,
      seed = seed,
      target_objective = central$objective,
      target_beta_l2_to_truth = sqrt(sum((central$beta - dat$beta)^2))
    )
    tr <- tr + 1

    scenario_rows[[sr]] <- data.frame(
      scenario = scenario,
      label = setup$label,
      seed = seed,
      tau = dat$tau,
      heterogeneity = dat$heterogeneity,
      n = nrow(dat$X),
      p = ncol(dat$X),
      n_clients = length(dat$client_indices),
      min_client_n = min(lengths(dat$client_indices)),
      median_client_n = median(lengths(dat$client_indices)),
      max_client_n = max(lengths(dat$client_indices)),
      K = setup$K,
      batch_size = setup$batch_size
    )
    sr <- sr + 1

    for (method in methods) {
      message(sprintf(
        "scenario=%s seed=%d method=%s",
        scenario, seed, method
      ))
      all_runs[[k]] <- run_method(
        method = method,
        dat = dat,
        scenario = scenario,
        K = setup$K,
        batch_size = setup$batch_size,
        target = target,
        seed = seed,
        rounds = rounds
      )
      k <- k + 1
    }
  }
}

summary_tbl <- do.call(rbind, lapply(all_runs, `[[`, "summary"))
trace_tbl <- do.call(rbind, lapply(all_runs, `[[`, "trace"))
target_tbl <- do.call(rbind, target_rows)
scenario_tbl <- do.call(rbind, scenario_rows)

agg_tbl <- aggregate(
  cbind(final_gap, beta_l2_to_target) ~ scenario + method,
  data = summary_tbl,
  FUN = mean
)
agg_sd <- aggregate(
  final_gap ~ scenario + method,
  data = summary_tbl,
  FUN = sd
)
names(agg_sd)[names(agg_sd) == "final_gap"] <- "final_gap_sd"
agg_tbl <- merge(agg_tbl, agg_sd, by = c("scenario", "method"), all.x = TRUE)

rank_tbl <- do.call(rbind, lapply(split(agg_tbl, agg_tbl$scenario), function(d) {
  d$rank <- rank(d$final_gap, ties.method = "min")
  d[order(d$rank), ]
}))

trace_agg <- aggregate(
  cbind(objective, gap) ~ scenario + method + round,
  data = trace_tbl,
  FUN = mean
)

write.csv(scenario_tbl,
          file.path(root, "results", "expanded_scenario_design.csv"),
          row.names = FALSE)
write.csv(target_tbl,
          file.path(root, "results", "expanded_scenario_targets.csv"),
          row.names = FALSE)
write.csv(summary_tbl,
          file.path(root, "results", "expanded_scenario_summary.csv"),
          row.names = FALSE)
write.csv(agg_tbl,
          file.path(root, "results", "expanded_scenario_aggregate.csv"),
          row.names = FALSE)
write.csv(rank_tbl,
          file.path(root, "results", "expanded_scenario_ranks.csv"),
          row.names = FALSE)
write.csv(trace_tbl,
          file.path(root, "results", "expanded_scenario_trace.csv"),
          row.names = FALSE)
write.csv(trace_agg,
          file.path(root, "results", "expanded_scenario_trace_aggregate.csv"),
          row.names = FALSE)

method_cols <- c(
  "QR box-dual" = "#2A6FBB",
  "FSPG-smooth" = "#984EA3",
  "FedSubGrad" = "#4DAF4A",
  "FedSPD-check" = "#D95F02",
  "FedQR-ADMM" = "#E6AB02"
)

scenario_labels <- aggregate(
  label ~ scenario,
  data = scenario_tbl,
  FUN = function(x) x[1]
)

plot_final_gap <- function() {
  png(file.path(root, "figures", "expanded_scenario_final_gap.png"),
      width = 2200, height = 1500, res = 180)
  op <- par(mar = c(7, 5, 4, 1), mfrow = c(2, 3))
  on.exit(par(op), add = TRUE)
  for (scenario in scenario_ids) {
    d <- agg_tbl[agg_tbl$scenario == scenario, ]
    d <- d[match(names(method_cols), d$method), ]
    vals <- pmax(d$final_gap, 1e-8)
    names(vals) <- d$method
    label <- scenario_labels$label[scenario_labels$scenario == scenario]
    barplot(
      log10(vals),
      col = method_cols[d$method],
      las = 2,
      ylab = "log10(final gap)",
      main = label
    )
  }
  dev.off()
}

plot_rank_heatmap <- function() {
  png(file.path(root, "figures", "expanded_scenario_rank_heatmap.png"),
      width = 1600, height = 1000, res = 180)
  rank_mat <- matrix(NA_real_, nrow = length(method_cols), ncol = length(scenario_ids),
                     dimnames = list(names(method_cols), scenario_ids))
  for (scenario in scenario_ids) {
    d <- rank_tbl[rank_tbl$scenario == scenario, ]
    rank_mat[d$method, scenario] <- d$rank
  }
  op <- par(mar = c(8, 8, 4, 2))
  on.exit(par(op), add = TRUE)
  image(
    x = seq_len(ncol(rank_mat)),
    y = seq_len(nrow(rank_mat)),
    z = t(rank_mat[nrow(rank_mat):1, ]),
    col = hcl.colors(5, "YlOrRd", rev = TRUE),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = "Method rank by expanded scenario (1 is best)"
  )
  axis(1, at = seq_len(ncol(rank_mat)), labels = scenario_labels$label[match(scenario_ids, scenario_labels$scenario)],
       las = 2, cex.axis = 0.75)
  axis(2, at = seq_len(nrow(rank_mat)), labels = rev(rownames(rank_mat)), las = 2)
  for (i in seq_len(ncol(rank_mat))) {
    for (j in seq_len(nrow(rank_mat))) {
      text(i, nrow(rank_mat) - j + 1, labels = rank_mat[j, i], cex = 0.9)
    }
  }
  dev.off()
}

plot_convergence <- function() {
  png(file.path(root, "figures", "expanded_scenario_convergence.png"),
      width = 2200, height = 1500, res = 180)
  op <- par(mar = c(5, 5, 4, 1), mfrow = c(2, 3))
  on.exit(par(op), add = TRUE)
  for (scenario in scenario_ids) {
    d <- trace_agg[trace_agg$scenario == scenario, ]
    y_lim <- quantile(pmax(d$gap, 0), 0.98, na.rm = TRUE)
    y_lim <- max(y_lim, 1e-4)
    label <- scenario_labels$label[scenario_labels$scenario == scenario]
    plot(
      NA,
      xlim = range(d$round),
      ylim = c(0, y_lim),
      xlab = "Communication round",
      ylab = "Objective gap",
      main = label
    )
    for (method in names(method_cols)) {
      dm <- d[d$method == method, ]
      lines(dm$round, pmax(dm$gap, 0), col = method_cols[[method]], lwd = 2)
    }
    legend("topright", legend = names(method_cols), col = method_cols, lwd = 2, bty = "n", cex = 0.7)
  }
  dev.off()
}

plot_final_gap()
plot_rank_heatmap()
plot_convergence()

cat("\nExpanded scenario ranks:\n")
print(rank_tbl[order(rank_tbl$scenario, rank_tbl$rank), ])

cat("\nFigures written to:\n")
cat(file.path(root, "figures", "expanded_scenario_final_gap.png"), "\n")
cat(file.path(root, "figures", "expanded_scenario_rank_heatmap.png"), "\n")
cat(file.path(root, "figures", "expanded_scenario_convergence.png"), "\n")
