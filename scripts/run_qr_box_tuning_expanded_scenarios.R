args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
} else {
  file.path(getwd(), "scripts", "run_qr_box_tuning_expanded_scenarios.R")
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

source(file.path(root, "R", "prox.R"))
source(file.path(root, "R", "objective.R"))
source(file.path(root, "R", "simulate.R"))
source(file.path(root, "R", "qr_pdhg.R"))
source(file.path(root, "R", "qr_box_fed_pdhg.R"))
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
    return(list(dat = dat, K = 4, batch_size = 10, label = "unbalanced clients"))
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

default_steps <- function(X, tau, step_rule = c("operator", "box")) {
  step_rule <- match.arg(step_rule)
  op_norm <- svd(X / sqrt(nrow(X)), nu = 0, nv = 0)$d[1]
  op_norm <- max(op_norm, .Machine$double.eps)
  base <- 0.95 / op_norm
  eta <- base
  sigma <- base
  if (step_rule == "box") {
    vmax <- max(tau, 1 - tau)
    eta <- eta / vmax
    sigma <- sigma * vmax
  }
  c(eta = eta, sigma = sigma)
}

make_variant <- function(name, method = "QR box-dual", rounds = 400,
                         step_rule = "operator", eta_scale = 1,
                         sigma_scale = 1, primal_clip = NULL) {
  list(
    name = name,
    method = method,
    rounds = rounds,
    step_rule = step_rule,
    eta_scale = eta_scale,
    sigma_scale = sigma_scale,
    primal_clip = primal_clip
  )
}

variants <- list(
  make_variant("QR default 400"),
  make_variant("QR long 1200", rounds = 1200),
  make_variant("QR damp both 800", rounds = 800, eta_scale = 0.5, sigma_scale = 0.5),
  make_variant("QR low sigma 800", rounds = 800, eta_scale = 1, sigma_scale = 0.5),
  make_variant("QR box-aware 800", rounds = 800, step_rule = "box"),
  make_variant("QR clip .25 800", rounds = 800, primal_clip = 0.25),
  make_variant("FSPG smooth 400", method = "FSPG-smooth", rounds = 400)
)

run_variant <- function(variant, dat, target, scenario, K, batch_size, seed) {
  if (variant$method == "FSPG-smooth") {
    fit <- fed_smooth_qr(
      dat$X, dat$y,
      client_indices = dat$client_indices,
      tau = dat$tau,
      rounds = variant$rounds,
      clients_per_round = K,
      batch_size = batch_size,
      smooth_mu = 0.05,
      beta_true = target$beta,
      trace_every = 5,
      seed = seed
    )
    trace <- data.frame(
      round = fit$trace$round,
      objective = fit$trace$objective,
      gap = fit$trace$objective - target$objective,
      variant = variant$name,
      scenario = scenario,
      seed = seed
    )
    return(list(
      beta = fit$beta,
      objective = fit$objective,
      trace = trace
    ))
  }

  base <- default_steps(dat$X, dat$tau, step_rule = variant$step_rule)
  fit <- qr_box_fed_pdhg(
    dat$X, dat$y,
    client_indices = dat$client_indices,
    tau = dat$tau,
    rounds = variant$rounds,
    clients_per_round = K,
    batch_size = batch_size,
    step_rule = variant$step_rule,
    eta = base[["eta"]] * variant$eta_scale,
    sigma = base[["sigma"]] * variant$sigma_scale,
    aggregation = "cached",
    primal_clip = variant$primal_clip,
    beta_ref = target$beta,
    trace_every = 5,
    seed = seed
  )
  trace <- data.frame(
    round = fit$trace$round,
    objective = fit$trace$objective,
    gap = fit$trace$objective - target$objective,
    variant = variant$name,
    scenario = scenario,
    seed = seed
  )
  list(beta = fit$beta, objective = fit$objective, trace = trace)
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

all_summary <- list()
all_trace <- list()
scenario_rows <- list()
k <- 1
tk <- 1
sk <- 1

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

    scenario_rows[[sk]] <- data.frame(
      scenario = scenario,
      label = setup$label,
      seed = seed,
      tau = dat$tau,
      heterogeneity = dat$heterogeneity,
      n = nrow(dat$X),
      p = ncol(dat$X),
      K = setup$K,
      batch_size = setup$batch_size
    )
    sk <- sk + 1

    for (variant in variants) {
      message(sprintf(
        "tune scenario=%s seed=%d variant=%s",
        scenario, seed, variant$name
      ))
      fit <- run_variant(
        variant = variant,
        dat = dat,
        target = target,
        scenario = scenario,
        K = setup$K,
        batch_size = setup$batch_size,
        seed = seed
      )

      all_summary[[k]] <- data.frame(
        scenario = scenario,
        seed = seed,
        variant = variant$name,
        method = variant$method,
        rounds = variant$rounds,
        final_objective = fit$objective,
        final_gap = fit$objective - target$objective,
        beta_l2_to_target = sqrt(sum((fit$beta - target$beta)^2))
      )
      k <- k + 1

      all_trace[[tk]] <- fit$trace
      tk <- tk + 1
    }
  }
}

summary_tbl <- do.call(rbind, all_summary)
trace_tbl <- do.call(rbind, all_trace)
scenario_tbl <- do.call(rbind, scenario_rows)

agg_tbl <- aggregate(
  cbind(final_gap, beta_l2_to_target) ~ scenario + variant + method + rounds,
  data = summary_tbl,
  FUN = mean
)
agg_sd <- aggregate(
  final_gap ~ scenario + variant + method + rounds,
  data = summary_tbl,
  FUN = sd
)
names(agg_sd)[names(agg_sd) == "final_gap"] <- "final_gap_sd"
agg_tbl <- merge(agg_tbl, agg_sd,
                 by = c("scenario", "variant", "method", "rounds"),
                 all.x = TRUE)

rank_tbl <- do.call(rbind, lapply(split(agg_tbl, agg_tbl$scenario), function(d) {
  d$rank <- rank(d$final_gap, ties.method = "min")
  d[order(d$rank), ]
}))

trace_agg <- aggregate(
  cbind(objective, gap) ~ scenario + variant + round,
  data = trace_tbl,
  FUN = mean
)

write.csv(scenario_tbl,
          file.path(root, "results", "qr_tuning_expanded_design.csv"),
          row.names = FALSE)
write.csv(summary_tbl,
          file.path(root, "results", "qr_tuning_expanded_summary.csv"),
          row.names = FALSE)
write.csv(agg_tbl,
          file.path(root, "results", "qr_tuning_expanded_aggregate.csv"),
          row.names = FALSE)
write.csv(rank_tbl,
          file.path(root, "results", "qr_tuning_expanded_ranks.csv"),
          row.names = FALSE)
write.csv(trace_tbl,
          file.path(root, "results", "qr_tuning_expanded_trace.csv"),
          row.names = FALSE)
write.csv(trace_agg,
          file.path(root, "results", "qr_tuning_expanded_trace_aggregate.csv"),
          row.names = FALSE)

variant_cols <- c(
  "QR default 400" = "#2A6FBB",
  "QR long 1200" = "#1B9E77",
  "QR damp both 800" = "#7570B3",
  "QR low sigma 800" = "#E7298A",
  "QR box-aware 800" = "#D95F02",
  "QR clip .25 800" = "#66A61E",
  "FSPG smooth 400" = "#984EA3"
)

scenario_labels <- aggregate(label ~ scenario, data = scenario_tbl, FUN = function(x) x[1])

plot_gap <- function() {
  png(file.path(root, "figures", "qr_tuning_expanded_gap.png"),
      width = 2200, height = 1500, res = 180)
  op <- par(mar = c(8, 5, 4, 1), mfrow = c(2, 3))
  on.exit(par(op), add = TRUE)
  for (scenario in scenario_ids) {
    d <- agg_tbl[agg_tbl$scenario == scenario, ]
    d <- d[match(names(variant_cols), d$variant), ]
    vals <- pmax(d$final_gap, 1e-8)
    names(vals) <- d$variant
    barplot(
      log10(vals),
      col = variant_cols[d$variant],
      las = 2,
      ylab = "log10(final gap)",
      main = scenario_labels$label[scenario_labels$scenario == scenario]
    )
  }
  dev.off()
}

plot_rank <- function() {
  png(file.path(root, "figures", "qr_tuning_expanded_rank_heatmap.png"),
      width = 1600, height = 1050, res = 180)
  rank_mat <- matrix(NA_real_, nrow = length(variant_cols), ncol = length(scenario_ids),
                     dimnames = list(names(variant_cols), scenario_ids))
  for (scenario in scenario_ids) {
    d <- rank_tbl[rank_tbl$scenario == scenario, ]
    rank_mat[d$variant, scenario] <- d$rank
  }
  op <- par(mar = c(8, 10, 4, 2))
  on.exit(par(op), add = TRUE)
  image(seq_len(ncol(rank_mat)), seq_len(nrow(rank_mat)),
        t(rank_mat[nrow(rank_mat):1, ]),
        col = hcl.colors(length(variant_cols), "YlOrRd", rev = TRUE),
        axes = FALSE, xlab = "", ylab = "",
        main = "QR tuning and FSPG reference ranks")
  axis(1, at = seq_len(ncol(rank_mat)),
       labels = scenario_labels$label[match(scenario_ids, scenario_labels$scenario)],
       las = 2, cex.axis = 0.75)
  axis(2, at = seq_len(nrow(rank_mat)), labels = rev(rownames(rank_mat)), las = 2, cex.axis = 0.75)
  for (i in seq_len(ncol(rank_mat))) {
    for (j in seq_len(nrow(rank_mat))) {
      text(i, nrow(rank_mat) - j + 1, labels = rank_mat[j, i], cex = 0.85)
    }
  }
  dev.off()
}

plot_gap()
plot_rank()

cat("\nQR tuning expanded ranks:\n")
print(rank_tbl[order(rank_tbl$scenario, rank_tbl$rank), ])

cat("\nFigures written to:\n")
cat(file.path(root, "figures", "qr_tuning_expanded_gap.png"), "\n")
cat(file.path(root, "figures", "qr_tuning_expanded_rank_heatmap.png"), "\n")
