args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
} else {
  file.path(getwd(), "scripts", "run_large_client_scale_stress_with_plots.R")
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

source(file.path(root, "scripts", "load_rfedqr.R"))
make_experiment_dirs(root)

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) as.integer(value) else default
}

env_num_vec <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }
  as.numeric(strsplit(value, ",", fixed = TRUE)[[1]])
}

rounds <- env_int("SCALE_STRESS_ROUNDS", 250)
trace_every <- env_int("SCALE_STRESS_TRACE_EVERY", 25)
target_iter <- env_int("SCALE_STRESS_TARGET_ITER", 5000)
seed_count <- env_int("SCALE_STRESS_SEEDS", 2)
client_counts <- env_num_vec("SCALE_STRESS_CLIENTS", c(50, 100))
taus <- env_num_vec("SCALE_STRESS_TAUS", c(0.9, 0.95))
heterogeneities <- strsplit(Sys.getenv("SCALE_STRESS_HETEROGENEITY", unset = "hard,extreme"),
                            ",", fixed = TRUE)[[1]]
sim_seeds <- seq(20260630, length.out = seed_count)

methods <- c(
  "QR box-dual",
  "QR box-dual stale",
  "QR box-dual robust",
  "QR box-dual stale+robust",
  "QR box-dual adaptive",
  "FSPG-smooth",
  "FedSPD-check"
)

method_controls <- list(
  "QR box-dual" = list(step_rule = "box"),
  "QR box-dual stale" = list(
    step_rule = "box",
    staleness_rate = 0.05,
    staleness_floor = 0.35
  ),
  "QR box-dual robust" = list(
    step_rule = "box",
    client_weighting = "uniform"
  ),
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
  "FedSPD-check" = list(Q = 3, gamma0 = 4)
)

target_for <- function(X, y, tau, seed) {
  fit <- qr_pdhg(
    X, y,
    tau = tau,
    max_iter = target_iter,
    step_rule = "box",
    trace_every = 100,
    seed = seed
  )
  list(beta = fit$beta, objective = fit$objective)
}

make_unbalanced_n <- function(n_clients, seed) {
  set.seed(seed)
  base <- round(exp(stats::rnorm(n_clients, log(70), 0.55)))
  pmin(220, pmax(25, base))
}

collect_run <- function(dat, target, client_count, heterogeneity, tau, seed,
                        clients_per_round, batch_size) {
  result <- run_fedqr_methods(
    methods,
    dat$X,
    dat$y,
    client_indices = dat$client_indices,
    tau = tau,
    rounds = rounds,
    clients_per_round = clients_per_round,
    batch_size = batch_size,
    seed = seed,
    trace_every = trace_every,
    method_controls = method_controls,
    term_names = colnames(dat$X),
    verbose = TRUE
  )

  summary <- result$summary
  summary$client_count <- client_count
  summary$clients_per_round <- clients_per_round
  summary$batch_size <- batch_size
  summary$heterogeneity <- heterogeneity
  summary$tau <- tau
  summary$seed <- seed
  summary$target_objective <- target$objective
  summary$target_gap <- summary$objective - target$objective
  summary <- summary[, c(
    "client_count", "clients_per_round", "batch_size", "heterogeneity",
    "tau", "seed", "method", "objective", "target_objective",
    "target_gap", "gap_to_best_observed", "beta_norm"
  )]

  trace <- result$trace
  trace$client_count <- client_count
  trace$clients_per_round <- clients_per_round
  trace$batch_size <- batch_size
  trace$heterogeneity <- heterogeneity
  trace$tau <- tau
  trace$seed <- seed
  trace$target_gap <- trace$objective - target$objective

  client_loss <- do.call(rbind, lapply(result$fits, function(fit) {
    out <- client_loss_summary(dat$X, dat$y, fit$beta, dat$client_indices, tau = tau)
    out$client_count <- client_count
    out$clients_per_round <- clients_per_round
    out$batch_size <- batch_size
    out$heterogeneity <- heterogeneity
    out$tau <- tau
    out$seed <- seed
    out$method <- fit$method
    out
  }))
  client_loss <- client_loss[, c(
    "client_count", "clients_per_round", "batch_size", "heterogeneity",
    "tau", "seed", "method", "global_mean_loss", "client_mean_loss",
    "worst_client_loss", "client_loss_sd", "client_q90_loss", "min_client_loss"
  )]

  design <- data.frame(
    client_count = client_count,
    clients_per_round = clients_per_round,
    participation_rate = clients_per_round / client_count,
    batch_size = batch_size,
    heterogeneity = heterogeneity,
    tau = tau,
    seed = seed,
    n = nrow(dat$X),
    p = ncol(dat$X),
    min_client_n = min(lengths(dat$client_indices)),
    median_client_n = stats::median(lengths(dat$client_indices)),
    max_client_n = max(lengths(dat$client_indices)),
    stringsAsFactors = FALSE
  )

  list(summary = summary, trace = trace, client_loss = client_loss, design = design)
}

runs <- list()
k <- 1
for (client_count in client_counts) {
  for (heterogeneity in heterogeneities) {
    for (tau in taus) {
      for (seed in sim_seeds) {
        n_per_client <- make_unbalanced_n(client_count, seed + client_count)
        dat <- make_hard_federated_qr_sim(
          n_clients = client_count,
          n_per_client = n_per_client,
          p = 12,
          tau = tau,
          heterogeneity = heterogeneity,
          seed = seed
        )
        target <- target_for(dat$X, dat$y, tau, seed)
        clients_per_round <- max(5, ceiling(0.10 * client_count))
        batch_size <- 25
        message(sprintf(
          "Scale stress: clients=%d, K=%d, tau=%.2f, heterogeneity=%s, seed=%d, n=%d",
          client_count, clients_per_round, tau, heterogeneity, seed, nrow(dat$X)
        ))
        runs[[k]] <- collect_run(
          dat = dat,
          target = target,
          client_count = client_count,
          heterogeneity = heterogeneity,
          tau = tau,
          seed = seed,
          clients_per_round = clients_per_round,
          batch_size = batch_size
        )
        k <- k + 1
      }
    }
  }
}

summary_tbl <- do.call(rbind, lapply(runs, `[[`, "summary"))
trace_tbl <- do.call(rbind, lapply(runs, `[[`, "trace"))
client_loss_tbl <- do.call(rbind, lapply(runs, `[[`, "client_loss"))
design_tbl <- do.call(rbind, lapply(runs, `[[`, "design"))

write.csv(summary_tbl, file.path(root, "results", "scale_stress_summary.csv"),
          row.names = FALSE)
write.csv(trace_tbl, file.path(root, "results", "scale_stress_trace.csv"),
          row.names = FALSE)
write.csv(client_loss_tbl, file.path(root, "results", "scale_stress_client_loss.csv"),
          row.names = FALSE)
write.csv(design_tbl, file.path(root, "results", "scale_stress_design.csv"),
          row.names = FALSE)

cols <- fedqr_default_method_cols(methods)

agg_gap <- aggregate(target_gap ~ client_count + heterogeneity + tau + method,
                     summary_tbl, mean)
agg_gap <- agg_gap[order(agg_gap$client_count, agg_gap$heterogeneity,
                         agg_gap$tau, agg_gap$target_gap), ]

agg_fair <- aggregate(worst_client_loss ~ client_count + heterogeneity + tau + method,
                      client_loss_tbl, mean)
agg_fair <- agg_fair[order(agg_fair$client_count, agg_fair$heterogeneity,
                           agg_fair$tau, agg_fair$worst_client_loss), ]

plot_panel_bars <- function(tbl, value, file, main_prefix, ylab) {
  panels <- unique(tbl[, c("client_count", "heterogeneity", "tau")])
  panels <- panels[order(panels$client_count, panels$heterogeneity, panels$tau), ]
  grDevices::png(file, width = 1800, height = 1200, res = 180)
  old <- par(no.readonly = TRUE)
  on.exit({
    par(old)
    grDevices::dev.off()
  })
  par(mfrow = c(2, ceiling(nrow(panels) / 2)), mar = c(8, 4, 3, 1))
  for (i in seq_len(nrow(panels))) {
    panel <- panels[i, ]
    d <- tbl[tbl$client_count == panel$client_count &
      tbl$heterogeneity == panel$heterogeneity &
      tbl$tau == panel$tau, ]
    d <- d[order(d[[value]]), ]
    barplot(
      height = d[[value]],
      names.arg = d$method,
      las = 2,
      col = cols[d$method],
      main = sprintf("%s: m=%d, %s, tau=%.2f",
                     main_prefix, panel$client_count,
                     panel$heterogeneity, panel$tau),
      ylab = ylab,
      cex.names = 0.55
    )
  }
}

plot_panel_bars(
  agg_gap,
  "target_gap",
  file.path(root, "figures", "scale_stress_final_gap.png"),
  "Final gap",
  "Mean target gap"
)
plot_panel_bars(
  agg_gap[agg_gap$method != "FedSPD-check", ],
  "target_gap",
  file.path(root, "figures", "scale_stress_final_gap_zoom.png"),
  "Final gap zoom",
  "Mean target gap"
)
plot_panel_bars(
  agg_fair,
  "worst_client_loss",
  file.path(root, "figures", "scale_stress_client_fairness.png"),
  "Fairness",
  "Mean worst-client loss"
)
plot_panel_bars(
  agg_fair[agg_fair$method != "FedSPD-check", ],
  "worst_client_loss",
  file.path(root, "figures", "scale_stress_client_fairness_zoom.png"),
  "Fairness zoom",
  "Mean worst-client loss"
)

adaptive_trace <- trace_tbl[trace_tbl$method == "QR box-dual adaptive" &
  !is.na(trace_tbl$adaptive_staleness_rate), ]
if (nrow(adaptive_trace) > 0) {
  adaptive_agg <- aggregate(
    cbind(adaptive_staleness_rate, adaptive_client_weight_sd) ~
      client_count + heterogeneity + tau + round,
    adaptive_trace,
    mean
  )
  grDevices::png(file.path(root, "figures", "scale_stress_adaptive_trace.png"),
                 width = 1600, height = 900, res = 180)
  old <- par(no.readonly = TRUE)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  groups <- interaction(adaptive_agg$client_count, adaptive_agg$heterogeneity,
                        adaptive_agg$tau, drop = TRUE)
  plot(NULL, xlim = range(adaptive_agg$round),
       ylim = range(adaptive_agg$adaptive_staleness_rate),
       xlab = "Round", ylab = "Adaptive staleness rate",
       main = "Staleness adaptation")
  palette_trace <- grDevices::rainbow(length(levels(groups)))
  for (g in seq_along(levels(groups))) {
    d <- adaptive_agg[groups == levels(groups)[g], ]
    lines(d$round, d$adaptive_staleness_rate, col = palette_trace[g], lwd = 1.8)
  }
  plot(NULL, xlim = range(adaptive_agg$round),
       ylim = range(adaptive_agg$adaptive_client_weight_sd),
       xlab = "Round", ylab = "Client weight sd",
       main = "Client-weight adaptation")
  for (g in seq_along(levels(groups))) {
    d <- adaptive_agg[groups == levels(groups)[g], ]
    lines(d$round, d$adaptive_client_weight_sd, col = palette_trace[g], lwd = 1.8)
  }
  par(old)
  grDevices::dev.off()
}

cat("\nScale stress target-gap aggregate:\n")
print(agg_gap, row.names = FALSE)

cat("\nScale stress fairness aggregate:\n")
print(agg_fair, row.names = FALSE)

cat("\nScale stress design summary:\n")
print(design_tbl, row.names = FALSE)
