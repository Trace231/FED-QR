args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
} else {
  file.path(getwd(), "scripts", "run_advanced_stabilization_experiment_with_plots.R")
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

source(file.path(root, "R", "prox.R"))
source(file.path(root, "R", "objective.R"))
source(file.path(root, "R", "simulate.R"))
source(file.path(root, "R", "qr_pdhg.R"))
source(file.path(root, "R", "qr_box_fed_pdhg.R"))

dir.create(file.path(root, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "figures"), recursive = TRUE, showWarnings = FALSE)

make_variant <- function(label, step_rule = "operator", aggregation = "cached",
                         dual_relaxation = 1, server_momentum = 0,
                         step_decay_power = 0, step_decay_offset = 100,
                         primal_clip = NULL) {
  list(
    label = label,
    step_rule = step_rule,
    aggregation = aggregation,
    dual_relaxation = dual_relaxation,
    server_momentum = server_momentum,
    step_decay_power = step_decay_power,
    step_decay_offset = step_decay_offset,
    primal_clip = primal_clip
  )
}

variants <- list(
  make_variant("A0 operator"),
  make_variant("M1 box-aware", step_rule = "box"),
  make_variant("M2 tau-adaptive", step_rule = "tau_adaptive"),
  make_variant("Dual relax .5", dual_relaxation = 0.5),
  make_variant("Server EMA .5", server_momentum = 0.5),
  make_variant("Step decay .25", step_decay_power = 0.25, step_decay_offset = 50),
  make_variant("Direction clip .25", primal_clip = 0.25),
  make_variant("Relax + EMA", dual_relaxation = 0.5, server_momentum = 0.5),
  make_variant("Relax + decay", dual_relaxation = 0.7, step_decay_power = 0.20, step_decay_offset = 50),
  make_variant("Selected reweighted", aggregation = "selected_reweighted")
)

variant_labels <- vapply(variants, `[[`, character(1), "label")

run_variant <- function(dat, target, tau, heterogeneity, scenario, K, batch_size,
                        variant, seed, rounds = 500) {
  fit <- qr_box_fed_pdhg(
    dat$X, dat$y,
    client_indices = dat$client_indices,
    tau = tau,
    rounds = rounds,
    clients_per_round = K,
    batch_size = batch_size,
    step_rule = variant$step_rule,
    aggregation = variant$aggregation,
    dual_relaxation = variant$dual_relaxation,
    server_momentum = variant$server_momentum,
    step_decay_power = variant$step_decay_power,
    step_decay_offset = variant$step_decay_offset,
    primal_clip = variant$primal_clip,
    beta_ref = target$beta,
    trace_every = 5,
    seed = seed
  )

  trace <- data.frame(
    heterogeneity = heterogeneity,
    tau = tau,
    scenario = scenario,
    seed = seed,
    variant = variant$label,
    round = fit$trace$round,
    objective = fit$trace$objective,
    gap = fit$trace$objective - target$objective,
    eta = fit$trace$eta,
    sigma = fit$trace$sigma,
    direction_norm = fit$trace$direction_norm
  )

  summary <- data.frame(
    heterogeneity = heterogeneity,
    tau = tau,
    scenario = scenario,
    seed = seed,
    variant = variant$label,
    step_rule = variant$step_rule,
    aggregation = variant$aggregation,
    dual_relaxation = variant$dual_relaxation,
    server_momentum = variant$server_momentum,
    step_decay_power = variant$step_decay_power,
    primal_clip = if (is.null(variant$primal_clip)) NA_real_ else variant$primal_clip,
    final_objective = fit$objective,
    final_gap = fit$objective - target$objective,
    beta_l2_to_target = sqrt(sum((fit$beta - target$beta)^2)),
    final_direction_norm = tail(fit$trace$direction_norm, 1)
  )

  list(summary = summary, trace = trace)
}

taus <- c(0.9, 0.95)
heterogeneity_levels <- c("hard", "extreme")
seeds <- c(20260526, 20260527, 20260528)
n_clients <- 20
n_per_client <- 80
p <- 15
rounds <- 500

scenarios <- data.frame(
  scenario = c("Full", "R1 + R2"),
  K = c(n_clients, 4),
  batch_size = c(100000, 10)
)

all_runs <- list()
target_rows <- list()
k <- 1
tr <- 1

for (seed in seeds) {
  for (heterogeneity in heterogeneity_levels) {
    for (tau in taus) {
      dat <- make_hard_federated_qr_sim(
        n_clients = n_clients,
        n_per_client = n_per_client,
        p = p,
        tau = tau,
        heterogeneity = heterogeneity,
        seed = seed
      )

      central <- qr_pdhg(
        dat$X, dat$y,
        tau = tau,
        max_iter = 5000,
        step_rule = "box",
        trace_every = 100,
        seed = seed
      )
      target <- list(beta = central$beta, objective = central$objective)
      target_rows[[tr]] <- data.frame(
        heterogeneity = heterogeneity,
        tau = tau,
        seed = seed,
        target_objective = central$objective
      )
      tr <- tr + 1

      full_batch <- max(lengths(dat$client_indices))
      for (s in seq_len(nrow(scenarios))) {
        scenario <- scenarios$scenario[s]
        K <- scenarios$K[s]
        batch_size <- if (scenarios$batch_size[s] >= 100000) {
          full_batch
        } else {
          scenarios$batch_size[s]
        }

        for (variant in variants) {
          message(sprintf(
            "seed=%d hetero=%s tau=%.2f scenario=%s variant=%s",
            seed, heterogeneity, tau, scenario, variant$label
          ))
          all_runs[[k]] <- run_variant(
            dat = dat,
            target = target,
            tau = tau,
            heterogeneity = heterogeneity,
            scenario = scenario,
            K = K,
            batch_size = batch_size,
            variant = variant,
            seed = seed,
            rounds = rounds
          )
          k <- k + 1
        }
      }
    }
  }
}

summary_tbl <- do.call(rbind, lapply(all_runs, `[[`, "summary"))
trace_tbl <- do.call(rbind, lapply(all_runs, `[[`, "trace"))
target_tbl <- do.call(rbind, target_rows)

agg_tbl <- aggregate(
  cbind(final_gap, beta_l2_to_target, final_direction_norm) ~
    heterogeneity + tau + scenario + variant,
  data = summary_tbl,
  FUN = mean
)
agg_sd <- aggregate(
  final_gap ~ heterogeneity + tau + scenario + variant,
  data = summary_tbl,
  FUN = sd
)
names(agg_sd)[names(agg_sd) == "final_gap"] <- "final_gap_sd"
agg_tbl <- merge(agg_tbl,
                 agg_sd,
                 by = c("heterogeneity", "tau", "scenario", "variant"),
                 all.x = TRUE)

base_tbl <- agg_tbl[agg_tbl$variant == "A0 operator",
                    c("heterogeneity", "tau", "scenario", "final_gap")]
names(base_tbl)[names(base_tbl) == "final_gap"] <- "base_final_gap"
agg_tbl <- merge(agg_tbl, base_tbl,
                 by = c("heterogeneity", "tau", "scenario"),
                 all.x = TRUE)
agg_tbl$gap_ratio_vs_base <- agg_tbl$final_gap / agg_tbl$base_final_gap

trace_agg <- aggregate(
  cbind(objective, gap, direction_norm) ~
    heterogeneity + tau + scenario + variant + round,
  data = trace_tbl,
  FUN = mean
)

write.csv(target_tbl,
          file.path(root, "results", "advanced_stabilization_targets.csv"),
          row.names = FALSE)
write.csv(summary_tbl,
          file.path(root, "results", "advanced_stabilization_summary.csv"),
          row.names = FALSE)
write.csv(agg_tbl,
          file.path(root, "results", "advanced_stabilization_aggregate.csv"),
          row.names = FALSE)
write.csv(trace_tbl,
          file.path(root, "results", "advanced_stabilization_trace.csv"),
          row.names = FALSE)
write.csv(trace_agg,
          file.path(root, "results", "advanced_stabilization_trace_aggregate.csv"),
          row.names = FALSE)

variant_cols <- setNames(
  c("#333333", "#2A6FBB", "#D95F02", "#1B9E77", "#7570B3",
    "#E7298A", "#66A61E", "#A6761D", "#E6AB02", "#666666"),
  variant_labels
)

plot_r1r2_gap <- function() {
  png(file.path(root, "figures", "advanced_stabilization_r1r2_gap.png"),
      width = 2200, height = 1400, res = 180)
  op <- par(mar = c(8.5, 5, 4, 1), mfrow = c(length(heterogeneity_levels), length(taus)))
  on.exit(par(op), add = TRUE)
  for (heterogeneity in heterogeneity_levels) {
    for (tau in taus) {
      d <- agg_tbl[
        agg_tbl$heterogeneity == heterogeneity &
          agg_tbl$tau == tau &
          agg_tbl$scenario == "R1 + R2",
      ]
      d <- d[match(variant_labels, d$variant), ]
      vals <- pmax(d$final_gap, 1e-8)
      names(vals) <- d$variant
      barplot(
        log10(vals),
        col = variant_cols[d$variant],
        las = 2,
        ylab = "log10(final gap)",
        main = sprintf("%s, tau=%.2f", heterogeneity, tau)
      )
      abline(h = log10(d$base_final_gap[1]), col = "#000000", lty = 3)
    }
  }
  dev.off()
}

plot_r1r2_ratio <- function() {
  png(file.path(root, "figures", "advanced_stabilization_r1r2_ratio.png"),
      width = 1800, height = 1200, res = 180)
  op <- par(mar = c(8.5, 5, 4, 1), mfrow = c(1, length(heterogeneity_levels)))
  on.exit(par(op), add = TRUE)
  for (heterogeneity in heterogeneity_levels) {
    d <- agg_tbl[agg_tbl$heterogeneity == heterogeneity &
                   agg_tbl$scenario == "R1 + R2" &
                   agg_tbl$variant != "A0 operator", ]
    d$key <- paste0(d$variant, "\n", "tau=", d$tau)
    vals <- d$gap_ratio_vs_base
    names(vals) <- d$key
    barplot(
      vals,
      col = variant_cols[d$variant],
      las = 2,
      ylab = "Final gap / A0 final gap",
      main = sprintf("R1 + R2 improvement ratio: %s", heterogeneity),
      ylim = c(0, max(vals, na.rm = TRUE) * 1.1)
    )
    abline(h = 1, col = "#000000", lty = 3)
  }
  dev.off()
}

plot_best_convergence <- function() {
  png(file.path(root, "figures", "advanced_stabilization_best_convergence.png"),
      width = 2200, height = 1400, res = 180)
  op <- par(mar = c(5, 5, 4, 1), mfrow = c(length(heterogeneity_levels), length(taus)))
  on.exit(par(op), add = TRUE)
  for (heterogeneity in heterogeneity_levels) {
    for (tau in taus) {
      d_rank <- agg_tbl[agg_tbl$heterogeneity == heterogeneity &
                          agg_tbl$tau == tau &
                          agg_tbl$scenario == "R1 + R2", ]
      keep <- unique(c(
        "A0 operator",
        as.character(d_rank$variant[order(d_rank$final_gap)][1:3])
      ))
      d <- trace_agg[trace_agg$heterogeneity == heterogeneity &
                       trace_agg$tau == tau &
                       trace_agg$scenario == "R1 + R2" &
                       trace_agg$variant %in% keep, ]
      y_lim <- quantile(pmax(d$gap, 0), 0.98, na.rm = TRUE)
      y_lim <- max(y_lim, 1e-4)
      plot(
        NA,
        xlim = range(d$round),
        ylim = c(0, y_lim),
        xlab = "Communication round",
        ylab = "Objective gap",
        main = sprintf("%s R1+R2, tau=%.2f", heterogeneity, tau)
      )
      for (variant in keep) {
        dm <- d[d$variant == variant, ]
        lines(dm$round, pmax(dm$gap, 0), col = variant_cols[[variant]], lwd = 2)
      }
      legend("topright", legend = keep, col = variant_cols[keep], lwd = 2, bty = "n", cex = 0.75)
    }
  }
  dev.off()
}

plot_r1r2_gap()
plot_r1r2_ratio()
plot_best_convergence()

cat("\nAdvanced stabilization aggregate, R1+R2 sorted by final gap:\n")
print(agg_tbl[agg_tbl$scenario == "R1 + R2", ][
  order(agg_tbl[agg_tbl$scenario == "R1 + R2", ]$heterogeneity,
        agg_tbl[agg_tbl$scenario == "R1 + R2", ]$tau,
        agg_tbl[agg_tbl$scenario == "R1 + R2", ]$final_gap),
])

cat("\nFigures written to:\n")
cat(file.path(root, "figures", "advanced_stabilization_r1r2_gap.png"), "\n")
cat(file.path(root, "figures", "advanced_stabilization_r1r2_ratio.png"), "\n")
cat(file.path(root, "figures", "advanced_stabilization_best_convergence.png"), "\n")
