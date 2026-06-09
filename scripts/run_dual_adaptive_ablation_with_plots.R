args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
} else {
  file.path(getwd(), "scripts", "run_dual_adaptive_ablation_with_plots.R")
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

source(file.path(root, "R", "prox.R"))
source(file.path(root, "R", "objective.R"))
source(file.path(root, "R", "simulate.R"))
source(file.path(root, "R", "qr_pdhg.R"))
source(file.path(root, "R", "qr_box_fed_pdhg.R"))

dir.create(file.path(root, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "figures"), recursive = TRUE, showWarnings = FALSE)

step_labels <- c(
  operator = "A0 operator",
  box = "M1 box-aware",
  tau_adaptive = "M1+M2 tau-adaptive"
)

round_to_threshold <- function(trace, target, tol = 1e-3) {
  hit <- trace$round[trace$objective <= target + tol]
  if (length(hit) == 0) NA_integer_ else min(hit)
}

run_step_case <- function(dat, tau, target, scenario, K, batch_size,
                          step_rule, seed, rounds = 500) {
  fit <- qr_box_fed_pdhg(
    dat$X, dat$y,
    client_indices = dat$client_indices,
    tau = tau,
    rounds = rounds,
    clients_per_round = K,
    batch_size = batch_size,
    step_rule = step_rule,
    aggregation = "cached",
    beta_ref = target$beta,
    trace_every = 5,
    seed = seed
  )

  trace <- data.frame(
    tau = tau,
    scenario = scenario,
    seed = seed,
    step_rule = step_rule,
    step_label = unname(step_labels[[step_rule]]),
    round = fit$trace$round,
    objective = fit$trace$objective,
    gap = fit$trace$objective - target$objective,
    eta = fit$eta,
    sigma = fit$sigma
  )

  summary <- data.frame(
    tau = tau,
    scenario = scenario,
    seed = seed,
    step_rule = step_rule,
    step_label = unname(step_labels[[step_rule]]),
    eta = fit$eta,
    sigma = fit$sigma,
    final_objective = fit$objective,
    final_gap = fit$objective - target$objective,
    beta_l2_to_target = sqrt(sum((fit$beta - target$beta)^2)),
    round_to_gap_1e_3 = round_to_threshold(fit$trace, target$objective, 1e-3),
    dual_min = min(unlist(lapply(fit$clients, `[[`, "v"), use.names = FALSE)),
    dual_max = max(unlist(lapply(fit$clients, `[[`, "v"), use.names = FALSE))
  )

  list(summary = summary, trace = trace)
}

taus <- c(0.5, 0.75, 0.9, 0.95)
seeds <- c(20260526, 20260527, 20260528)
n_clients <- 20
n_per_client <- 80
p <- 15
rounds <- 500
heterogeneity <- "hard"
step_rules <- c("operator", "box", "tau_adaptive")

scenarios <- data.frame(
  scenario = c("Full", "R1 client", "R2 user", "R1 + R2"),
  K = c(n_clients, 4, n_clients, 4),
  batch_size = c(100000, 100000, 10, 10)
)

all_runs <- list()
target_rows <- list()
k <- 1
tr <- 1

for (seed in seeds) {
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
      tau = tau,
      seed = seed,
      target_objective = central$objective,
      target_beta_l2_to_truth = sqrt(sum((central$beta - dat$beta)^2))
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

      for (step_rule in step_rules) {
        message(sprintf(
          "seed=%d tau=%.2f scenario=%s step=%s",
          seed, tau, scenario, step_rule
        ))
        all_runs[[k]] <- run_step_case(
          dat = dat,
          tau = tau,
          target = target,
          scenario = scenario,
          K = K,
          batch_size = batch_size,
          step_rule = step_rule,
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
target_tbl <- do.call(rbind, target_rows)

agg_tbl <- aggregate(
  cbind(final_gap, beta_l2_to_target, eta, sigma) ~ tau + scenario + step_rule + step_label,
  data = summary_tbl,
  FUN = mean
)
agg_sd <- aggregate(
  cbind(final_gap, beta_l2_to_target) ~ tau + scenario + step_rule + step_label,
  data = summary_tbl,
  FUN = sd
)
names(agg_sd)[names(agg_sd) == "final_gap"] <- "final_gap_sd"
names(agg_sd)[names(agg_sd) == "beta_l2_to_target"] <- "beta_l2_to_target_sd"
agg_tbl <- merge(agg_tbl, agg_sd,
                 by = c("tau", "scenario", "step_rule", "step_label"),
                 all.x = TRUE)

hit_tbl <- aggregate(
  round_to_gap_1e_3 ~ tau + scenario + step_rule + step_label,
  data = summary_tbl,
  FUN = function(x) mean(x, na.rm = TRUE)
)
names(hit_tbl)[names(hit_tbl) == "round_to_gap_1e_3"] <- "mean_round_to_gap_1e_3"
agg_tbl <- merge(agg_tbl, hit_tbl,
                 by = c("tau", "scenario", "step_rule", "step_label"),
                 all.x = TRUE)

trace_agg <- aggregate(
  cbind(objective, gap) ~ tau + scenario + step_rule + step_label + round,
  data = trace_tbl,
  FUN = mean
)

write.csv(target_tbl,
          file.path(root, "results", "dual_adaptive_ablation_targets.csv"),
          row.names = FALSE)
write.csv(summary_tbl,
          file.path(root, "results", "dual_adaptive_ablation_summary.csv"),
          row.names = FALSE)
write.csv(agg_tbl,
          file.path(root, "results", "dual_adaptive_ablation_aggregate.csv"),
          row.names = FALSE)
write.csv(trace_tbl,
          file.path(root, "results", "dual_adaptive_ablation_trace.csv"),
          row.names = FALSE)
write.csv(trace_agg,
          file.path(root, "results", "dual_adaptive_ablation_trace_aggregate.csv"),
          row.names = FALSE)

step_cols <- c(
  operator = "#7F7F7F",
  box = "#2A6FBB",
  tau_adaptive = "#D95F02"
)

plot_final_gap_r1r2 <- function() {
  png(file.path(root, "figures", "dual_adaptive_ablation_r1r2_gap.png"),
      width = 1800, height = 1100, res = 180)
  op <- par(mar = c(7.5, 5, 4, 1), mfrow = c(1, length(taus)))
  on.exit(par(op), add = TRUE)
  for (tau in taus) {
    d <- agg_tbl[agg_tbl$tau == tau & agg_tbl$scenario == "R1 + R2", ]
    d <- d[match(step_rules, d$step_rule), ]
    vals <- pmax(d$final_gap, 1e-8)
    names(vals) <- d$step_label
    barplot(
      log10(vals),
      col = step_cols[d$step_rule],
      las = 2,
      ylab = "log10(final gap)",
      main = sprintf("R1 + R2, tau=%.2f", tau)
    )
    abline(h = -3, col = "#555555", lty = 3)
  }
  dev.off()
}

plot_convergence_r1r2 <- function() {
  png(file.path(root, "figures", "dual_adaptive_ablation_r1r2_convergence.png"),
      width = 1800, height = 1100, res = 180)
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
    for (step_rule in step_rules) {
      dm <- d[d$step_rule == step_rule, ]
      lines(dm$round, pmax(dm$gap, 0), col = step_cols[[step_rule]], lwd = 2)
    }
    legend("topright", legend = unname(step_labels[step_rules]),
           col = step_cols[step_rules], lwd = 2, bty = "n", cex = 0.75)
  }
  dev.off()
}

plot_gap_grid <- function() {
  png(file.path(root, "figures", "dual_adaptive_ablation_all_scenarios.png"),
      width = 2200, height = 1500, res = 180)
  op <- par(mar = c(6.5, 4.5, 3.2, 0.5), mfrow = c(length(taus), nrow(scenarios)))
  on.exit(par(op), add = TRUE)
  for (tau in taus) {
    for (scenario in scenarios$scenario) {
      d <- agg_tbl[agg_tbl$tau == tau & agg_tbl$scenario == scenario, ]
      d <- d[match(step_rules, d$step_rule), ]
      vals <- pmax(d$final_gap, 1e-8)
      names(vals) <- d$step_label
      barplot(
        log10(vals),
        col = step_cols[d$step_rule],
        las = 2,
        ylab = "log10(gap)",
        main = sprintf("tau=%.2f | %s", tau, scenario)
      )
      abline(h = -3, col = "#555555", lty = 3)
    }
  }
  dev.off()
}

plot_final_gap_r1r2()
plot_convergence_r1r2()
plot_gap_grid()

cat("\nDual adaptive ablation aggregate:\n")
print(agg_tbl[order(agg_tbl$tau, agg_tbl$scenario, agg_tbl$final_gap), ])

cat("\nFigures written to:\n")
cat(file.path(root, "figures", "dual_adaptive_ablation_r1r2_gap.png"), "\n")
cat(file.path(root, "figures", "dual_adaptive_ablation_r1r2_convergence.png"), "\n")
cat(file.path(root, "figures", "dual_adaptive_ablation_all_scenarios.png"), "\n")
