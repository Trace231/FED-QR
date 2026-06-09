args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
} else {
  file.path(getwd(), "scripts", "run_hard_noniid_with_plots.R")
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

source(file.path(root, "R", "prox.R"))
source(file.path(root, "R", "objective.R"))
source(file.path(root, "R", "simulate.R"))
source(file.path(root, "R", "qr_pdhg.R"))
source(file.path(root, "R", "fedspd_dp_qr.R"))
source(file.path(root, "R", "qr_box_fed_pdhg.R"))

dir.create(file.path(root, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "figures"), recursive = TRUE, showWarnings = FALSE)

run_box <- function(dat, scenario, K, batch_size, target_obj, seed, rounds = 600) {
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
    method = "QR box-dual",
    scenario = scenario,
    heterogeneity = dat$heterogeneity,
    seed = seed,
    gap = objective - target_obj
  )

  summary <- data.frame(
    method = "QR box-dual",
    scenario = scenario,
    heterogeneity = dat$heterogeneity,
    seed = seed,
    final_objective = fit$objective,
    final_gap = fit$objective - target_obj,
    beta_l2_error = sqrt(sum((fit$beta - dat$beta)^2)),
    dual_min = min(unlist(lapply(fit$clients, `[[`, "v"), use.names = FALSE)),
    dual_max = max(unlist(lapply(fit$clients, `[[`, "v"), use.names = FALSE))
  )

  list(summary = summary, trace = trace)
}

run_fedspd <- function(dat, scenario, K, batch_size, target_obj, seed, rounds = 600) {
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
    method = "FedSPD-check",
    scenario = scenario,
    heterogeneity = dat$heterogeneity,
    seed = seed,
    gap = fit$trace$qr_objective - target_obj
  )

  summary <- data.frame(
    method = "FedSPD-check",
    scenario = scenario,
    heterogeneity = dat$heterogeneity,
    seed = seed,
    final_objective = fit$qr_objective,
    final_gap = fit$qr_objective - target_obj,
    beta_l2_error = sqrt(sum((fit$beta - dat$beta)^2)),
    dual_min = NA_real_,
    dual_max = NA_real_
  )

  list(summary = summary, trace = trace)
}

seeds <- c(20260526, 20260527, 20260528)
heterogeneity_levels <- c("mild", "hard", "extreme")
n_clients <- 20

scenarios <- data.frame(
  scenario = c("Deterministic", "R1 client sampling", "R2 user mini-batch", "R1 + R2"),
  K = c(n_clients, 4, n_clients, 4),
  batch_size = c(100000, 100000, 10, 10)
)

all_runs <- list()
client_info <- list()
central_rows <- list()
k <- 1
ci <- 1
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
    dat$client_info$seed <- seed
    client_info[[ci]] <- dat$client_info
    ci <- ci + 1

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
      beta_l2_error = sqrt(sum((central$beta - dat$beta)^2)),
      dual_min = NA_real_,
      dual_max = NA_real_
    )
    cr <- cr + 1

    full_batch <- max(lengths(dat$client_indices))
    for (s in seq_len(nrow(scenarios))) {
      scenario <- scenarios$scenario[s]
      K <- scenarios$K[s]
      batch_size <- if (scenarios$batch_size[s] >= 100000) full_batch else scenarios$batch_size[s]

      message(sprintf("seed=%d hetero=%s scenario=%s box", seed, heterogeneity, scenario))
      all_runs[[k]] <- run_box(dat, scenario, K, batch_size, target_obj, seed)
      k <- k + 1

      if (scenario %in% c("Deterministic", "R1 + R2")) {
        message(sprintf("seed=%d hetero=%s scenario=%s fedspd", seed, heterogeneity, scenario))
        all_runs[[k]] <- run_fedspd(dat, scenario, K, batch_size, target_obj, seed)
        k <- k + 1
      }
    }
  }
}

summary_tbl <- do.call(rbind, lapply(all_runs, `[[`, "summary"))
trace_tbl <- do.call(rbind, lapply(all_runs, `[[`, "trace"))
client_tbl <- do.call(rbind, client_info)
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
          file.path(root, "results", "hard_noniid_summary.csv"),
          row.names = FALSE)
write.csv(agg_tbl,
          file.path(root, "results", "hard_noniid_aggregate.csv"),
          row.names = FALSE)
write.csv(trace_tbl,
          file.path(root, "results", "hard_noniid_trace.csv"),
          row.names = FALSE)
write.csv(trace_agg,
          file.path(root, "results", "hard_noniid_trace_aggregate.csv"),
          row.names = FALSE)
write.csv(client_tbl,
          file.path(root, "results", "hard_noniid_client_info.csv"),
          row.names = FALSE)

plot_final_gap <- function() {
  png(file.path(root, "figures", "hard_noniid_final_gap.png"),
      width = 1800, height = 1100, res = 180)
  op <- par(mar = c(8, 5, 4, 1), mfrow = c(1, 3))
  on.exit(par(op), add = TRUE)
  for (h in heterogeneity_levels) {
    d <- subset(agg_tbl, heterogeneity == h & method != "Central QR-PDHG")
    scenarios_order <- unique(d$scenario)
    methods <- unique(d$method)
    mat <- sapply(methods, function(m) {
      vals <- d$final_gap[match(paste(m, scenarios_order), paste(d$method, d$scenario))]
      vals
    })
    mat <- t(mat)
    barplot(
      mat,
      beside = TRUE,
      names.arg = scenarios_order,
      las = 2,
      col = c("#2A6FBB", "#D95F02")[seq_len(nrow(mat))],
      ylab = "Final objective gap",
      main = paste("Heterogeneity:", h),
      ylim = c(0, max(mat, na.rm = TRUE) * 1.2)
    )
    legend("topright", legend = rownames(mat), fill = c("#2A6FBB", "#D95F02"), bty = "n", cex = 0.8)
  }
  dev.off()
}

plot_r1r2_curve <- function() {
  png(file.path(root, "figures", "hard_noniid_r1r2_convergence.png"),
      width = 1800, height = 1100, res = 180)
  op <- par(mar = c(5, 5, 4, 1), mfrow = c(1, 3))
  on.exit(par(op), add = TRUE)
  for (h in heterogeneity_levels) {
    d <- subset(trace_agg, heterogeneity == h & scenario == "R1 + R2")
    plot(
      NA,
      xlim = range(d$round),
      ylim = range(pmax(d$gap, 0), na.rm = TRUE),
      xlab = "Communication round",
      ylab = "Objective gap",
      main = paste("R1 + R2:", h)
    )
    cols <- c("QR box-dual" = "#2A6FBB", "FedSPD-check" = "#D95F02")
    for (m in names(cols)) {
      dm <- subset(d, method == m)
      lines(dm$round, pmax(dm$gap, 0), col = cols[[m]], lwd = 2)
    }
    legend("topright", legend = names(cols), col = cols, lwd = 2, bty = "n", cex = 0.8)
  }
  dev.off()
}

plot_client_heterogeneity <- function() {
  png(file.path(root, "figures", "hard_noniid_client_heterogeneity.png"),
      width = 1800, height = 1100, res = 180)
  op <- par(mar = c(5, 5, 4, 1), mfrow = c(1, 3))
  on.exit(par(op), add = TRUE)
  for (h in heterogeneity_levels) {
    d <- subset(client_tbl, heterogeneity == h & seed == seeds[1])
    plot(
      d$client,
      d$y_q90,
      type = "b",
      pch = 19,
      col = "#2A6FBB",
      xlab = "Client",
      ylab = "Client response quantiles",
      main = paste("Client shifts:", h),
      ylim = range(c(d$y_q50, d$y_q90))
    )
    lines(d$client, d$y_q50, type = "b", pch = 17, col = "#D95F02")
    legend("topleft", legend = c("y q90", "y median"), col = c("#2A6FBB", "#D95F02"),
           pch = c(19, 17), lwd = 1, bty = "n")
  }
  dev.off()
}

plot_final_gap()
plot_r1r2_curve()
plot_client_heterogeneity()

cat("\nHard non-IID aggregate:\n")
print(agg_tbl)

cat("\nFigures written to:\n")
cat(file.path(root, "figures", "hard_noniid_final_gap.png"), "\n")
cat(file.path(root, "figures", "hard_noniid_r1r2_convergence.png"), "\n")
cat(file.path(root, "figures", "hard_noniid_client_heterogeneity.png"), "\n")

