args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
} else {
  file.path(getwd(), "scripts", "run_r1_r2_stress_test.R")
}
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

source(file.path(root, "R", "prox.R"))
source(file.path(root, "R", "objective.R"))
source(file.path(root, "R", "simulate.R"))
source(file.path(root, "R", "qr_pdhg.R"))
source(file.path(root, "R", "fedspd_dp_qr.R"))
source(file.path(root, "R", "qr_box_fed_pdhg.R"))

dir.create(file.path(root, "results"), recursive = TRUE, showWarnings = FALSE)

partition_diagnostics <- function(dat, client_indices, partition_label, alpha) {
  do.call(rbind, lapply(seq_along(client_indices), function(j) {
    idx <- client_indices[[j]]
    data.frame(
      partition = partition_label,
      alpha = alpha,
      client = j,
      n = length(idx),
      y_mean = mean(dat$y[idx]),
      y_sd = sd(dat$y[idx]),
      y_q10 = as.numeric(quantile(dat$y[idx], 0.10)),
      y_q50 = as.numeric(quantile(dat$y[idx], 0.50)),
      y_q90 = as.numeric(quantile(dat$y[idx], 0.90))
    )
  }))
}

round_to_threshold <- function(trace, target, tol = 1e-3) {
  hit <- trace$round[trace$objective <= target + tol]
  if (length(hit) == 0) NA_integer_ else min(hit)
}

summarize_run <- function(trace, target, checkpoints = c(25, 50, 100, 250, 500)) {
  vals <- sapply(checkpoints, function(cp) {
    idx <- max(which(trace$round <= cp))
    trace$objective[idx]
  })
  names(vals) <- paste0("obj_round_", checkpoints)
  c(vals, round_to_target_1e_3 = round_to_threshold(trace, target, 1e-3))
}

run_box_case <- function(dat, client_indices, tau, scenario, alpha, K, batch_size,
                         rounds, target_obj, seed) {
  fit <- qr_box_fed_pdhg(
    dat$X, dat$y,
    client_indices = client_indices,
    tau = tau,
    rounds = rounds,
    clients_per_round = K,
    batch_size = batch_size,
    step_rule = "operator",
    aggregation = "cached",
    beta_ref = dat$beta,
    trace_every = 5,
    seed = seed
  )

  trace <- fit$trace
  trace$method <- "qr_box_fed_pdhg"
  trace$scenario <- scenario
  trace$alpha <- alpha
  trace$tau <- tau
  trace$K <- K
  trace$batch_size <- batch_size
  trace$seed <- seed

  extra <- summarize_run(trace, target_obj)
  summary <- data.frame(
    method = "qr_box_fed_pdhg",
    scenario = scenario,
    alpha = alpha,
    tau = tau,
    K = K,
    batch_size = batch_size,
    seed = seed,
    rounds = rounds,
    final_objective = fit$objective,
    final_gap = fit$objective - target_obj,
    beta_l2_error = sqrt(sum((fit$beta - dat$beta)^2)),
    dual_min = min(unlist(lapply(fit$clients, `[[`, "v"), use.names = FALSE)),
    dual_max = max(unlist(lapply(fit$clients, `[[`, "v"), use.names = FALSE)),
    t(extra),
    check.names = FALSE
  )

  list(summary = summary, trace = trace)
}

run_fedspd_check_case <- function(dat, client_indices, tau, scenario, alpha, K,
                                  batch_size, rounds, target_obj, seed) {
  fit <- fedspd_dp_qr(
    dat$X, dat$y,
    client_indices = client_indices,
    tau = tau,
    loss = "check",
    rounds = min(rounds, 500),
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
    method = "fedspd_check",
    scenario = scenario,
    alpha = alpha,
    tau = tau,
    K = K,
    batch_size = batch_size,
    seed = seed
  )

  extra <- summarize_run(trace, target_obj)
  summary <- data.frame(
    method = "fedspd_check",
    scenario = scenario,
    alpha = alpha,
    tau = tau,
    K = K,
    batch_size = batch_size,
    seed = seed,
    rounds = min(rounds, 500),
    final_objective = fit$qr_objective,
    final_gap = fit$qr_objective - target_obj,
    beta_l2_error = sqrt(sum((fit$beta - dat$beta)^2)),
    dual_min = NA_real_,
    dual_max = NA_real_,
    t(extra),
    check.names = FALSE
  )

  list(summary = summary, trace = trace)
}

tau <- 0.9
n_clients <- 20
n <- 1200
p <- 15
rounds <- 500
seeds <- c(20260526, 20260527, 20260528)
alphas <- c(10, 1, 0.25, 0.05)

scenarios <- data.frame(
  scenario = c("deterministic", "R1_client_sampling", "R2_user_minibatch", "R1_plus_R2"),
  K = c(n_clients, 4, n_clients, 4),
  batch_size = c(100000, 100000, 10, 10)
)

all_runs <- list()
all_partition <- list()
central_rows <- list()
k <- 1
pd <- 1
cr <- 1

for (seed in seeds) {
  dat <- make_qr_sim(n = n, p = p, tau = tau, noise = "asymmetric", seed = seed)
  central <- qr_pdhg(
    dat$X, dat$y,
    tau = tau,
    max_iter = 3000,
    step_rule = "generic",
    trace_every = 100,
    seed = seed
  )
  target_obj <- central$objective

  central_rows[[cr]] <- data.frame(
    method = "central_qr_pdhg",
    scenario = "central",
    alpha = NA_real_,
    tau = tau,
    K = NA_integer_,
    batch_size = NA_integer_,
    seed = seed,
    rounds = 3000,
    final_objective = target_obj,
    final_gap = 0,
    beta_l2_error = sqrt(sum((central$beta - dat$beta)^2)),
    dual_min = NA_real_,
    dual_max = NA_real_,
    obj_round_25 = NA_real_,
    obj_round_50 = NA_real_,
    obj_round_100 = NA_real_,
    obj_round_250 = NA_real_,
    obj_round_500 = NA_real_,
    round_to_target_1e_3 = NA_integer_
  )
  cr <- cr + 1

  for (alpha in alphas) {
    partition_label <- paste0("dirichlet_alpha_", alpha)
    client_indices <- dirichlet_partition(dat$y, n_clients = n_clients, alpha = alpha, seed = seed)
    full_batch <- max(lengths(client_indices))

    all_partition[[pd]] <- partition_diagnostics(dat, client_indices, partition_label, alpha)
    pd <- pd + 1

    for (s in seq_len(nrow(scenarios))) {
      scenario <- scenarios$scenario[s]
      K <- scenarios$K[s]
      batch_size <- if (scenarios$batch_size[s] >= 100000) full_batch else scenarios$batch_size[s]

      message(sprintf(
        "seed=%d alpha=%.2f scenario=%s box",
        seed, alpha, scenario
      ))
      all_runs[[k]] <- run_box_case(
        dat, client_indices, tau,
        scenario = scenario,
        alpha = alpha,
        K = K,
        batch_size = batch_size,
        rounds = rounds,
        target_obj = target_obj,
        seed = seed
      )
      k <- k + 1

      if (scenario %in% c("deterministic", "R1_plus_R2")) {
        message(sprintf(
          "seed=%d alpha=%.2f scenario=%s fedspd_check",
          seed, alpha, scenario
        ))
        all_runs[[k]] <- run_fedspd_check_case(
          dat, client_indices, tau,
          scenario = scenario,
          alpha = alpha,
          K = K,
          batch_size = batch_size,
          rounds = rounds,
          target_obj = target_obj,
          seed = seed
        )
        k <- k + 1
      }
    }
  }
}

summary_tbl <- do.call(rbind, lapply(all_runs, `[[`, "summary"))
trace_tbl <- do.call(rbind, lapply(all_runs, `[[`, "trace"))
partition_tbl <- do.call(rbind, all_partition)
summary_tbl <- rbind(summary_tbl, do.call(rbind, central_rows))

aggregate_tbl <- aggregate(
  cbind(
    final_gap,
    beta_l2_error,
    obj_round_25,
    obj_round_50,
    obj_round_100,
    obj_round_250,
    obj_round_500
  ) ~ method + scenario + alpha + tau + K + batch_size,
  data = subset(summary_tbl, method != "central_qr_pdhg"),
  FUN = mean,
  na.rm = TRUE
)

write.csv(summary_tbl,
          file.path(root, "results", "r1_r2_stress_summary.csv"),
          row.names = FALSE)
write.csv(trace_tbl,
          file.path(root, "results", "r1_r2_stress_trace.csv"),
          row.names = FALSE)
write.csv(partition_tbl,
          file.path(root, "results", "r1_r2_partition_diagnostics.csv"),
          row.names = FALSE)
write.csv(aggregate_tbl,
          file.path(root, "results", "r1_r2_stress_aggregate.csv"),
          row.names = FALSE)

cat("\nR1/R2 stress aggregate:\n")
print(aggregate_tbl)

