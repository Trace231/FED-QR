qr_box_fed_pdhg <- function(X, y, client_indices = NULL, n_clients = 4,
                            tau = 0.5, lambda = 0,
                            penalty = c("none", "l1", "mcp", "scad"),
                            rounds = 1000, clients_per_round = NULL,
                            batch_size = NULL, eta = NULL, sigma = NULL,
                            theta = 1,
                            step_rule = c("operator", "box", "tau_adaptive"),
                            aggregation = c("cached", "selected_reweighted"),
                            intercept = FALSE, beta_ref = NULL,
                            seed = 1, trace_every = 10,
                            m2_alpha = 0.5, m2_beta = 0.5,
                            dual_relaxation = 1,
                            server_momentum = 0,
                            step_decay_power = 0,
                            step_decay_offset = 100,
                            primal_clip = NULL,
                            staleness = c("none", "exponential", "inverse", "adaptive"),
                            staleness_rate = 0.03,
                            staleness_floor = 0.25,
                            staleness_normalize = TRUE,
                            adaptive_staleness_target = NULL,
                            client_weighting = c("sample", "uniform", "sqrt_size", "custom", "adaptive"),
                            client_weights = NULL,
                            adaptive_client_power = 1,
                            adaptive_client_blend = 0.75,
                            adaptive_client_floor = 0.05,
                            adaptive_client_smooth = 0.5,
                            verbose = FALSE) {
  penalty <- match.arg(penalty)
  step_rule <- match.arg(step_rule)
  aggregation <- match.arg(aggregation)
  staleness <- match.arg(staleness)
  client_weighting <- match.arg(client_weighting)

  X <- as.matrix(X)
  y <- as.numeric(y)
  stopifnot(
    nrow(X) == length(y), tau > 0, tau < 1, rounds >= 1,
    dual_relaxation > 0, dual_relaxation <= 1,
    server_momentum >= 0, server_momentum < 1,
    step_decay_power >= 0, step_decay_offset > 0,
    staleness_rate >= 0, staleness_floor > 0, staleness_floor <= 1,
    adaptive_client_power >= 0,
    adaptive_client_blend >= 0, adaptive_client_blend <= 1,
    adaptive_client_floor >= 0, adaptive_client_floor < 1,
    adaptive_client_smooth >= 0, adaptive_client_smooth <= 1
  )

  if (intercept) {
    X <- cbind("(Intercept)" = 1, X)
  }

  n <- nrow(X)
  p <- ncol(X)

  if (is.null(client_indices)) {
    client_indices <- iid_partition(n, n_clients = n_clients, seed = seed)
  }
  client_indices <- lapply(client_indices, as.integer)
  client_indices <- client_indices[lengths(client_indices) > 0]
  n_clients <- length(client_indices)
  stopifnot(n_clients >= 1)

  if (is.null(clients_per_round)) {
    clients_per_round <- n_clients
  }
  clients_per_round <- min(clients_per_round, n_clients)

  client_sizes <- lengths(client_indices)
  if (is.null(batch_size)) {
    batch_size <- max(client_sizes)
  }
  fixed_client_weights <- client_weight_vector(
    client_sizes,
    client_weighting = if (client_weighting == "adaptive") "sample" else client_weighting,
    client_weights = client_weights
  )
  client_weights <- fixed_client_weights
  if (is.null(adaptive_staleness_target)) {
    participation_rate <- clients_per_round / n_clients
    adaptive_staleness_target <- max(
      1,
      (1 - participation_rate) / max(participation_rate, .Machine$double.eps)
    )
  }

  penalty_factor <- rep(1, p)
  if (!is.null(colnames(X)) && colnames(X)[1] == "(Intercept)") {
    penalty_factor[1] <- 0
  }

  op_norm <- svd(X / sqrt(n), nu = 0, nv = 0)$d[1]
  op_norm <- max(op_norm, .Machine$double.eps)
  base_step <- 0.95 / op_norm

  vmax <- max(tau, 1 - tau)
  tau_width <- max(0.05, 1 - 2 * abs(tau - 0.5))

  if (is.null(eta) || is.null(sigma)) {
    eta0 <- base_step
    sigma0 <- base_step

    if (step_rule %in% c("box", "tau_adaptive")) {
      eta0 <- eta0 / vmax
      sigma0 <- sigma0 * vmax
    }

    if (step_rule == "tau_adaptive") {
      eta0 <- eta0 * tau_width^(-m2_beta)
      sigma0 <- sigma0 * tau_width^(m2_alpha)
    }

    if (is.null(eta)) eta <- eta0
    if (is.null(sigma)) sigma <- sigma0
  }

  lower <- tau - 1
  upper <- tau

  clients <- lapply(client_indices, function(idx) {
    list(
      X = X[idx, , drop = FALSE],
      y = y[idx],
      v = numeric(length(idx)),
      direction = numeric(p),
      raw_direction = numeric(p),
      n = length(idx)
    )
  })

  set.seed(seed)
  beta <- numeric(p)
  beta_bar <- beta
  direction_ema <- numeric(p)
  client_age <- integer(n_clients)
  stale_weight <- rep(1, n_clients)

  n_trace <- floor(rounds / trace_every) + 1
  trace <- data.frame(
    round = integer(n_trace),
    objective = numeric(n_trace),
    beta_l2_error = numeric(n_trace),
    selected_clients = integer(n_trace),
    mean_selected_n = numeric(n_trace),
    dual_min = numeric(n_trace),
    dual_max = numeric(n_trace),
    eta = numeric(n_trace),
    sigma = numeric(n_trace),
    direction_norm = numeric(n_trace),
    client_objective = numeric(n_trace),
    mean_staleness = numeric(n_trace),
    max_staleness = numeric(n_trace),
    mean_stale_weight = numeric(n_trace),
    adaptive_staleness_rate = numeric(n_trace),
    adaptive_client_weight_sd = numeric(n_trace),
    max_client_weight = numeric(n_trace)
  )

  trace_pos <- 1
  trace[trace_pos, ] <- c(
    0,
    qr_objective(X, y, beta, tau, lambda, penalty, penalty_factor),
    if (is.null(beta_ref)) NA_real_ else sqrt(sum((beta - beta_ref)^2)),
    0,
    0,
    0,
    0,
    eta,
    sigma,
    0,
    client_qr_objective(
      X, y, beta, client_indices,
      tau = tau, lambda = lambda, penalty = penalty,
      client_weighting = if (client_weighting == "adaptive") "custom" else client_weighting,
      client_weights = client_weights,
      penalty_factor = penalty_factor
    ),
    0,
    0,
    1,
    staleness_rate,
    stats::sd(client_weights),
    max(client_weights)
  )

  for (round in seq_len(rounds)) {
    step_scale <- (1 + round / step_decay_offset)^(-step_decay_power)
    eta_t <- eta * step_scale
    sigma_t <- sigma * step_scale

    selected <- if (clients_per_round == n_clients) {
      seq_len(n_clients)
    } else {
      sample.int(n_clients, clients_per_round)
    }
    selected_n <- client_sizes[selected]

    if (client_weighting == "adaptive") {
      local_losses <- vapply(clients, function(client) {
        residual <- as.numeric(client$y - client$X %*% beta_bar)
        mean(check_loss(residual, tau))
      }, numeric(1))
      loss_center <- sum(fixed_client_weights * local_losses)
      if (is.finite(loss_center) && loss_center > .Machine$double.eps) {
        relative_loss <- pmax(local_losses / loss_center, .Machine$double.eps)
        raw_weights <- fixed_client_weights * relative_loss^adaptive_client_power
        raw_weights <- raw_weights / sum(raw_weights)
        floored_weights <- adaptive_client_floor / n_clients +
          (1 - adaptive_client_floor) * raw_weights
        target_weights <- (1 - adaptive_client_blend) * fixed_client_weights +
          adaptive_client_blend * floored_weights
        client_weights <- (1 - adaptive_client_smooth) * client_weights +
          adaptive_client_smooth * target_weights
        client_weights <- client_weights / sum(client_weights)
      }
    }

    selected_direction <- numeric(p)
    for (j in selected) {
      client <- clients[[j]]
      b <- min(batch_size, client$n)
      local_idx <- if (b == client$n) seq_len(client$n) else sample.int(client$n, b)

      residual <- as.numeric(client$y[local_idx] -
        client$X[local_idx, , drop = FALSE] %*% beta_bar)
      v_candidate <- pmin(upper, pmax(lower, client$v[local_idx] + sigma_t * residual))
      client$v[local_idx] <- (1 - dual_relaxation) * client$v[local_idx] +
        dual_relaxation * v_candidate

      client$raw_direction <- as.numeric(crossprod(client$X, client$v) / client$n)
      client$direction <- client_weights[j] * client$raw_direction
      selected_direction <- selected_direction +
        (client_weights[j] / sum(client_weights[selected])) *
          client$raw_direction
      clients[[j]] <- client
    }

    client_age <- client_age + 1L
    client_age[selected] <- 0L
    adaptive_rate_t <- staleness_rate
    if (staleness == "adaptive") {
      age_pressure <- mean(client_age) / max(adaptive_staleness_target, .Machine$double.eps)
      adaptive_rate_t <- staleness_rate * pmax(0.25, age_pressure)
    }
    stale_weight <- switch(
      staleness,
      none = rep(1, n_clients),
      exponential = pmax(staleness_floor, exp(-staleness_rate * client_age)),
      inverse = pmax(staleness_floor, 1 / (1 + staleness_rate * client_age)),
      adaptive = pmax(staleness_floor, exp(-adaptive_rate_t * client_age))
    )
    aggregation_weight <- client_weights * stale_weight
    if (staleness != "none" && isTRUE(staleness_normalize)) {
      denom <- sum(aggregation_weight)
      if (denom > .Machine$double.eps) {
        aggregation_weight <- aggregation_weight / denom
      }
    }

    primal_direction <- if (aggregation == "cached") {
      Reduce("+", Map(function(client, weight) {
        weight * client$raw_direction
      }, clients, aggregation_weight))
    } else {
      selected_direction
    }

    if (server_momentum > 0) {
      direction_ema <- server_momentum * direction_ema +
        (1 - server_momentum) * primal_direction
      primal_direction <- direction_ema
    }

    direction_norm <- sqrt(sum(primal_direction^2))
    if (!is.null(primal_clip) && is.finite(primal_clip) &&
        primal_clip > 0 && direction_norm > primal_clip) {
      primal_direction <- primal_direction * (primal_clip / direction_norm)
      direction_norm <- primal_clip
    }

    beta_old <- beta
    beta <- prox_penalty(
      beta + eta_t * primal_direction,
      eta_t * lambda,
      penalty,
      penalty_factor,
      step_size = eta_t,
      lambda_value = lambda
    )
    beta_bar <- beta + theta * (beta - beta_old)

    if (round %% trace_every == 0 || round == rounds) {
      all_v <- unlist(lapply(clients, `[[`, "v"), use.names = FALSE)
      trace_pos <- trace_pos + 1
      trace[trace_pos, ] <- c(
        round,
        qr_objective(X, y, beta, tau, lambda, penalty, penalty_factor),
        if (is.null(beta_ref)) NA_real_ else sqrt(sum((beta - beta_ref)^2)),
        length(selected),
        mean(selected_n),
        min(all_v),
        max(all_v),
        eta_t,
        sigma_t,
        direction_norm,
          client_qr_objective(
            X, y, beta, client_indices,
            tau = tau, lambda = lambda, penalty = penalty,
            client_weighting = if (client_weighting == "adaptive") "custom" else client_weighting,
            client_weights = client_weights,
            penalty_factor = penalty_factor
          ),
          mean(client_age),
          max(client_age),
          mean(stale_weight),
          adaptive_rate_t,
          stats::sd(client_weights),
          max(client_weights)
        )
      if (verbose) {
        message(sprintf(
          "round=%d objective=%.6f dual=[%.3f, %.3f]",
          round,
          trace[trace_pos, "objective"],
          trace[trace_pos, "dual_min"],
          trace[trace_pos, "dual_max"]
        ))
      }
    }
  }

  trace <- trace[seq_len(trace_pos), , drop = FALSE]

  list(
    beta = beta,
    clients = clients,
    objective = tail(trace$objective, 1),
    trace = trace,
    tau = tau,
    lambda = lambda,
    penalty = penalty,
    eta = eta,
    sigma = sigma,
    theta = theta,
    step_rule = step_rule,
    aggregation = aggregation,
    clients_per_round = clients_per_round,
    batch_size = batch_size,
    client_sizes = client_sizes,
    client_weights = client_weights,
    fixed_client_weights = fixed_client_weights,
    client_weighting = client_weighting,
    dual_box = c(lower, upper),
    dual_relaxation = dual_relaxation,
    server_momentum = server_momentum,
    step_decay_power = step_decay_power,
    step_decay_offset = step_decay_offset,
    primal_clip = primal_clip,
    staleness = staleness,
    staleness_rate = staleness_rate,
    staleness_floor = staleness_floor,
    staleness_normalize = staleness_normalize,
    adaptive_staleness_target = adaptive_staleness_target,
    adaptive_client_power = adaptive_client_power,
    adaptive_client_blend = adaptive_client_blend,
    adaptive_client_floor = adaptive_client_floor,
    adaptive_client_smooth = adaptive_client_smooth
  )
}
