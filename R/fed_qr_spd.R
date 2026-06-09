fed_qr_spd <- function(X, y, client_indices = NULL, n_clients = 4,
                       tau = 0.5, lambda = 0, penalty = c("none", "l1"),
                       rounds = 1000, clients_per_round = NULL,
                       batch_size = NULL, eta = NULL, sigma = NULL, theta = 1,
                       step_rule = c("generic", "box", "tau_adaptive"),
                       direction_mode = c("batch", "memory"),
                       client_weighting = c("renormalized", "unbiased", "global"),
                       intercept = FALSE, beta_true = NULL, seed = 1,
                       trace_every = 10, m2_alpha = 0.5, m2_beta = 0.5,
                       verbose = FALSE) {
  penalty <- match.arg(penalty)
  step_rule <- match.arg(step_rule)
  direction_mode <- match.arg(direction_mode)
  client_weighting <- match.arg(client_weighting)

  X <- as.matrix(X)
  y <- as.numeric(y)
  stopifnot(nrow(X) == length(y), tau > 0, tau < 1, rounds >= 1)

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

  penalty_factor <- rep(1, p)
  if (!is.null(colnames(X)) && colnames(X)[1] == "(Intercept)") {
    penalty_factor[1] <- 0
  }

  singular_max <- svd(X / sqrt(n), nu = 0, nv = 0)$d[1]
  singular_max <- max(singular_max, .Machine$double.eps)
  base_step <- 0.9 / singular_max

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

  clients <- lapply(client_indices, function(idx) {
    list(
      X = X[idx, , drop = FALSE],
      y = y[idx],
      v = numeric(length(idx)),
      n = length(idx)
    )
  })

  set.seed(seed)
  beta <- numeric(p)
  beta_bar <- beta
  lower <- tau - 1
  upper <- tau

  n_trace <- floor(rounds / trace_every) + 1
  trace <- data.frame(
    round = integer(n_trace),
    objective = numeric(n_trace),
    beta_l2_error = numeric(n_trace),
    selected_clients = integer(n_trace),
    mean_selected_n = numeric(n_trace)
  )
  trace_pos <- 1
  trace[trace_pos, ] <- c(
    0,
    qr_objective(X, y, beta, tau, lambda, penalty, penalty_factor),
    if (is.null(beta_true)) NA_real_ else sqrt(sum((beta - beta_true)^2)),
    0,
    0
  )

  for (round in seq_len(rounds)) {
    selected <- if (clients_per_round == n_clients) {
      seq_len(n_clients)
    } else {
      sample.int(n_clients, clients_per_round)
    }

    selected_n <- client_sizes[selected]
    direction <- numeric(p)

    for (j in selected) {
      client <- clients[[j]]
      b <- min(batch_size, client$n)
      local_idx <- if (b == client$n) seq_len(client$n) else sample.int(client$n, b)

      residual_bar <- as.numeric(client$y[local_idx] -
        client$X[local_idx, , drop = FALSE] %*% beta_bar)
      client$v[local_idx] <- pmin(upper, pmax(lower, client$v[local_idx] + sigma * residual_bar))

      local_direction <- if (direction_mode == "memory") {
        as.numeric(crossprod(client$X, client$v) / client$n)
      } else {
        as.numeric(crossprod(
          client$X[local_idx, , drop = FALSE],
          client$v[local_idx]
        ) / length(local_idx))
      }

      weight <- switch(
        client_weighting,
        renormalized = client$n / sum(selected_n),
        unbiased = (n_clients / clients_per_round) * client$n / n,
        global = client$n / n
      )

      direction <- direction + weight * local_direction
      clients[[j]] <- client
    }

    beta_old <- beta
    beta <- prox_penalty(beta + eta * direction, eta * lambda, penalty, penalty_factor)
    beta_bar <- beta + theta * (beta - beta_old)

    if (round %% trace_every == 0 || round == rounds) {
      trace_pos <- trace_pos + 1
      trace[trace_pos, ] <- c(
        round,
        qr_objective(X, y, beta, tau, lambda, penalty, penalty_factor),
        if (is.null(beta_true)) NA_real_ else sqrt(sum((beta - beta_true)^2)),
        length(selected),
        mean(selected_n)
      )
      if (verbose) {
        message(sprintf(
          "round=%d objective=%.6f",
          round,
          trace[trace_pos, "objective"]
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
    step_rule = step_rule,
    direction_mode = direction_mode,
    client_weighting = client_weighting,
    eta = eta,
    sigma = sigma,
    clients_per_round = clients_per_round,
    batch_size = batch_size,
    client_sizes = client_sizes
  )
}
