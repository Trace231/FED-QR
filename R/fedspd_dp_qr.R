smooth_check_loss <- function(u, tau, mu = 0.1) {
  tau * u + mu * log1p(exp(-u / mu))
}

smooth_qr_objective <- function(X, y, beta, tau = 0.5, mu = 0.1,
                                lambda = 0, penalty = c("none", "l1"),
                                penalty_factor = NULL) {
  penalty <- match.arg(penalty)
  residual <- as.numeric(y - X %*% beta)
  loss <- mean(smooth_check_loss(residual, tau, mu))

  if (is.null(penalty_factor)) {
    penalty_factor <- rep(1, length(beta))
  }

  pen <- switch(
    penalty,
    none = 0,
    l1 = sum(penalty_factor * abs(beta))
  )

  loss + lambda * pen
}

smooth_qr_gradient <- function(X, y, beta, tau = 0.5, mu = 0.1) {
  residual <- as.numeric(y - X %*% beta)
  psi <- tau - plogis(-residual / mu)
  -as.numeric(crossprod(X, psi) / nrow(X))
}

check_qr_subgradient <- function(X, y, beta, tau = 0.5) {
  residual <- as.numeric(y - X %*% beta)
  psi <- tau - as.numeric(residual < 0)
  -as.numeric(crossprod(X, psi) / nrow(X))
}

fedspd_gamma <- function(t, Q, batch_size, p_client, d, rho,
                         epsilon = 1, delta = 1e-4, G = 1,
                         phi = 1, d_lambda = 1, d_x = 1) {
  if (Q > 1) {
    c_q <- G^2 + 2 * d_lambda^2 + 2 * phi^2 / batch_size +
      16 * rho * d * G^2 * log(1.25 / delta) / ((Q - 1)^2 * epsilon^2)
  } else {
    c_q <- G^2 + 2 * d_lambda^2 + 2 * phi^2 / batch_size +
      16 * rho * d * G^2 * log(1.25 / delta) / (epsilon^2)
  }

  2 * sqrt(Q * p_client * c_q) * sqrt(t) / d_x
}

fedspd_dp_sigma <- function(Q, gamma_t, rho, epsilon = 1, delta = 1e-4, G = 1) {
  if (Q > 1) {
    4 * Q * G * sqrt(2 * log(1.25 / delta)) /
      ((Q - 1) * epsilon * (rho + gamma_t))
  } else {
    4 * G * sqrt(2 * log(1.25 / delta)) /
      (epsilon * (rho + gamma_t))
  }
}

fedspd_dp_qr <- function(X, y, client_indices = NULL, n_clients = 4,
                         tau = 0.5, smooth_mu = 0.1,
                         lambda = 0, penalty = c("none", "l1"),
                         loss = c("smooth", "check"),
                         rounds = 100, Q = 5, clients_per_round = NULL,
                         batch_size = 10, rho = 20,
                         gamma_rule = c("paper", "sqrt", "constant"),
                         gamma0 = NULL, dp = FALSE, epsilon = 1,
                         delta = 1e-4, G = 1, phi = 1,
                         d_lambda = 1, d_x = 1,
                         intercept = FALSE, beta_ref = NULL,
                         seed = 1, trace_every = 1,
                         sample_with_replacement = FALSE,
                         verbose = FALSE) {
  penalty <- match.arg(penalty)
  loss <- match.arg(loss)
  gamma_rule <- match.arg(gamma_rule)

  X <- as.matrix(X)
  y <- as.numeric(y)
  stopifnot(nrow(X) == length(y), tau > 0, tau < 1, rounds >= 1, Q >= 1)

  if (intercept) {
    X <- cbind("(Intercept)" = 1, X)
  }

  n <- nrow(X)
  d <- ncol(X)

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
  p_client <- clients_per_round / n_clients

  penalty_factor <- rep(1, d)
  if (!is.null(colnames(X)) && colnames(X)[1] == "(Intercept)") {
    penalty_factor[1] <- 0
  }

  if (is.null(gamma0)) {
    singular_max <- svd(X / sqrt(n), nu = 0, nv = 0)$d[1]
    gamma0 <- max(singular_max, .Machine$double.eps)
  }

  clients <- lapply(client_indices, function(idx) {
    list(
      X = X[idx, , drop = FALSE],
      y = y[idx],
      x = numeric(d),
      lambda = numeric(d),
      y_tilde = numeric(d),
      n = length(idx)
    )
  })

  x0 <- Reduce("+", lapply(clients, `[[`, "y_tilde")) / n_clients

  set.seed(seed)

  n_trace <- floor(rounds / trace_every) + 1
  trace <- data.frame(
    round = integer(n_trace),
    qr_objective = numeric(n_trace),
    smooth_objective = numeric(n_trace),
    consensus_gap = numeric(n_trace),
    beta_l2_error = numeric(n_trace),
    gamma = numeric(n_trace),
    dp_sigma = numeric(n_trace),
    selected_clients = integer(n_trace)
  )

  trace_pos <- 1
  trace[trace_pos, ] <- c(
    0,
    qr_objective(X, y, x0, tau, lambda, penalty, penalty_factor),
    smooth_qr_objective(X, y, x0, tau, smooth_mu, lambda, penalty, penalty_factor),
    mean(vapply(clients, function(cl) sqrt(sum((cl$x - x0)^2)), numeric(1))),
    if (is.null(beta_ref)) NA_real_ else sqrt(sum((x0 - beta_ref)^2)),
    NA_real_,
    0,
    0
  )

  for (t in seq_len(rounds)) {
    x0 <- Reduce("+", lapply(clients, `[[`, "y_tilde")) / n_clients

    gamma_t <- switch(
      gamma_rule,
      paper = fedspd_gamma(
        t, Q, batch_size, p_client, d, rho,
        epsilon = epsilon, delta = delta, G = G,
        phi = phi, d_lambda = d_lambda, d_x = d_x
      ),
      sqrt = gamma0 * sqrt(t),
      constant = gamma0
    )

    sigma_t <- if (dp) {
      fedspd_dp_sigma(Q, gamma_t, rho, epsilon = epsilon, delta = delta, G = G)
    } else {
      0
    }

    selected <- if (clients_per_round == n_clients) {
      seq_len(n_clients)
    } else {
      sample.int(n_clients, clients_per_round)
    }

    for (i in selected) {
      client <- clients[[i]]
      local_x <- client$x
      local_path <- matrix(NA_real_, nrow = Q, ncol = d)

      for (r in seq_len(Q)) {
        b <- min(batch_size, client$n)
        local_idx <- if (sample_with_replacement) {
          sample.int(client$n, b, replace = TRUE)
        } else if (b == client$n) {
          seq_len(client$n)
        } else {
          sample.int(client$n, b)
        }

        grad <- if (loss == "smooth") {
          smooth_qr_gradient(
            client$X[local_idx, , drop = FALSE],
            client$y[local_idx],
            local_x,
            tau = tau,
            mu = smooth_mu
          )
        } else {
          check_qr_subgradient(
            client$X[local_idx, , drop = FALSE],
            client$y[local_idx],
            local_x,
            tau = tau
          )
        }

        z <- (gamma_t * local_x + rho * x0 + client$lambda - grad) /
          (gamma_t + rho)
        local_x <- prox_penalty(
          z,
          lambda / (gamma_t + rho),
          penalty = penalty,
          penalty_factor = penalty_factor
        )
        local_path[r, ] <- local_x
      }

      client$x <- as.numeric(colMeans(local_path))
      client$lambda <- client$lambda + rho * (x0 - client$x)
      noise <- if (dp) rnorm(d, mean = 0, sd = sigma_t) else numeric(d)
      client$y_tilde <- client$x - client$lambda / rho + noise
      clients[[i]] <- client
    }

    if (t %% trace_every == 0 || t == rounds) {
      x0_trace <- Reduce("+", lapply(clients, `[[`, "y_tilde")) / n_clients
      trace_pos <- trace_pos + 1
      trace[trace_pos, ] <- c(
        t,
        qr_objective(X, y, x0_trace, tau, lambda, penalty, penalty_factor),
        smooth_qr_objective(X, y, x0_trace, tau, smooth_mu, lambda, penalty, penalty_factor),
        mean(vapply(clients, function(cl) sqrt(sum((cl$x - x0_trace)^2)), numeric(1))),
        if (is.null(beta_ref)) NA_real_ else sqrt(sum((x0_trace - beta_ref)^2)),
        gamma_t,
        sigma_t,
        length(selected)
      )
      if (verbose) {
        message(sprintf(
          "round=%d qr_objective=%.6f consensus_gap=%.6f",
          t,
          trace[trace_pos, "qr_objective"],
          trace[trace_pos, "consensus_gap"]
        ))
      }
    }
  }

  x0 <- Reduce("+", lapply(clients, `[[`, "y_tilde")) / n_clients
  trace <- trace[seq_len(trace_pos), , drop = FALSE]

  list(
    x0 = x0,
    clients = clients,
    beta = x0,
    qr_objective = tail(trace$qr_objective, 1),
    smooth_objective = tail(trace$smooth_objective, 1),
    consensus_gap = tail(trace$consensus_gap, 1),
    trace = trace,
    tau = tau,
    smooth_mu = smooth_mu,
    lambda = lambda,
    penalty = penalty,
    loss = loss,
    rounds = rounds,
    Q = Q,
    clients_per_round = clients_per_round,
    batch_size = batch_size,
    rho = rho,
    gamma_rule = gamma_rule,
    dp = dp,
    epsilon = epsilon,
    delta = delta
  )
}
