prox_check_loss <- function(z, kappa, tau) {
  z - pmin(kappa * tau, pmax(kappa * (tau - 1), z))
}

lasso_cd <- function(X, y, lambda, penalty_factor = NULL, beta_init = NULL,
                     max_iter = 200, tol = 1e-8) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  p <- ncol(X)

  if (is.null(penalty_factor)) {
    penalty_factor <- rep(1, p)
  }
  if (is.null(beta_init)) {
    beta <- numeric(p)
  } else {
    beta <- as.numeric(beta_init)
  }

  x_norm2 <- colSums(X^2)
  x_norm2[x_norm2 == 0] <- 1
  residual <- y - as.numeric(X %*% beta)

  for (iter in seq_len(max_iter)) {
    beta_old <- beta
    for (j in seq_len(p)) {
      residual <- residual + X[, j] * beta[j]
      z <- sum(X[, j] * residual)
      beta[j] <- soft_threshold(z, lambda * penalty_factor[j]) / x_norm2[j]
      if (penalty_factor[j] == 0) {
        beta[j] <- z / x_norm2[j]
      }
      residual <- residual - X[, j] * beta[j]
    }
    if (sqrt(sum((beta - beta_old)^2)) < tol * (1 + sqrt(sum(beta_old^2)))) {
      break
    }
  }

  beta
}

solve_beta_ls <- function(X, target, ridge = 0, rhs_extra = NULL, beta_init = NULL,
                          lambda = 0, penalty = c("none", "l1"),
                          penalty_factor = NULL, cd_max_iter = 200) {
  penalty <- match.arg(penalty)
  X <- as.matrix(X)
  target <- as.numeric(target)
  p <- ncol(X)

  if (is.null(rhs_extra)) {
    rhs_extra <- numeric(p)
  }

  if (penalty == "l1" && lambda > 0) {
    if (ridge > 0) {
      X_aug <- rbind(X, sqrt(ridge) * diag(p))
      y_aug <- c(target, rhs_extra / sqrt(ridge))
    } else {
      X_aug <- X
      y_aug <- target
    }
    return(lasso_cd(
      X_aug, y_aug,
      lambda = lambda,
      penalty_factor = penalty_factor,
      beta_init = beta_init,
      max_iter = cd_max_iter
    ))
  }

  lhs <- crossprod(X)
  if (ridge > 0) {
    lhs <- lhs + ridge * diag(p)
  }
  rhs <- as.numeric(crossprod(X, target) + rhs_extra)

  chol2inv_safe <- tryCatch(chol(lhs), error = function(e) NULL)
  if (is.null(chol2inv_safe)) {
    lhs <- lhs + 1e-8 * diag(p)
    chol2inv_safe <- chol(lhs)
  }
  as.numeric(backsolve(chol2inv_safe, forwardsolve(t(chol2inv_safe), rhs)))
}

qr_admm <- function(X, y, tau = 0.5, lambda = 0,
                    penalty = c("none", "l1"), rho = 0.3,
                    max_iter = 2000, intercept = FALSE,
                    beta_init = NULL, trace_every = 10,
                    abstol = 1e-5, reltol = 1e-4,
                    cd_max_iter = 200, check_convergence = FALSE,
                    verbose = FALSE) {
  penalty <- match.arg(penalty)
  X <- as.matrix(X)
  y <- as.numeric(y)
  stopifnot(nrow(X) == length(y), tau > 0, tau < 1, rho > 0)

  if (intercept) {
    X <- cbind("(Intercept)" = 1, X)
  }

  n <- nrow(X)
  p <- ncol(X)

  penalty_factor <- rep(1, p)
  if (!is.null(colnames(X)) && colnames(X)[1] == "(Intercept)") {
    penalty_factor[1] <- 0
  }

  beta <- if (is.null(beta_init)) numeric(p) else as.numeric(beta_init)
  r <- y - as.numeric(X %*% beta)
  u <- numeric(n)

  n_trace <- floor(max_iter / trace_every) + 1
  trace <- data.frame(
    iter = integer(n_trace),
    objective = numeric(n_trace),
    primal_residual = numeric(n_trace),
    dual_residual = numeric(n_trace)
  )
  trace_pos <- 1
  trace[trace_pos, ] <- c(
    0,
    qr_objective(X, y, beta, tau, lambda, penalty, penalty_factor),
    sqrt(sum((y - as.numeric(X %*% beta) - r)^2)),
    NA_real_
  )

  for (iter in seq_len(max_iter)) {
    r_old <- r
    target <- y - r + u
    beta <- solve_beta_ls(
      X, target,
      lambda = lambda / rho,
      penalty = penalty,
      penalty_factor = penalty_factor,
      beta_init = beta,
      cd_max_iter = cd_max_iter
    )

    z <- y - as.numeric(X %*% beta) + u
    r <- prox_check_loss(z, kappa = 1 / (n * rho), tau = tau)
    primal_res <- y - as.numeric(X %*% beta) - r
    u <- u + primal_res

    dual_res <- rho * sqrt(sum((r - r_old)^2))

    if (iter %% trace_every == 0 || iter == max_iter) {
      trace_pos <- trace_pos + 1
      trace[trace_pos, ] <- c(
        iter,
        qr_objective(X, y, beta, tau, lambda, penalty, penalty_factor),
        sqrt(sum(primal_res^2)),
        dual_res
      )
      if (verbose) {
        message(sprintf(
          "iter=%d objective=%.6f primal=%.3e dual=%.3e",
          iter,
          trace[trace_pos, "objective"],
          trace[trace_pos, "primal_residual"],
          trace[trace_pos, "dual_residual"]
        ))
      }
    }

    if (check_convergence) {
      eps_pri <- sqrt(n) * abstol + reltol * max(
        sqrt(sum((y - as.numeric(X %*% beta))^2)),
        sqrt(sum(r^2))
      )
      eps_dual <- sqrt(p) * abstol + reltol * sqrt(sum(as.numeric(crossprod(X, rho * u))^2))
      if (sqrt(sum(primal_res^2)) < eps_pri && dual_res < eps_dual) {
        break
      }
    }
  }

  trace <- trace[seq_len(trace_pos), , drop = FALSE]

  list(
    beta = beta,
    residual = r,
    dual = u,
    objective = tail(trace$objective, 1),
    trace = trace,
    tau = tau,
    lambda = lambda,
    penalty = penalty,
    rho = rho,
    iterations = iter
  )
}

local_qr_consensus_admm <- function(X, y, z, u_consensus, tau = 0.5,
                                    rho_consensus = 10, rho_residual = 1,
                                    n_total = nrow(X), beta_init = NULL,
                                    inner_iter = 50, batch_idx = NULL) {
  X_full <- as.matrix(X)
  y_full <- as.numeric(y)
  p <- ncol(X_full)

  if (is.null(batch_idx)) {
    Xb <- X_full
    yb <- y_full
  } else {
    Xb <- X_full[batch_idx, , drop = FALSE]
    yb <- y_full[batch_idx]
  }

  beta <- if (is.null(beta_init)) z else as.numeric(beta_init)
  r <- yb - as.numeric(Xb %*% beta)
  u_residual <- numeric(length(yb))

  lhs <- rho_residual * crossprod(Xb) + rho_consensus * diag(p)
  chol_lhs <- tryCatch(chol(lhs), error = function(e) NULL)
  if (is.null(chol_lhs)) {
    lhs <- lhs + 1e-8 * diag(p)
    chol_lhs <- chol(lhs)
  }

  for (inner in seq_len(inner_iter)) {
    rhs <- rho_residual * as.numeric(crossprod(Xb, yb - r + u_residual)) +
      rho_consensus * (z - u_consensus)
    beta <- as.numeric(backsolve(chol_lhs, forwardsolve(t(chol_lhs), rhs)))

    q <- yb - as.numeric(Xb %*% beta) + u_residual
    r <- prox_check_loss(q, kappa = 1 / (n_total * rho_residual), tau = tau)
    u_residual <- u_residual + yb - as.numeric(Xb %*% beta) - r
  }

  beta
}

fed_qr_admm <- function(X, y, client_indices = NULL, n_clients = 4,
                        tau = 0.5, lambda = 0, penalty = c("none", "l1"),
                        rounds = 500, clients_per_round = NULL,
                        batch_size = NULL, rho_consensus = 10,
                        rho_residual = 1, inner_iter = 30,
                        intercept = FALSE, beta_ref = NULL,
                        seed = 1, trace_every = 5,
                        verbose = FALSE) {
  penalty <- match.arg(penalty)
  X <- as.matrix(X)
  y <- as.numeric(y)
  stopifnot(nrow(X) == length(y), tau > 0, tau < 1)

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

  clients <- lapply(client_indices, function(idx) {
    list(
      X = X[idx, , drop = FALSE],
      y = y[idx],
      beta = numeric(p),
      u = numeric(p),
      n = length(idx)
    )
  })

  z <- numeric(p)
  set.seed(seed)

  n_trace <- floor(rounds / trace_every) + 1
  trace <- data.frame(
    round = integer(n_trace),
    objective = numeric(n_trace),
    consensus_gap = numeric(n_trace),
    beta_l2_error = numeric(n_trace),
    selected_clients = integer(n_trace)
  )
  trace_pos <- 1
  trace[trace_pos, ] <- c(
    0,
    qr_objective(X, y, z, tau, lambda, penalty, penalty_factor),
    mean(vapply(clients, function(cl) sqrt(sum((cl$beta - z)^2)), numeric(1))),
    if (is.null(beta_ref)) NA_real_ else sqrt(sum((z - beta_ref)^2)),
    0
  )

  for (round in seq_len(rounds)) {
    selected <- if (clients_per_round == n_clients) {
      seq_len(n_clients)
    } else {
      sample.int(n_clients, clients_per_round)
    }

    for (j in selected) {
      client <- clients[[j]]
      b <- min(batch_size, client$n)
      batch_idx <- if (b == client$n) seq_len(client$n) else sample.int(client$n, b)
      client$beta <- local_qr_consensus_admm(
        client$X, client$y,
        z = z,
        u_consensus = client$u,
        tau = tau,
        rho_consensus = rho_consensus,
        rho_residual = rho_residual,
        n_total = n,
        beta_init = client$beta,
        inner_iter = inner_iter,
        batch_idx = batch_idx
      )
      clients[[j]] <- client
    }

    z_bar <- Reduce("+", lapply(clients, function(cl) cl$beta + cl$u)) / n_clients
    z <- prox_penalty(
      z_bar,
      gamma = lambda / (n_clients * rho_consensus),
      penalty = penalty,
      penalty_factor = penalty_factor
    )

    for (j in selected) {
      clients[[j]]$u <- clients[[j]]$u + clients[[j]]$beta - z
    }

    if (round %% trace_every == 0 || round == rounds) {
      trace_pos <- trace_pos + 1
      trace[trace_pos, ] <- c(
        round,
        qr_objective(X, y, z, tau, lambda, penalty, penalty_factor),
        mean(vapply(clients, function(cl) sqrt(sum((cl$beta - z)^2)), numeric(1))),
        if (is.null(beta_ref)) NA_real_ else sqrt(sum((z - beta_ref)^2)),
        length(selected)
      )
      if (verbose) {
        message(sprintf(
          "round=%d objective=%.6f gap=%.3e",
          round,
          trace[trace_pos, "objective"],
          trace[trace_pos, "consensus_gap"]
        ))
      }
    }
  }

  trace <- trace[seq_len(trace_pos), , drop = FALSE]

  list(
    beta = z,
    clients = clients,
    objective = tail(trace$objective, 1),
    trace = trace,
    tau = tau,
    lambda = lambda,
    penalty = penalty,
    rho_consensus = rho_consensus,
    rho_residual = rho_residual,
    inner_iter = inner_iter,
    clients_per_round = clients_per_round,
    batch_size = batch_size
  )
}
