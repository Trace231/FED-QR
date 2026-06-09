qr_pdhg <- function(X, y, tau = 0.5, lambda = 0,
                    penalty = c("none", "l1", "mcp", "scad"),
                    max_iter = 2000, eta = NULL, sigma = NULL, theta = 1,
                    batch_size = NULL, step_rule = c("generic", "box", "tau_adaptive"),
                    intercept = FALSE, seed = 1, trace_every = 25,
                    m2_alpha = 0.5, m2_beta = 0.5, verbose = FALSE) {
  penalty <- match.arg(penalty)
  step_rule <- match.arg(step_rule)

  X <- as.matrix(X)
  y <- as.numeric(y)
  stopifnot(nrow(X) == length(y), tau > 0, tau < 1, max_iter >= 1)

  if (intercept) {
    X <- cbind("(Intercept)" = 1, X)
  }

  n <- nrow(X)
  p <- ncol(X)

  if (is.null(batch_size)) {
    batch_size <- n
  }
  batch_size <- min(batch_size, n)

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

  set.seed(seed)
  beta <- numeric(p)
  beta_bar <- beta
  v <- numeric(n)

  n_trace <- floor(max_iter / trace_every) + 1
  trace <- data.frame(iter = integer(n_trace), objective = numeric(n_trace))
  trace_pos <- 1
  trace[trace_pos, ] <- c(0, qr_objective(X, y, beta, tau, lambda, penalty, penalty_factor))

  lower <- tau - 1
  upper <- tau

  for (iter in seq_len(max_iter)) {
    idx <- if (batch_size == n) seq_len(n) else sample.int(n, batch_size)

    residual_bar <- as.numeric(y[idx] - X[idx, , drop = FALSE] %*% beta_bar)
    v[idx] <- pmin(upper, pmax(lower, v[idx] + sigma * residual_bar))

    primal_direction <- as.numeric(crossprod(X[idx, , drop = FALSE], v[idx]) / length(idx))

    beta_old <- beta
    beta <- prox_penalty(
      beta + eta * primal_direction,
      eta * lambda,
      penalty,
      penalty_factor,
      step_size = eta,
      lambda_value = lambda
    )
    beta_bar <- beta + theta * (beta - beta_old)

    if (iter %% trace_every == 0 || iter == max_iter) {
      trace_pos <- trace_pos + 1
      trace[trace_pos, ] <- c(iter, qr_objective(X, y, beta, tau, lambda, penalty, penalty_factor))
      if (verbose) {
        message(sprintf("iter=%d objective=%.6f", iter, trace[trace_pos, "objective"]))
      }
    }
  }

  trace <- trace[seq_len(trace_pos), , drop = FALSE]

  list(
    beta = beta,
    dual = v,
    objective = tail(trace$objective, 1),
    trace = trace,
    tau = tau,
    lambda = lambda,
    penalty = penalty,
    step_rule = step_rule,
    eta = eta,
    sigma = sigma
  )
}
