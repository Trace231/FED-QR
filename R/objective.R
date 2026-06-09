check_loss <- function(u, tau) {
  ifelse(u >= 0, tau * u, (tau - 1) * u)
}

penalty_value <- function(beta, lambda = 0, penalty = c("none", "l1", "mcp", "scad"),
                          penalty_factor = NULL, mcp_gamma = 3, scad_a = 3.7) {
  penalty <- match.arg(penalty)
  if (is.null(penalty_factor)) {
    penalty_factor <- rep(1, length(beta))
  }

  abs_beta <- abs(beta)
  lambda_vec <- lambda * penalty_factor

  switch(
    penalty,
    none = 0,
    l1 = sum(lambda_vec * abs_beta),
    mcp = sum(ifelse(
      abs_beta <= mcp_gamma * lambda_vec,
      lambda_vec * abs_beta - abs_beta^2 / (2 * mcp_gamma),
      mcp_gamma * lambda_vec^2 / 2
    )),
    scad = {
      t <- abs_beta
      out <- numeric(length(t))
      r1 <- t <= lambda_vec
      r2 <- t > lambda_vec & t <= scad_a * lambda_vec
      r3 <- t > scad_a * lambda_vec
      out[r1] <- lambda_vec[r1] * t[r1]
      out[r2] <- (-t[r2]^2 + 2 * scad_a * lambda_vec[r2] * t[r2] -
        lambda_vec[r2]^2) / (2 * (scad_a - 1))
      out[r3] <- (scad_a + 1) * lambda_vec[r3]^2 / 2
      sum(out)
    }
  )
}

qr_objective <- function(X, y, beta, tau = 0.5, lambda = 0,
                         penalty = c("none", "l1", "mcp", "scad"),
                         penalty_factor = NULL, mcp_gamma = 3, scad_a = 3.7) {
  penalty <- match.arg(penalty)
  residual <- as.numeric(y - X %*% beta)
  loss <- mean(check_loss(residual, tau))

  if (is.null(penalty_factor)) {
    penalty_factor <- rep(1, length(beta))
  }

  loss + penalty_value(
    beta,
    lambda = lambda,
    penalty = penalty,
    penalty_factor = penalty_factor,
    mcp_gamma = mcp_gamma,
    scad_a = scad_a
  )
}

client_weight_vector <- function(client_sizes,
                                 client_weighting = c("sample", "uniform", "sqrt_size", "custom", "adaptive"),
                                 client_weights = NULL) {
  client_weighting <- match.arg(client_weighting)
  client_sizes <- as.numeric(client_sizes)
  stopifnot(length(client_sizes) > 0, all(client_sizes > 0))

  weights <- switch(
    client_weighting,
    sample = client_sizes / sum(client_sizes),
    uniform = rep(1 / length(client_sizes), length(client_sizes)),
    sqrt_size = sqrt(client_sizes) / sum(sqrt(client_sizes)),
    adaptive = {
      if (is.null(client_weights)) {
        client_sizes / sum(client_sizes)
      } else {
        client_weights <- as.numeric(client_weights)
        if (length(client_weights) != length(client_sizes)) {
          stop("client_weights must have one entry per client.", call. = FALSE)
        }
        if (any(!is.finite(client_weights)) || any(client_weights < 0) || sum(client_weights) <= 0) {
          stop("client_weights must be finite, nonnegative, and have positive sum.", call. = FALSE)
        }
        client_weights / sum(client_weights)
      }
    },
    custom = {
      if (is.null(client_weights)) {
        stop("client_weights must be supplied when client_weighting = 'custom'.", call. = FALSE)
      }
      client_weights <- as.numeric(client_weights)
      if (length(client_weights) != length(client_sizes)) {
        stop("client_weights must have one entry per client.", call. = FALSE)
      }
      if (any(!is.finite(client_weights)) || any(client_weights < 0) || sum(client_weights) <= 0) {
        stop("client_weights must be finite, nonnegative, and have positive sum.", call. = FALSE)
      }
      client_weights / sum(client_weights)
    }
  )
  weights
}

client_qr_objective <- function(X, y, beta, client_indices, tau = 0.5, lambda = 0,
                                penalty = c("none", "l1", "mcp", "scad"),
                                client_weighting = c("sample", "uniform", "sqrt_size", "custom", "adaptive"),
                                client_weights = NULL, penalty_factor = NULL,
                                mcp_gamma = 3, scad_a = 3.7) {
  penalty <- match.arg(penalty)
  X <- as.matrix(X)
  y <- as.numeric(y)
  client_indices <- lapply(client_indices, as.integer)
  sizes <- lengths(client_indices)
  weights <- client_weight_vector(sizes, client_weighting, client_weights)
  losses <- vapply(client_indices, function(idx) {
    residual <- as.numeric(y[idx] - X[idx, , drop = FALSE] %*% beta)
    mean(check_loss(residual, tau))
  }, numeric(1))
  if (is.null(penalty_factor)) {
    penalty_factor <- rep(1, length(beta))
  }
  sum(weights * losses) + penalty_value(
    beta,
    lambda = lambda,
    penalty = penalty,
    penalty_factor = penalty_factor,
    mcp_gamma = mcp_gamma,
    scad_a = scad_a
  )
}

client_loss_summary <- function(X, y, beta, client_indices, tau = 0.5) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  client_indices <- lapply(client_indices, as.integer)
  losses <- vapply(client_indices, function(idx) {
    residual <- as.numeric(y[idx] - X[idx, , drop = FALSE] %*% beta)
    mean(check_loss(residual, tau))
  }, numeric(1))
  data.frame(
    global_mean_loss = mean(check_loss(as.numeric(y - X %*% beta), tau)),
    client_mean_loss = mean(losses),
    worst_client_loss = max(losses),
    client_loss_sd = stats::sd(losses),
    client_q90_loss = as.numeric(stats::quantile(losses, 0.9)),
    min_client_loss = min(losses),
    stringsAsFactors = FALSE
  )
}

quantile_coverage <- function(X, y, beta, tau = 0.5, offsets = NULL, client_indices = NULL) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  pred <- as.numeric(X %*% beta)
  if (!is.null(offsets)) {
    if (is.null(client_indices)) {
      pred <- pred + as.numeric(offsets)
    } else {
      client_offsets <- as.numeric(offsets)
      for (j in seq_along(client_indices)) {
        pred[client_indices[[j]]] <- pred[client_indices[[j]]] + client_offsets[j]
      }
    }
  }
  mean(y <= pred)
}

intercept_column <- function(X) {
  if (!is.null(colnames(X)) && "(Intercept)" %in% colnames(X)) {
    return(which(colnames(X) == "(Intercept)")[1])
  }
  is_intercept <- vapply(seq_len(ncol(X)), function(j) {
    all(abs(X[, j] - 1) < 1e-12)
  }, logical(1))
  which(is_intercept)[1]
}

calibrate_quantile_intercept <- function(X, y, beta, tau = 0.5, client_indices = NULL,
                                         mode = c("global", "client_offset")) {
  mode <- match.arg(mode)
  X <- as.matrix(X)
  y <- as.numeric(y)
  beta <- as.numeric(beta)
  residual <- as.numeric(y - X %*% beta)

  if (mode == "global") {
    j0 <- intercept_column(X)
    if (is.na(j0)) {
      stop("Global intercept calibration requires an intercept column.", call. = FALSE)
    }
    delta <- as.numeric(stats::quantile(residual, tau, names = FALSE))
    beta_cal <- beta
    beta_cal[j0] <- beta_cal[j0] + delta
    return(list(beta = beta_cal, offset = delta, mode = mode))
  }

  if (is.null(client_indices)) {
    stop("client_indices must be supplied when mode = 'client_offset'.", call. = FALSE)
  }
  offsets <- vapply(client_indices, function(idx) {
    as.numeric(stats::quantile(residual[idx], tau, names = FALSE))
  }, numeric(1))
  list(beta = beta, offsets = offsets, mode = mode)
}

adaptive_calibrate_quantile <- function(X, y, beta, tau = 0.5, client_indices = NULL,
                                        metric = c("mean_client", "worst_client", "global")) {
  metric <- match.arg(metric)
  X <- as.matrix(X)
  y <- as.numeric(y)
  beta <- as.numeric(beta)

  raw <- calibration_summary(X, y, beta, tau = tau, client_indices = client_indices)
  candidates <- list(raw = list(beta = beta, offsets = NULL, selected_mode = "raw"))

  j0 <- intercept_column(X)
  if (!is.na(j0)) {
    global <- calibrate_quantile_intercept(X, y, beta, tau = tau, mode = "global")
    candidates$global_intercept <- list(
      beta = global$beta,
      offsets = NULL,
      selected_mode = "global_intercept"
    )
  }

  if (!is.null(client_indices)) {
    client <- calibrate_quantile_intercept(
      X, y, beta,
      tau = tau,
      client_indices = client_indices,
      mode = "client_offset"
    )
    candidates$client_offset <- list(
      beta = client$beta,
      offsets = client$offsets,
      selected_mode = "client_offset"
    )
  }

  scores <- vapply(candidates, function(candidate) {
    summary <- calibration_summary(
      X, y, candidate$beta,
      tau = tau,
      client_indices = client_indices,
      offsets = candidate$offsets
    )
    switch(
      metric,
      global = summary$global_coverage_error,
      mean_client = summary$mean_client_coverage_error,
      worst_client = summary$worst_client_coverage_error
    )
  }, numeric(1))

  best_name <- names(which.min(scores))[1]
  best <- candidates[[best_name]]
  best$mode <- "adaptive"
  best$metric <- metric
  best$score <- scores[[best_name]]
  best$scores <- scores
  best$raw_summary <- raw
  best
}

calibration_summary <- function(X, y, beta, tau = 0.5, client_indices = NULL,
                                offsets = NULL) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  pred <- as.numeric(X %*% beta)
  if (!is.null(offsets) && !is.null(client_indices)) {
    for (j in seq_along(client_indices)) {
      pred[client_indices[[j]]] <- pred[client_indices[[j]]] + offsets[j]
    }
  } else if (!is.null(offsets)) {
    pred <- pred + as.numeric(offsets)
  }
  global_cov <- mean(y <= pred)

  if (is.null(client_indices)) {
    client_cov <- global_cov
  } else {
    client_cov <- vapply(client_indices, function(idx) mean(y[idx] <= pred[idx]), numeric(1))
  }

  data.frame(
    global_coverage = global_cov,
    global_coverage_error = abs(global_cov - tau),
    mean_client_coverage_error = mean(abs(client_cov - tau)),
    worst_client_coverage_error = max(abs(client_cov - tau)),
    coverage_sd = stats::sd(client_cov),
    min_client_coverage = min(client_cov),
    max_client_coverage = max(client_cov),
    stringsAsFactors = FALSE
  )
}
