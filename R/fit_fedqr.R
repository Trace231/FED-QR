#' List registered federated quantile-regression methods
#'
#' @return Character vector of method names accepted by [fit_fedqr()].
#' @export
fedqr_methods <- function() {
  c(
    "QR box-dual",
    "QR box-dual long",
    "QR box-dual stale",
    "QR box-dual robust",
    "QR box-dual stale+robust",
    "QR box-dual adaptive",
    "FSPG-smooth",
    "FedSubGrad",
    "FedSPD-check",
    "FedSPD-smooth",
    "FedQR-ADMM"
  )
}

match_fedqr_method <- function(method) {
  aliases <- c(
    "qr_box_dual" = "QR box-dual",
    "qr-box-dual" = "QR box-dual",
    "QR box-dual" = "QR box-dual",
    "qr_box_dual_long" = "QR box-dual long",
    "qr-box-dual-long" = "QR box-dual long",
    "QR box-dual long" = "QR box-dual long",
    "qr_box_dual_stale" = "QR box-dual stale",
    "qr-box-dual-stale" = "QR box-dual stale",
    "QR box-dual stale" = "QR box-dual stale",
    "qr_box_dual_robust" = "QR box-dual robust",
    "qr-box-dual-robust" = "QR box-dual robust",
    "QR box-dual robust" = "QR box-dual robust",
    "qr_box_dual_stale_robust" = "QR box-dual stale+robust",
    "qr-box-dual-stale-robust" = "QR box-dual stale+robust",
    "QR box-dual stale+robust" = "QR box-dual stale+robust",
    "qr_box_dual_adaptive" = "QR box-dual adaptive",
    "qr-box-dual-adaptive" = "QR box-dual adaptive",
    "QR box-dual adaptive" = "QR box-dual adaptive",
    "fspg" = "FSPG-smooth",
    "fspg_smooth" = "FSPG-smooth",
    "FSPG-smooth" = "FSPG-smooth",
    "fed_subgrad" = "FedSubGrad",
    "FedSubGrad" = "FedSubGrad",
    "fedspd_check" = "FedSPD-check",
    "FedSPD-check" = "FedSPD-check",
    "fedspd_smooth" = "FedSPD-smooth",
    "FedSPD-smooth" = "FedSPD-smooth",
    "fedqr_admm" = "FedQR-ADMM",
    "FedQR-ADMM" = "FedQR-ADMM"
  )
  method <- as.character(method)[1]
  if (!method %in% names(aliases)) {
    stop("Unknown method: ", method, call. = FALSE)
  }
  unname(aliases[[method]])
}

merge_control <- function(defaults, control) {
  if (is.null(control)) {
    return(defaults)
  }
  modifyList(defaults, control)
}

normalize_fedqr_trace <- function(trace, method) {
  if (is.null(trace)) {
    return(data.frame(round = integer(), objective = numeric(), method = character()))
  }
  trace <- as.data.frame(trace)
  if (!"round" %in% names(trace) && "iter" %in% names(trace)) {
    names(trace)[names(trace) == "iter"] <- "round"
  }
  if (!"objective" %in% names(trace) && "qr_objective" %in% names(trace)) {
    trace$objective <- trace$qr_objective
  }
  if (!"round" %in% names(trace)) {
    trace$round <- seq_len(nrow(trace)) - 1
  }
  if (!"objective" %in% names(trace)) {
    stop("Trace for method ", method, " does not contain an objective column.", call. = FALSE)
  }
  trace$method <- method
  first <- c("round", "objective", "method")
  trace[, c(first, setdiff(names(trace), first)), drop = FALSE]
}

rbind_fill <- function(tables) {
  cols <- unique(unlist(lapply(tables, names), use.names = FALSE))
  filled <- lapply(tables, function(tbl) {
    missing <- setdiff(cols, names(tbl))
    for (col in missing) {
      tbl[[col]] <- NA
    }
    tbl[, cols, drop = FALSE]
  })
  do.call(rbind, filled)
}

#' Fit one federated quantile-regression method
#'
#' This is the stable public wrapper used by experiment scripts. It preserves
#' all lower-level algorithm implementations and normalizes their output shape.
#'
#' @param method Method name. See [fedqr_methods()].
#' @param X Numeric design matrix.
#' @param y Numeric response vector.
#' @param client_indices Optional list of row indices for each client.
#' @param tau Quantile level in `(0, 1)`.
#' @param lambda Penalty level.
#' @param penalty Penalty type passed to lower-level solvers.
#' @param rounds Communication rounds for federated methods.
#' @param clients_per_round Number of active clients per round.
#' @param batch_size Local mini-batch size.
#' @param seed Random seed.
#' @param trace_every Trace interval.
#' @param control Named list of method-specific overrides.
#' @return List with normalized `beta`, `objective`, `trace`, and raw `fit`.
#' @export
fit_fedqr <- function(method, X, y, client_indices = NULL,
                      tau = 0.5, lambda = 0, penalty = "none",
                      rounds = 1000, clients_per_round = NULL,
                      batch_size = NULL, seed = 1, trace_every = 10,
                      control = list()) {
  method <- match_fedqr_method(method)
  X <- as.matrix(X)
  y <- as.numeric(y)
  stopifnot(nrow(X) == length(y), tau > 0, tau < 1, rounds >= 1)

  base <- list(
    X = X,
    y = y,
    client_indices = client_indices,
    tau = tau,
    lambda = lambda,
    penalty = penalty,
    rounds = rounds,
    clients_per_round = clients_per_round,
    batch_size = batch_size,
    seed = seed,
    trace_every = trace_every
  )

  fit <- switch(
    method,
    "QR box-dual" = {
      args <- merge_control(c(base, list(step_rule = "operator", aggregation = "cached")), control)
      do.call(qr_box_fed_pdhg, args)
    },
    "QR box-dual long" = {
      long_multiplier <- if (is.null(control$long_multiplier)) 2 else control$long_multiplier
      control$long_multiplier <- NULL
      base$rounds <- as.integer(rounds * long_multiplier)
      args <- merge_control(c(base, list(step_rule = "operator", aggregation = "cached")), control)
      do.call(qr_box_fed_pdhg, args)
    },
    "QR box-dual stale" = {
      args <- merge_control(c(base, list(
        step_rule = "operator",
        aggregation = "cached",
        staleness = "exponential"
      )), control)
      do.call(qr_box_fed_pdhg, args)
    },
    "QR box-dual robust" = {
      args <- merge_control(c(base, list(
        step_rule = "operator",
        aggregation = "cached",
        client_weighting = "uniform"
      )), control)
      do.call(qr_box_fed_pdhg, args)
    },
    "QR box-dual stale+robust" = {
      args <- merge_control(c(base, list(
        step_rule = "operator",
        aggregation = "cached",
        staleness = "exponential",
        client_weighting = "uniform"
      )), control)
      do.call(qr_box_fed_pdhg, args)
    },
    "QR box-dual adaptive" = {
      args <- merge_control(c(base, list(
        step_rule = "operator",
        aggregation = "cached",
        staleness = "adaptive",
        client_weighting = "adaptive"
      )), control)
      do.call(qr_box_fed_pdhg, args)
    },
    "FSPG-smooth" = {
      args <- merge_control(c(base, list(smooth_mu = 0.05)), control)
      do.call(fed_smooth_qr, args)
    },
    "FedSubGrad" = {
      args <- merge_control(base, control)
      do.call(fed_subgrad_qr, args)
    },
    "FedSPD-check" = {
      args <- merge_control(c(base, list(
        loss = "check",
        Q = 3,
        rho = 20,
        gamma_rule = "sqrt",
        gamma0 = 4,
        dp = FALSE
      )), control)
      do.call(fedspd_dp_qr, args)
    },
    "FedSPD-smooth" = {
      args <- merge_control(c(base, list(
        loss = "smooth",
        smooth_mu = 0.1,
        Q = 3,
        rho = 20,
        gamma_rule = "sqrt",
        gamma0 = 4,
        dp = FALSE
      )), control)
      do.call(fedspd_dp_qr, args)
    },
    "FedQR-ADMM" = {
      args <- merge_control(c(base, list(
        rho_consensus = 10,
        rho_residual = 1,
        inner_iter = 30
      )), control)
      do.call(fed_qr_admm, args)
    }
  )

  objective <- if (!is.null(fit$objective)) {
    fit$objective
  } else if (!is.null(fit$qr_objective)) {
    fit$qr_objective
  } else {
    tail(normalize_fedqr_trace(fit$trace, method)$objective, 1)
  }

  list(
    method = method,
    beta = as.numeric(fit$beta),
    objective = as.numeric(objective),
    trace = normalize_fedqr_trace(fit$trace, method),
    fit = fit
  )
}

#' Fit several federated quantile-regression methods
#'
#' @param methods Character vector of method names.
#' @param method_controls Optional named list of per-method control lists.
#' @inheritParams fit_fedqr
#' @param term_names Optional coefficient names.
#' @return List with `summary`, `trace`, `coefficients`, and raw `fits`.
#' @export
run_fedqr_methods <- function(methods, X, y, client_indices = NULL,
                              tau = 0.5, lambda = 0, penalty = "none",
                              rounds = 1000, clients_per_round = NULL,
                              batch_size = NULL, seed = 1,
                              trace_every = 10, control = list(),
                              method_controls = list(),
                              term_names = colnames(X),
                              verbose = TRUE) {
  methods <- vapply(methods, match_fedqr_method, character(1))
  fits <- lapply(methods, function(method) {
    if (isTRUE(verbose)) {
      message(sprintf("Running method=%s", method))
    }
    method_control <- method_controls[[method]]
    if (is.null(method_control)) {
      method_control <- method_controls[[gsub("[ -]", "_", method)]]
    }
    fit_fedqr(
      method = method,
      X = X,
      y = y,
      client_indices = client_indices,
      tau = tau,
      lambda = lambda,
      penalty = penalty,
      rounds = rounds,
      clients_per_round = clients_per_round,
      batch_size = batch_size,
      seed = seed,
      trace_every = trace_every,
      control = merge_control(control, method_control)
    )
  })
  names(fits) <- methods

  summary <- do.call(rbind, lapply(fits, function(fit) {
    data.frame(
      method = fit$method,
      objective = fit$objective,
      beta_norm = sqrt(sum(fit$beta^2))
    )
  }))
  rownames(summary) <- NULL
  best_obj <- min(summary$objective)
  summary$gap_to_best_observed <- summary$objective - best_obj

  trace <- rbind_fill(lapply(fits, `[[`, "trace"))
  trace$gap_to_best_observed <- trace$objective - best_obj

  betas <- do.call(cbind, lapply(fits, `[[`, "beta"))
  colnames(betas) <- methods
  if (is.null(term_names)) {
    term_names <- paste0("beta_", seq_len(nrow(betas)))
  }
  coefficients <- data.frame(term = term_names, betas, check.names = FALSE)

  result <- list(
    summary = summary,
    trace = trace,
    coefficients = coefficients,
    fits = fits
  )
  validate_result_table(result$summary)
  result
}
