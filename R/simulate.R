make_qr_sim <- function(n = 500, p = 10, tau = 0.5, signal = 1,
                        noise = c("normal", "t", "asymmetric"),
                        seed = 1, intercept = TRUE) {
  noise <- match.arg(noise)
  set.seed(seed)

  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  beta <- signal * c(1, -0.8, 0.6, -0.4, rep(0, max(0, p - 4)))
  beta <- beta[seq_len(p)]

  eps <- switch(
    noise,
    normal = rnorm(n),
    t = rt(n, df = 3),
    asymmetric = rexp(n) - log(2)
  )

  eps_q <- switch(
    noise,
    normal = qnorm(tau),
    t = qt(tau, df = 3),
    asymmetric = qexp(tau) - log(2)
  )

  y <- as.numeric(X %*% beta + eps)

  if (intercept) {
    X <- cbind("(Intercept)" = 1, X)
    beta <- c("(Intercept)" = eps_q, beta)
  }

  colnames(X) <- make.names(colnames(X), unique = TRUE)

  list(X = X, y = y, beta = beta, tau = tau)
}

iid_partition <- function(n, n_clients = 4, seed = 1) {
  set.seed(seed)
  shuffled <- sample.int(n)
  split(shuffled, rep(seq_len(n_clients), length.out = n))
}

dirichlet_partition <- function(y, n_clients = 4, alpha = 10, seed = 1) {
  set.seed(seed)
  n <- length(y)
  ranks <- rank(y, ties.method = "random")
  strata <- cut(ranks, breaks = quantile(ranks, probs = seq(0, 1, length.out = n_clients + 1)),
                include.lowest = TRUE, labels = FALSE)

  client_ids <- integer(n)
  for (s in sort(unique(strata))) {
    idx <- which(strata == s)
    weights <- rgamma(n_clients, shape = alpha, rate = 1)
    probs <- weights / sum(weights)
    client_ids[idx] <- sample(seq_len(n_clients), length(idx), replace = TRUE, prob = probs)
  }

  out <- split(seq_len(n), client_ids)
  out[as.character(seq_len(n_clients))]
}

make_hard_federated_qr_sim <- function(n_clients = 20, n_per_client = 80, p = 15,
                                       tau = 0.9, heterogeneity = c("mild", "hard", "extreme"),
                                       sparsity = min(5, p),
                                       seed = 1, intercept = TRUE) {
  heterogeneity <- match.arg(heterogeneity)
  set.seed(seed)

  level <- switch(
    heterogeneity,
    mild = list(mean_sd = 0.3, intercept_sd = 0.2, scale_min = 0.8, scale_max = 1.4),
    hard = list(mean_sd = 0.8, intercept_sd = 0.8, scale_min = 0.5, scale_max = 2.5),
    extreme = list(mean_sd = 1.4, intercept_sd = 1.5, scale_min = 0.3, scale_max = 4.0)
  )

  signs <- rep(c(1, -1), length.out = sparsity)
  beta_nonzero <- signs * seq(1.2, 0.4, length.out = sparsity)
  beta <- c(beta_nonzero, rep(0, max(0, p - sparsity)))
  beta <- beta[seq_len(p)]

  X_list <- vector("list", n_clients)
  y_list <- vector("list", n_clients)
  client_indices <- vector("list", n_clients)
  client_info <- data.frame()

  start <- 1
  for (j in seq_len(n_clients)) {
    n_j <- if (length(n_per_client) == 1) n_per_client else n_per_client[j]
    mean_shift <- rnorm(p, sd = level$mean_sd)
    mean_shift[seq(2, p, by = 2)] <- 0

    scale_j <- runif(1, level$scale_min, level$scale_max)
    intercept_j <- rnorm(1, sd = level$intercept_sd)
    tail_mix <- rbinom(n_j, size = 1, prob = min(0.7, 0.15 + 0.03 * j))

    X_j <- matrix(rnorm(n_j * p), nrow = n_j, ncol = p)
    X_j <- sweep(X_j, 2, mean_shift, "+")

    eps_j <- scale_j * (rexp(n_j) - log(2))
    eps_j <- eps_j + tail_mix * rexp(n_j, rate = 0.4)
    y_j <- as.numeric(X_j %*% beta + intercept_j + eps_j)

    X_list[[j]] <- X_j
    y_list[[j]] <- y_j
    idx <- start:(start + n_j - 1)
    client_indices[[j]] <- idx
    start <- start + n_j

    client_info <- rbind(client_info, data.frame(
      client = j,
      n = n_j,
      heterogeneity = heterogeneity,
      intercept_shift = intercept_j,
      noise_scale = scale_j,
      mean_shift_l2 = sqrt(sum(mean_shift^2)),
      y_mean = mean(y_j),
      y_q50 = as.numeric(quantile(y_j, 0.5)),
      y_q90 = as.numeric(quantile(y_j, 0.9))
    ))
  }

  X <- do.call(rbind, X_list)
  y <- unlist(y_list, use.names = FALSE)

  eps_q <- qexp(tau) - log(2)
  if (intercept) {
    X <- cbind("(Intercept)" = 1, X)
    beta <- c("(Intercept)" = eps_q, beta)
  }
  colnames(X) <- make.names(colnames(X), unique = TRUE)

  list(
    X = X,
    y = y,
    beta = beta,
    tau = tau,
    client_indices = client_indices,
    client_info = client_info,
    heterogeneity = heterogeneity
  )
}
