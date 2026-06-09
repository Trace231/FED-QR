test_that("check loss and QR objective are numerically correct", {
  u <- c(-2, -1, 0, 1, 3)
  expect_equal(check_loss(u, 0.75), c(0.5, 0.25, 0, 0.75, 2.25))

  X <- cbind(1, c(0, 1, 2))
  y <- c(1, 2, 4)
  beta <- c(1, 1)
  residual <- y - as.numeric(X %*% beta)
  expect_equal(qr_objective(X, y, beta, tau = 0.5), mean(abs(residual)) / 2)
})

test_that("prox operators return finite vectors with expected dimensions", {
  z <- c(-2, -0.2, 0, 0.4, 2)
  for (penalty in c("none", "l1", "mcp", "scad")) {
    out <- prox_penalty(
      z,
      gamma = 0.1,
      penalty = penalty,
      step_size = 0.5,
      lambda_value = 0.2
    )
    expect_length(out, length(z))
    expect_true(all(is.finite(out)))
  }
})

test_that("partitions cover samples without duplicates", {
  iid <- iid_partition(50, n_clients = 5, seed = 1)
  expect_equal(sort(unlist(iid, use.names = FALSE)), seq_len(50))

  y <- rnorm(60)
  noniid <- dirichlet_partition(y, n_clients = 6, alpha = 1, seed = 2)
  expect_equal(sort(unlist(noniid, use.names = FALSE)), seq_len(60))
})

test_that("QR box-dual keeps sample duals inside exact QR box", {
  dat <- make_qr_sim(n = 120, p = 5, tau = 0.9, seed = 10)
  clients <- iid_partition(nrow(dat$X), n_clients = 4, seed = 11)
  fit <- qr_box_fed_pdhg(
    dat$X,
    dat$y,
    client_indices = clients,
    tau = 0.9,
    rounds = 20,
    clients_per_round = 2,
    batch_size = 20,
    trace_every = 5,
    seed = 12
  )
  expect_true(all(fit$trace$dual_min >= 0.9 - 1 - 1e-12))
  expect_true(all(fit$trace$dual_max <= 0.9 + 1e-12))
})
