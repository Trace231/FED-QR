test_that("fit_fedqr returns normalized output for registered methods", {
  dat <- make_qr_sim(n = 160, p = 5, tau = 0.75, seed = 20)
  clients <- iid_partition(nrow(dat$X), n_clients = 4, seed = 21)
  initial <- qr_objective(dat$X, dat$y, numeric(ncol(dat$X)), tau = 0.75)

  controls <- list(
    "QR box-dual" = list(),
    "FSPG-smooth" = list(smooth_mu = 0.05),
    "FedSubGrad" = list(),
    "FedSPD-check" = list(Q = 2, gamma0 = 2),
    "FedQR-ADMM" = list(rho_consensus = 1, rho_residual = 0.3, inner_iter = 5)
  )

  for (method in names(controls)) {
    fit <- fit_fedqr(
      method,
      dat$X,
      dat$y,
      client_indices = clients,
      tau = 0.75,
      rounds = 8,
      clients_per_round = 2,
      batch_size = 20,
      trace_every = 4,
      seed = 22,
      control = controls[[method]]
    )
    expect_type(fit$beta, "double")
    expect_length(fit$beta, ncol(dat$X))
    expect_true(is.finite(fit$objective))
    expect_true(fit$objective < initial + 10)
    expect_true(all(c("round", "objective", "method") %in% names(fit$trace)))
  }
})

test_that("run_fedqr_methods builds summary, trace, and coefficients", {
  dat <- make_qr_sim(n = 140, p = 4, tau = 0.5, seed = 30)
  clients <- iid_partition(nrow(dat$X), n_clients = 4, seed = 31)
  result <- run_fedqr_methods(
    c("QR box-dual", "FSPG-smooth", "FedSubGrad"),
    dat$X,
    dat$y,
    client_indices = clients,
    tau = 0.5,
    rounds = 10,
    clients_per_round = 2,
    batch_size = 20,
    trace_every = 5,
    seed = 32,
    verbose = FALSE
  )
  expect_true(all(c("summary", "trace", "coefficients", "fits") %in% names(result)))
  expect_equal(nrow(result$summary), 3)
  expect_equal(nrow(result$coefficients), ncol(dat$X))
  expect_true(all(result$summary$gap_to_best_observed >= -1e-12))
  expect_silent(validate_result_table(result$summary))
})
