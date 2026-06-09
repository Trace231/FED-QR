test_that("default QR box-dual behavior is unchanged by new options", {
  dat <- make_qr_sim(n = 120, p = 4, tau = 0.9, seed = 101)
  clients <- iid_partition(nrow(dat$X), n_clients = 4, seed = 102)
  old_default <- qr_box_fed_pdhg(
    dat$X, dat$y,
    client_indices = clients,
    tau = 0.9,
    rounds = 12,
    clients_per_round = 2,
    batch_size = 20,
    trace_every = 4,
    seed = 103
  )
  explicit_default <- qr_box_fed_pdhg(
    dat$X, dat$y,
    client_indices = clients,
    tau = 0.9,
    rounds = 12,
    clients_per_round = 2,
    batch_size = 20,
    trace_every = 4,
    seed = 103,
    staleness = "none",
    client_weighting = "sample"
  )
  expect_equal(old_default$beta, explicit_default$beta)
  expect_equal(old_default$objective, explicit_default$objective)
})

test_that("staleness weights and client weights are well behaved", {
  dat <- make_qr_sim(n = 150, p = 5, tau = 0.9, seed = 111)
  clients <- iid_partition(nrow(dat$X), n_clients = 5, seed = 112)
  fit <- qr_box_fed_pdhg(
    dat$X, dat$y,
    client_indices = clients,
    tau = 0.9,
    rounds = 15,
    clients_per_round = 2,
    batch_size = 20,
    trace_every = 5,
    seed = 113,
    staleness = "exponential",
    staleness_floor = 0.4,
    client_weighting = "uniform"
  )
  expect_true(all(fit$trace$mean_stale_weight >= 0.4 - 1e-12))
  expect_true(all(fit$trace$mean_stale_weight <= 1 + 1e-12))
  expect_equal(sum(fit$client_weights), 1)
  expect_equal(as.numeric(fit$client_weights), rep(1 / 5, 5))
})

test_that("client objective and calibration summaries are finite", {
  dat <- make_qr_sim(n = 180, p = 5, tau = 0.9, seed = 121)
  clients <- iid_partition(nrow(dat$X), n_clients = 6, seed = 122)
  beta <- dat$beta + rnorm(length(dat$beta), sd = 0.25)

  obj_sample <- client_qr_objective(
    dat$X, dat$y, beta,
    client_indices = clients,
    tau = 0.9,
    client_weighting = "sample"
  )
  obj_global <- qr_objective(dat$X, dat$y, beta, tau = 0.9)
  expect_equal(obj_sample, obj_global)

  loss <- client_loss_summary(dat$X, dat$y, beta, clients, tau = 0.9)
  expect_true(all(is.finite(unlist(loss))))

  before <- calibration_summary(dat$X, dat$y, beta, tau = 0.9, client_indices = clients)
  cal <- calibrate_quantile_intercept(dat$X, dat$y, beta, tau = 0.9, mode = "global")
  after <- calibration_summary(dat$X, dat$y, cal$beta, tau = 0.9, client_indices = clients)
  expect_lte(after$global_coverage_error, before$global_coverage_error + 1e-12)

  client_cal <- calibrate_quantile_intercept(
    dat$X, dat$y, beta,
    tau = 0.9,
    client_indices = clients,
    mode = "client_offset"
  )
  client_after <- calibration_summary(
    dat$X, dat$y, client_cal$beta,
    tau = 0.9,
    client_indices = clients,
    offsets = client_cal$offsets
  )
  expect_lte(client_after$mean_client_coverage_error,
             before$mean_client_coverage_error + 1e-12)
})

test_that("advanced QR box-dual method aliases run through fit_fedqr", {
  dat <- make_qr_sim(n = 120, p = 4, tau = 0.9, seed = 131)
  clients <- iid_partition(nrow(dat$X), n_clients = 4, seed = 132)
  methods <- c("QR box-dual stale", "QR box-dual robust", "QR box-dual stale+robust")
  for (method in methods) {
    fit <- fit_fedqr(
      method,
      dat$X,
      dat$y,
      client_indices = clients,
      tau = 0.9,
      rounds = 8,
      clients_per_round = 2,
      batch_size = 20,
      trace_every = 4,
      seed = 133
    )
    expect_true(is.finite(fit$objective))
    expect_true(all(c("mean_staleness", "mean_stale_weight", "client_objective") %in% names(fit$trace)))
  }
})
