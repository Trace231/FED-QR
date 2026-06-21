# rfedqr 0.1.0

* Initial public release.
* Provides federated penalized quantile regression solvers with a unified
  `fit_fedqr()` interface.
* Adds QR box-dual primal-dual splitting, smoothing, subgradient, FedSPD-style,
  and ADMM baselines.
* Includes client-level loss diagnostics, quantile calibration utilities, and
  simulation and partition helpers.
* Includes tests for objectives, proximal operators, data partitions, solver
  interfaces, calibration utilities, and saved-result regressions when local
  experiment outputs are available.
