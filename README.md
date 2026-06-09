# FED-QR

`rfedqr` is a research implementation of federated penalized quantile regression, centered on a QR-specific box-dual primal-dual method and a set of non-toy baselines.

The project is organized as both:

- an installable R package exposing reusable algorithm interfaces;
- a reproducible experiment bundle for simulations, Heart Disease, NYC Taxi, and UCI Household Power.

## Main Algorithm

The main method is `QR box-dual`, implemented by `qr_box_fed_pdhg()`.

It uses the exact quantile check-loss dual representation:

```text
rho_tau(u) = max_{v in [tau - 1, tau]} v u
```

and maintains sample-level dual variables on clients with cached inactive-client directions.

## Minimal Example

```r
library(rfedqr)

dat <- make_qr_sim(n = 300, p = 8, tau = 0.9, seed = 1)
clients <- iid_partition(nrow(dat$X), n_clients = 6, seed = 2)

fit <- fit_fedqr(
  "QR box-dual",
  dat$X,
  dat$y,
  client_indices = clients,
  tau = 0.9,
  rounds = 100,
  clients_per_round = 3,
  batch_size = 30,
  seed = 3
)

fit$objective
```

To compare methods:

```r
result <- run_fedqr_methods(
  c("QR box-dual", "QR box-dual stale", "QR box-dual robust", "FSPG-smooth"),
  dat$X,
  dat$y,
  client_indices = clients,
  tau = 0.9,
  rounds = 100,
  clients_per_round = 3,
  batch_size = 30
)

result$summary
```

## Public Interfaces

- `fit_fedqr()`: normalized interface for one method.
- `run_fedqr_methods()`: runs multiple methods and returns `summary`, `trace`, and `coefficients`.
- `qr_box_fed_pdhg()`: main QR box-dual federated PDHG solver.
- `QR box-dual stale`, `QR box-dual robust`, `QR box-dual stale+robust`: optional advanced variants with staleness-aware cached directions and client-robust weighting.
- `client_loss_summary()` and `client_qr_objective()`: client-level fairness and robust-objective diagnostics.
- `quantile_coverage()`, `calibrate_quantile_intercept()`, `calibration_summary()`: high-quantile coverage and calibration refinement tools.
- `fedspd_dp_qr()`: FedSPD-DP-style consensus baseline.
- `fed_subgrad_qr()` and `fed_smooth_qr()`: subgradient and smoothing baselines.
- `fed_qr_admm()`: federated ADMM baseline.
- `make_experiment_dirs()`, `save_experiment_outputs()`, `plot_convergence()`, `plot_final_gap()`: experiment utilities.

## Reproducing Experiments

Fast simulation and validation:

```bash
Rscript --vanilla scripts/run_box_dual_comparison.R
Rscript --vanilla scripts/run_hd_baseline_comparison_with_plots.R
Rscript --vanilla scripts/run_advanced_innovation_experiment_with_plots.R
Rscript --vanilla scripts/run_large_realdata_coefficient_evaluation.R
```

Large real-data experiments:

```bash
scripts/prepare_nyc_taxi_data.sh
NYC_ROUNDS=500 NYC_BATCH_SIZE=5000 Rscript --vanilla scripts/run_nyc_taxi_large_scale_experiment.R

scripts/prepare_household_power_data.sh
POWER_ROUNDS=400 POWER_BATCH_SIZE=5000 Rscript --vanilla scripts/run_household_power_large_scale_experiment.R
```

The large processed data and generated results are kept locally in `data/`, `results/`, and `figures/`. They are intentionally excluded from R package builds via `.Rbuildignore`.

## Testing

```bash
Rscript --vanilla -e 'if (requireNamespace("testthat", quietly=TRUE)) testthat::test_dir("tests/testthat")'
R CMD check .
```

The automated tests are layered:

- unit tests for objectives, proximal operators, partitions, and dual-box validity;
- small smoke tests for the unified method interface;
- regression checks against saved large-experiment summaries when those files are present locally.
