# Experiment Manifest

This manifest maps the main result files to the scripts and settings that generate them. Large data files and generated outputs are retained locally but excluded from R package builds.

## Core Method Validation

| Purpose | Script | Main outputs |
|---|---|---|
| Central and federated QR box-dual validation | `scripts/run_box_dual_comparison.R` | `results/box_dual_comparison_summary.csv`, `results/box_dual_comparison_trace.csv` |
| FedSPD-DP reproduction baseline | `scripts/run_fedspd_dp_reproduction.R` | `results/fedspd_dp_reproduction_summary.csv`, `results/fedspd_dp_reproduction_trace.csv` |
| ADMM validation | `scripts/run_admm_validation.R` | `results/admm_validation_summary.csv` |

## Simulation Stress Tests

| Purpose | Script | Main outputs |
|---|---|---|
| R1/R2 isolation | `scripts/run_r1_r2_stress_test.R` | `results/r1_r2_stress_summary.csv`, `results/r1_r2_stress_trace.csv` |
| Hard non-IID test | `scripts/run_hard_noniid_with_plots.R` | `results/hard_noniid_summary.csv`, `figures/hard_noniid_final_gap.png` |
| Baseline comparison under hard non-IID | `scripts/run_hard_baseline_comparison_with_plots.R` | `results/hard_baseline_summary.csv`, `figures/hard_baseline_final_gap.png` |
| Expanded scenario matrix | `scripts/run_expanded_scenario_suite_with_plots.R` | `results/expanded_scenario_summary.csv`, `figures/expanded_scenario_rank_heatmap.png` |
| QR tuning across scenarios | `scripts/run_qr_box_tuning_expanded_scenarios.R` | `results/qr_tuning_expanded_summary.csv`, `figures/qr_tuning_expanded_gap.png` |
| Advanced innovation variants | `scripts/run_advanced_innovation_experiment_with_plots.R` | `results/advanced_innovation_summary.csv`, `results/advanced_innovation_client_loss.csv`, `results/advanced_innovation_calibration.csv`, `figures/advanced_innovation_final_gap.png` |

## Penalty Experiments

| Purpose | Script | Main outputs |
|---|---|---|
| L1 penalized QR | `scripts/run_l1_penalized_experiment_with_plots.R` | `results/l1_penalized_summary.csv`, `figures/l1_penalized_selection.png` |
| MCP/SCAD nonconvex penalties | `scripts/run_nonconvex_penalty_experiment_with_plots.R` | `results/nonconvex_penalty_summary.csv`, `figures/nonconvex_penalty_selection.png` |

## Real Data

| Dataset | Script | Setting | Main outputs |
|---|---|---|---|
| UCI Heart Disease | `scripts/run_hd_baseline_comparison_with_plots.R` | 4 natural centers, `tau = 0.5, 0.75, 0.9` | `results/hd_baseline_aggregate.csv`, `figures/hd_baseline_r1r2_final_gap.png` |
| NYC Taxi | `scripts/run_nyc_taxi_large_scale_experiment.R` | 1,000,000 modeled trips, pickup-zone clients, `tau = 0.9` | `results/nyc_taxi_large_scale_summary_batch5000.csv`, `figures/nyc_taxi_large_scale_final_gap_batch5000.png` |
| UCI Household Power | `scripts/run_household_power_large_scale_experiment.R` | 2,049,280 modeled rows, calendar-month clients, `tau = 0.9` | `results/household_power_large_scale_summary.csv`, `figures/household_power_large_scale_final_gap.png` |
| Large-data coefficient diagnostics | `scripts/run_large_realdata_coefficient_evaluation.R` | Reuses saved NYC and Household coefficients, no retraining | `results/advanced_innovation_large_client_loss.csv`, `results/advanced_innovation_large_calibration.csv` |

## Current Large-Data Conclusions

- NYC Taxi batch-size 5,000 run: `QR box-dual long` is best observed.
- Household Power full-data run: `QR box-dual long` is best observed, with `FSPG-smooth` very close.
- Heart Disease R1+R2 aggregate: `QR box-dual` beats `FedSPD-check` across the tested quantiles.
- Staleness-aware cached directions improve the QR box-dual optimization gap in the extreme non-IID high-quantile simulation.
- Client-robust weighting improves worst-client loss on Heart Disease and several unbalanced simulation settings, even when the global mean objective is not the winner.
- Intercept and client-offset calibration sharply reduce high-quantile coverage error on simulation, Heart Disease, NYC Taxi, and Household Power coefficient diagnostics.
- The `QR box-dual adaptive` variant is registered for future advanced-experiment reruns; it adapts stale decay and client weights from the current cache-age and per-client loss state.

## Fast Verification

Use the test suite for quick checks:

```bash
Rscript --vanilla -e 'if (requireNamespace("testthat", quietly=TRUE)) testthat::test_dir("tests/testthat")'
```

Use full experiment scripts only when refreshing tables or final figures.
