# UCI Household Power Large-Scale Real-Data Experiment

## Purpose

This experiment adds a second non-NYC large-scale real-data benchmark for the federated quantile regression project. NYC Taxi tests spatially distributed transportation data; this benchmark tests a household electricity time series with strong temporal drift.

Data source:

- UCI Individual Household Electric Power Consumption dataset;
- minute-level measurements from one household;
- about 2.08 million raw records before removing missing values.

## Files

Preparation:

- `scripts/prepare_household_power_data.sh`

Experiment:

- `scripts/run_household_power_large_scale_experiment.R`

Outputs:

- `data/processed/household_power_consumption.csv`
- `results/household_power_large_scale_design.csv`
- `results/household_power_large_scale_clients.csv`
- `results/household_power_large_scale_summary.csv`
- `results/household_power_large_scale_trace.csv`
- `results/household_power_large_scale_coefficients.csv`
- `figures/household_power_client_heterogeneity.png`
- `figures/household_power_large_scale_convergence.png`
- `figures/household_power_large_scale_final_gap.png`

## Modeling Setup

Response:

- `log1p(Global_active_power)`

Quantile:

- `tau = 0.9`

Clients:

- calendar months;
- 48 monthly clients after cleaning.

Features:

- standardized global reactive power;
- standardized voltage;
- standardized sub-metering 1, 2, and 3;
- hour-of-day sine/cosine;
- day-of-week sine/cosine;
- month factor.

We intentionally exclude `Global_intensity` because it is too directly tied to active power and would make the prediction task less informative.

## Distributed Stochastic Setting

Final full-data run:

| quantity | value |
|---|---:|
| modeled rows | 2,049,280 |
| features | 21 |
| clients | 48 |
| min client size | 21,992 |
| median client size | 43,474.5 |
| max client size | 44,640 |
| clients per round | 10 |
| local user batch size | 5,000 |
| base rounds | 400 |
| QR long rounds | 800 |

This is a real distributed stochastic optimization setting:

- R1: only 10 of 48 monthly clients participate per communication round;
- R2: each active client updates only 5,000 local observations per round;
- clients are non-IID because electricity usage changes across seasons and months.

## Results

| method | objective | gap to best observed |
|---|---:|---:|
| QR box-dual long | 0.0419800 | 0 |
| FSPG-smooth | 0.0423711 | 0.0003911 |
| FedSubGrad | 0.0456561 | 0.0036761 |
| QR box-dual | 0.0500239 | 0.0080439 |
| FedSPD-check | 0.0550988 | 0.0131188 |

## Interpretation

The result is deliberately more nuanced than a toy win.

- QR box-dual long gives the best observed objective on the full 2.05-million-row task.
- FSPG-smooth is very close under the shorter 400-round budget, confirming again that smoothing is a strong short-horizon baseline.
- The 400-round QR box-dual run underfits relative to the long run, which is consistent with the sample-level dual-refresh interpretation: the method needs enough local dual coverage when there are millions of observations.
- FedSPD-check remains clearly behind because it is a direct nonsmooth consensus adaptation and does not use the exact QR box-dual geometry.

Current report wording:

> On a second large-scale real dataset, UCI Household Power, cached QR box-dual again achieves the best observed objective after a longer communication horizon. The gap between QR box-dual long and FSPG-smooth is small, which is useful evidence rather than a weakness: smoothing is a genuinely competitive short-budget baseline, while QR box-dual provides the best final nonsmoothed check-loss objective after enough dual refresh.

## Reproduction

```bash
scripts/prepare_household_power_data.sh

POWER_ROUNDS=400 \
POWER_FEDSPD_ROUNDS=400 \
POWER_TRACE_EVERY=25 \
POWER_CLIENTS_PER_ROUND=10 \
POWER_BATCH_SIZE=5000 \
POWER_TAU=0.9 \
Rscript --vanilla scripts/run_household_power_large_scale_experiment.R
```
