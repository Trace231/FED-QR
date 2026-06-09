# NYC Taxi Large-Scale Federated QR Experiment

## Purpose

This file records the first real large-scale distributed stochastic optimization experiment for the project.

HD is a real non-IID medical-center dataset, but it is small. NYC Taxi is used to test whether the proposed federated QR methods remain meaningful when the raw data scale is in the millions.

## Data Source

Official NYC TLC Yellow Taxi trip records:

- `yellow_tripdata_2024-01.parquet`
- `yellow_tripdata_2024-02.parquet`
- `yellow_tripdata_2024-03.parquet`
- `taxi_zone_lookup.csv`

Preparation script:

- `scripts/prepare_nyc_taxi_data.sh`

The script uses DuckDB CLI to query Parquet files without loading the raw files into R.

Generated files:

- `data/raw/nyc_taxi/yellow_tripdata_2024-01.parquet`
- `data/raw/nyc_taxi/yellow_tripdata_2024-02.parquet`
- `data/raw/nyc_taxi/yellow_tripdata_2024-03.parquet`
- `data/raw/nyc_taxi/taxi_zone_lookup.csv`
- `data/processed/nyc_taxi_yellow_qr_2024q1_sample.csv`
- `data/processed/nyc_taxi_yellow_qr_2024q1_design_summary.csv`

Raw valid trips after filtering:

| raw clean rows | clean pickup-zone clients | eligible clients | modeled rows |
|---:|---:|---:|---:|
| 8,417,330 | 257 | 59 | 1,000,000 |

The modeling sample has:

| n | p | clients | min client n | median client n | max client n |
|---:|---:|---:|---:|---:|---:|
| 1,000,000 | 20 | 59 | 639 | 14,755 | 50,432 |

## Modeling Task

Response:

- `log_total_amount`

Quantile:

- `tau = 0.9`

Features:

- standardized log trip distance;
- standardized log duration;
- standardized passenger count;
- standardized congestion surcharge and airport fee;
- pickup hour sine/cosine;
- pickup day-of-week sine/cosine;
- pickup borough;
- vendor ID;
- rate code;
- payment type.

Federated clients:

- pickup zones, using `PULocationID`.

This is a natural non-IID split: pickup zones differ strongly in airport traffic, Manhattan density, trip length, and fare distribution.

## Experiment

Main script:

- `scripts/run_nyc_taxi_large_scale_experiment.R`

Final saved setting:

- 1,000,000 modeled rows;
- 59 clients;
- 12 clients per communication round;
- local user batch size = 5,000;
- base rounds = 500;
- QR box-dual long = 1,000 rounds;
- `tau = 0.9`.

Outputs:

- `results/nyc_taxi_large_scale_design_batch5000.csv`
- `results/nyc_taxi_large_scale_clients.csv`
- `results/nyc_taxi_large_scale_summary_batch5000.csv`
- `results/nyc_taxi_large_scale_trace_batch5000.csv`
- `results/nyc_taxi_large_scale_coefficients.csv`
- `figures/nyc_taxi_large_scale_convergence_batch5000.png`
- `figures/nyc_taxi_large_scale_final_gap_batch5000.png`

## Results

| method | objective | gap to best observed |
|---|---:|---:|
| QR box-dual long | 0.0206365 | 0 |
| QR box-dual | 0.0220496 | 0.0014131 |
| FSPG-smooth | 0.0240271 | 0.0033905 |
| FedSubGrad | 0.1043675 | 0.0837310 |
| FedSPD-check | 0.3452998 | 0.3246633 |

## Interpretation

The first small-batch NYC run showed that default QR box-dual is not automatically stable at large scale. With local user batch size 200, FSPG-smooth was much better.

After increasing the local user batch size to 5,000, the result changed:

- QR box-dual long achieved the best observed full-sample objective.
- FSPG-smooth remained a strong short-horizon baseline, but was no longer best.
- FedSubGrad and FedSPD-check were much worse at this scale.

This gives a sharper and more realistic conclusion:

> Large-scale federated QR requires enough local user-level dual updates per communication round. In the million-row NYC Taxi experiment, cached QR box-dual becomes the best method after increasing the local mini-batch size, supporting its role as a real distributed stochastic optimization method rather than only a small-simulation algorithm.

Important caveat:

- The NYC experiment uses best observed objective among distributed methods, not a centralized exact optimum. A centralized exact QR solve at this scale is not the target and would be computationally expensive.
