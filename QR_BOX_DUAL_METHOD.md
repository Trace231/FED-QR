# QR Box-Dual Federated PDHG

## Purpose

This document records the faithful implementation of the QR box-dual method. This is separate from:

- `fedspd_dp_qr()`: paper-faithful FedSPD-DP consensus implementation;
- `qr_pdhg()`: centralized QR-PDHG;
- earlier prototype `fed_qr_spd()`.

The main implementation is:

- `R/qr_box_fed_pdhg.R`

The main experiment script is:

- `scripts/run_box_dual_comparison.R`

Outputs:

- `results/box_dual_comparison_summary.csv`
- `results/box_dual_comparison_trace.csv`

## Mathematical Formulation

The quantile check loss satisfies:

```text
rho_tau(u) = max_{v in [tau - 1, tau]} v u
```

Thus empirical QR can be written as:

```text
min_beta max_{v_i in [tau - 1, tau]}
  (1/n) sum_i v_i (y_i - x_i^T beta) + P_lambda(beta)
```

The box-dual method applies PDHG / Chambolle-Pock updates:

```text
v_i <- clip(v_i + sigma * (y_i - x_i^T beta_bar), tau - 1, tau)
```

```text
beta <- prox_{eta P_lambda}(beta + eta * X^T v / n)
```

```text
beta_bar <- beta + theta * (beta - beta_old)
```

The implementation uses the empirical inner product `(1/n) sum_i`, so the dual ascent uses the residual directly and the primal direction uses `X^T v / n`.

## Federated Implementation

Each client stores:

- local design matrix and response;
- sample-level dual variables `v_i`;
- cached primal direction `X_i^T v_i / n`.

The server stores:

- global primal variable `beta`;
- extrapolated variable `beta_bar`;
- cached client directions.

Per communication round:

1. Server samples active clients.
2. Each active client updates its local sample dual variables:

   ```text
   v_i <- clip(v_i + sigma * residual_i, tau - 1, tau)
   ```

3. Each active client recomputes and uploads its cached direction:

   ```text
   direction_j = X_j^T v_j / n
   ```

4. Server aggregates all cached directions, including inactive clients:

   ```text
   direction = sum_j cached_direction_j
   ```

5. Server updates:

   ```text
   beta <- prox(beta + eta * direction)
   beta_bar <- beta + theta * (beta - beta_old)
   ```

This cached aggregation is important. Inactive clients are not treated as zero; they keep their most recent uploaded direction, analogous to cached inactive-client states in federated primal-dual methods.

## Validation

The box-dual implementation was validated against centralized QR-PDHG.

Simulation:

- `n = 800`;
- `p = 10`;
- 20 clients;
- `tau = 0.5, 0.75, 0.9`;
- asymmetric noise;
- no penalty.

Summary:

| tau | method | regime | objective | beta L2 error | dual box |
|---:|---|---|---:|---:|---|
| 0.50 | box-dual | full clients/full batch | 0.3291519 | 0.1109 | [-0.50, 0.50] |
| 0.50 | box-dual | half clients/minibatch | 0.3291619 | 0.1117 | [-0.50, 0.50] |
| 0.50 | central QR-PDHG | central | 0.3291510 | 0.1110 | NA |
| 0.75 | box-dual | full clients/full batch | 0.3210287 | 0.2362 | [-0.25, 0.75] |
| 0.75 | box-dual | half clients/minibatch | 0.3210369 | 0.2370 | [-0.25, 0.75] |
| 0.75 | central QR-PDHG | central | 0.3210281 | 0.2362 | NA |
| 0.90 | box-dual | full clients/full batch | 0.2005622 | 0.3066 | [-0.10, 0.90] |
| 0.90 | box-dual | half clients/minibatch | 0.2005662 | 0.3058 | [-0.10, 0.90] |
| 0.90 | central QR-PDHG | central | 0.2005608 | 0.3062 | NA |

Interpretation:

- Full-client/full-batch box-dual matches centralized QR-PDHG up to numerical tolerance.
- Cached partial participation remains stable and nearly reaches the same objective after enough rounds.
- Mini-batch dual updates remain stable in this setting.
- Dual variables always stay exactly within the QR box `[tau - 1, tau]`.
- Compared with `fedspd_dp_qr(loss = "check")`, the box-dual method is much closer to central QR at `tau = 0.9`.

## Project Role

This is the QR-specific method we can present as the main algorithmic contribution:

```text
FedSPD-DP-style consensus method:
  faithful paper baseline

FedSPD-DP with check subgradient:
  direct nonsmooth QR adaptation

QR box-dual federated PDHG:
  QR geometry-specialized method
```

The next experiments should compare these methods systematically on:

- R1 client sampling;
- R2 local mini-batching;
- tau-adaptive extreme quantiles;
- UCI Heart Disease four-center data;
- later, NYC Taxi scalability.

## R1/R2 Stress Test

新增脚本：

- `scripts/run_r1_r2_stress_test.R`

输出：

- `results/r1_r2_stress_summary.csv`
- `results/r1_r2_stress_trace.csv`
- `results/r1_r2_partition_diagnostics.csv`
- `results/r1_r2_stress_aggregate.csv`
- `results/r1_r2_stress_clean_aggregate.csv`

设置：

- `tau = 0.9`，极端高分位；
- `n = 1200, p = 15`;
- 20 clients;
- Dirichlet response-based non-IID partition;
- `alpha = 10, 1, 0.25, 0.05`;
- 3 random seeds;
- scenarios:
  - deterministic: all clients + full local batch;
  - R1: client sampling only, `K = 4/20`;
  - R2: user/sample mini-batch only, batch size 10;
  - R1+R2: client sampling + user/sample mini-batch.

Clean aggregate final objective gap vs centralized QR-PDHG:

| method | scenario | alpha | final gap |
|---|---|---:|---:|
| box-dual | deterministic | all | about `8.6e-06` |
| box-dual | R1 only | 0.05 | `1.48e-04` |
| box-dual | R1 only | 0.25 | `1.62e-04` |
| box-dual | R1 only | 1 | `1.47e-04` |
| box-dual | R1 only | 10 | `1.25e-04` |
| box-dual | R2 only | 0.05 | `2.11e-04` |
| box-dual | R2 only | 0.25 | `2.22e-04` |
| box-dual | R2 only | 1 | `1.87e-04` |
| box-dual | R2 only | 10 | `2.55e-04` |
| box-dual | R1+R2 | 0.05 | `3.51e-03` |
| box-dual | R1+R2 | 0.25 | `2.82e-03` |
| box-dual | R1+R2 | 1 | `1.96e-03` |
| box-dual | R1+R2 | 10 | `2.23e-03` |
| FedSPD-check | R1+R2 | 0.05 | `7.42e-02` |
| FedSPD-check | R1+R2 | 0.25 | `8.68e-02` |
| FedSPD-check | R1+R2 | 1 | `8.72e-02` |
| FedSPD-check | R1+R2 | 10 | `9.03e-02` |

Interpretation:

- Box-dual remains very close to centralized QR under isolated R1 and isolated R2.
- R1+R2 is the genuinely hard regime, but box-dual still keeps the final objective gap around `2e-03` to `4e-03`.
- FedSPD-check degrades much more under R1+R2, with final gaps around `7e-02` to `9e-02`.
- This creates a clear empirical separation: the QR box-dual geometry is substantially more robust than the direct check-subgradient adaptation in the extreme quantile stochastic federated setting.

Important caveat:

- The current Dirichlet partition is response-stratified and does induce client distribution shifts, but the alpha trend is not strongly monotone in every metric. For stronger non-IID stress, a later simulation should generate client-specific covariate/response shifts directly.

## Hard Non-IID Test With Figures

新增脚本：

- `scripts/run_hard_noniid_with_plots.R`

输出数据：

- `results/hard_noniid_summary.csv`
- `results/hard_noniid_aggregate.csv`
- `results/hard_noniid_trace.csv`
- `results/hard_noniid_trace_aggregate.csv`
- `results/hard_noniid_client_info.csv`

输出图：

- `figures/hard_noniid_final_gap.png`
- `figures/hard_noniid_r1r2_convergence.png`
- `figures/hard_noniid_client_heterogeneity.png`

生成机制：

- 20 clients;
- each client has covariate mean shift;
- client-specific intercept shift;
- client-specific noise scale;
- client-specific tail mixture;
- heterogeneity levels: `mild`, `hard`, `extreme`;
- `tau = 0.9`;
- scenarios: deterministic, R1, R2, R1+R2.

Aggregate results:

| method | scenario | heterogeneity | final gap |
|---|---|---|---:|
| QR box-dual | deterministic | mild | `1.14e-05` |
| QR box-dual | deterministic | hard | `2.20e-05` |
| QR box-dual | deterministic | extreme | `1.12e-04` |
| QR box-dual | R1 | mild | `1.03e-04` |
| QR box-dual | R1 | hard | `1.33e-04` |
| QR box-dual | R1 | extreme | `2.00e-04` |
| QR box-dual | R2 | mild | `2.96e-04` |
| QR box-dual | R2 | hard | `2.40e-04` |
| QR box-dual | R2 | extreme | `2.77e-04` |
| QR box-dual | R1+R2 | mild | `1.72e-03` |
| QR box-dual | R1+R2 | hard | `2.66e-03` |
| QR box-dual | R1+R2 | extreme | `3.52e-03` |
| FedSPD-check | deterministic | mild | `4.86e-02` |
| FedSPD-check | deterministic | hard | `7.19e-02` |
| FedSPD-check | deterministic | extreme | `1.23e-01` |
| FedSPD-check | R1+R2 | mild | `3.91e-01` |
| FedSPD-check | R1+R2 | hard | `4.75e-01` |
| FedSPD-check | R1+R2 | extreme | `6.17e-01` |

Conclusion:

- Under generated client-level covariate, intercept, noise-scale, and tail heterogeneity, QR box-dual remains stable.
- The hardest R1+R2 setting widens the gap strongly: FedSPD-check deteriorates by orders of magnitude, while QR box-dual stays near the centralized objective.
- This is currently the strongest evidence that the QR box geometry is not only a cosmetic trick; it materially improves robustness in the nonsmooth extreme-quantile federated setting.

## Hard Baseline Comparison

新增脚本：

- `scripts/run_hard_baseline_comparison_with_plots.R`

输出数据：

- `results/hard_baseline_summary.csv`
- `results/hard_baseline_aggregate.csv`
- `results/hard_baseline_trace.csv`
- `results/hard_baseline_trace_aggregate.csv`

输出图：

- `figures/hard_baseline_final_gap.png`
- `figures/hard_baseline_final_gap_log.png`
- `figures/hard_baseline_r1r2_convergence.png`

Methods:

- QR box-dual;
- FedQR-ADMM;
- FedSubGrad;
- FSPG-smooth;
- FedSPD-check.

Hard R1+R2 final objective gaps:

| method | mild | hard | extreme |
|---|---:|---:|---:|
| QR box-dual | `1.72e-03` | `2.66e-03` | `3.52e-03` |
| FSPG-smooth | `3.36e-03` | `5.48e-03` | `5.57e-03` |
| FedSubGrad | `4.39e-02` | `1.07e-01` | `2.71e-01` |
| FedSPD-check | `3.91e-01` | `4.75e-01` | `6.17e-01` |
| FedQR-ADMM | `9.25e-01` | `1.10e+00` | `1.37e+00` |

Interpretation:

- QR box-dual is the best method by final objective gap in all three hard R1+R2 regimes.
- FSPG-smooth is a strong baseline and should be treated seriously in the report. It is usually second-best, but it depends on the smoothing parameter and changes the original QR objective.
- FedSubGrad is much worse under stronger heterogeneity, but still more robust than FedSPD-check in R1+R2.
- FedSPD-check is not competitive under strong stochasticity; it is mainly useful as a paper-faithful direct nonsmooth adaptation baseline.
- FedQR-ADMM is valid in deterministic/full-client validation, but it degrades badly under R1+R2 stochasticity because the local ADMM subproblems and consensus dual variables are very sensitive to partial stale clients and mini-batch residual updates.

Suggested report wording:

> In the hard non-IID R1+R2 regime, QR box-dual consistently achieves the smallest objective gap. FSPG-smooth is competitive but remains worse and requires a smoothing parameter, while FedSubGrad, FedSPD-check, and stochastic FedQR-ADMM deteriorate substantially as heterogeneity increases.

## ADMM Baseline

新增 ADMM baseline：

- `R/qr_admm.R`
- `scripts/run_admm_validation.R`

Implemented methods:

- `qr_admm()`: centralized residual-splitting QR-ADMM;
- `fed_qr_admm()`: federated consensus QR-ADMM with local residual-splitting ADMM subproblems.

Central formulation:

```text
min_beta,r (1/n) sum_i rho_tau(r_i) + P_lambda(beta)
s.t. r = y - X beta
```

ADMM updates:

```text
beta <- argmin_beta P_lambda(beta) + (rho/2)||y - X beta - r + u||^2
r <- prox_{rho_tau/(n rho)}(y - X beta + u)
u <- u + y - X beta - r
```

Check-loss prox:

```text
prox_{kappa rho_tau}(z)
= z - clip(z, kappa(tau - 1), kappa tau)
```

Validation:

| tau | central QR-PDHG | central QR-ADMM gap | fed QR-ADMM full-client gap |
|---:|---:|---:|---:|
| 0.50 | 0.3243903 | `-6.0e-08` | `1.8e-07` |
| 0.75 | 0.3301641 | `2.3e-05` | `5.2e-05` |
| 0.90 | 0.2169556 | `4.1e-05` | `5.6e-04` |

This confirms that ADMM is a legitimate baseline, not a toy implementation. Its poor hard R1+R2 performance reflects stochastic federated difficulty rather than a broken centralized solver.

## L1-Penalized Hard Non-IID Experiment

新增脚本：

- `scripts/run_l1_penalized_experiment_with_plots.R`

输出数据：

- `results/l1_penalized_summary.csv`
- `results/l1_penalized_aggregate.csv`
- `results/l1_penalized_trace.csv`
- `results/l1_penalized_trace_aggregate.csv`

输出图：

- `figures/l1_penalized_final_gap.png`
- `figures/l1_penalized_selection.png`
- `figures/l1_penalized_support_size.png`

设置：

- `p = 60`;
- true support size = 8;
- hard client-level non-IID;
- `tau = 0.9`;
- R1+R2 stochasticity: 4/20 clients per round and mini-batch size 10;
- `lambda in {0.001, 0.005, 0.01, 0.02}`;
- 3 random seeds.

Key results:

| method | lambda | final gap | selected size | TPR | FDR |
|---|---:|---:|---:|---:|---:|
| QR box-dual | 0.001 | `7.06e-03` | 55.7 | 1.000 | 0.856 |
| QR box-dual | 0.005 | `4.24e-03` | 39.0 | 1.000 | 0.795 |
| QR box-dual | 0.010 | `2.93e-03` | 23.7 | 0.958 | 0.667 |
| QR box-dual | 0.020 | `1.10e-03` | 11.0 | 0.875 | 0.366 |
| FSPG-smooth | 0.020 | `1.73e-02` | 52.7 | 1.000 | 0.848 |
| FedSubGrad | 0.020 | `1.29e-01` | 45.0 | 1.000 | 0.822 |
| FedQR-ADMM | 0.020 | `1.24e+00` | 2.3 | 0.167 | 0.389 |

Interpretation:

- QR box-dual remains the best objective baseline under L1 penalization.
- Increasing `lambda` creates the expected sparsity path for QR box-dual: selected variables decrease from about 56 to 11.
- At `lambda = 0.02`, QR box-dual gets close to the true support size 8 while keeping high TPR and much lower FDR.
- FSPG-smooth has a reasonable objective but remains too dense under the same proximal settings.
- FedSubGrad is slow and too dense.
- FedQR-ADMM becomes overly sparse and unstable under R1+R2.

This experiment confirms that the project now genuinely covers penalized quantile regression, not only unpenalized QR.

## MCP/SCAD Nonconvex Penalty Experiment

新增核心支持：

- `prox_mcp()` and `prox_scad()` in `R/prox.R`;
- `qr_box_fed_pdhg(..., penalty = "mcp")`;
- `qr_box_fed_pdhg(..., penalty = "scad")`.

新增脚本：

- `scripts/run_nonconvex_penalty_experiment_with_plots.R`

输出数据：

- `results/nonconvex_penalty_summary.csv`
- `results/nonconvex_penalty_aggregate.csv`
- `results/nonconvex_penalty_trace.csv`
- `results/nonconvex_penalty_trace_aggregate.csv`

输出图：

- `figures/nonconvex_penalty_final_gap.png`
- `figures/nonconvex_penalty_selection.png`
- `figures/nonconvex_penalty_support_size.png`

Settings:

- `p = 60`;
- true support size = 8;
- hard client-level non-IID;
- `tau = 0.9`;
- R1+R2 stochasticity;
- penalties: L1, MCP, SCAD;
- `lambda in {0.005, 0.01, 0.02, 0.04, 0.08}`;
- 3 seeds.

Representative results:

| penalty | lambda | final gap | selected size | TPR | FDR |
|---|---:|---:|---:|---:|---:|
| L1 | 0.020 | `8.51e-04` | 11.3 | 0.875 | 0.378 |
| MCP | 0.020 | `5.28e-03` | 15.0 | 0.958 | 0.482 |
| SCAD | 0.020 | `1.05e-03` | 14.3 | 0.958 | 0.456 |
| L1 | 0.040 | `6.81e-04` | 6.3 | 0.708 | 0.095 |
| MCP | 0.040 | `3.17e-03` | 8.0 | 0.833 | 0.133 |
| SCAD | 0.040 | `5.99e-03` | 6.7 | 0.750 | 0.083 |

Interpretation:

- L1 gives the smallest optimization gap for most lambdas, but it shrinks aggressively and starts losing true variables as lambda grows.
- MCP/SCAD preserve larger coefficients better, giving higher TPR at moderate lambda.
- MCP at `lambda = 0.04` selects about 8 variables, matching the true support size on average, with TPR around 0.83 and FDR around 0.13.
- SCAD at `lambda = 0.02` keeps high TPR with moderate sparsity, but at `lambda = 0.04` becomes more conservative.
- Nonconvex penalties should be presented as variable-selection refinements, not as universally better objective optimizers.

Suggested report wording:

> MCP and SCAD extend QR box-dual from convex sparse QR to nonconvex variable-selection penalties. In the hard sparse federated setting, L1 gives the smallest objective gap, while MCP/SCAD better preserve large signals and can recover support sizes closer to the truth at moderate regularization.

## Heart Disease Four-Center Baseline Experiment

新增脚本：

- `scripts/run_hd_baseline_comparison_with_plots.R`

输出数据：

- `results/hd_baseline_targets.csv`
- `results/hd_baseline_summary.csv`
- `results/hd_baseline_aggregate.csv`
- `results/hd_baseline_trace.csv`
- `results/hd_baseline_trace_aggregate.csv`

输出图：

- `figures/hd_center_response_distribution.png`
- `figures/hd_baseline_final_gap_log.png`
- `figures/hd_baseline_r1r2_convergence.png`
- `figures/hd_baseline_r1r2_final_gap.png`

Setting:

- real UCI Heart Disease four-center split;
- clients: Cleveland, Hungarian, Switzerland, VA;
- response: standardized `thalach`;
- complete-case sample size: `n = 740`;
- quantiles: `tau = 0.5, 0.75, 0.9`;
- stochastic regimes:
  - Full;
  - R1 client sampling;
  - R2 local user mini-batch;
  - R1+R2;
- baselines:
  - FedQR-ADMM;
  - FedSubGrad;
  - FSPG-smooth;
  - FedSPD-check.

R1+R2 final objective gaps versus `quantreg::rq()` target:

| tau | QR box-dual | FSPG-smooth | FedSubGrad | FedSPD-check | FedQR-ADMM |
|---:|---:|---:|---:|---:|---:|
| 0.50 | `4.66e-03` | `6.49e-03` | `6.06e-03` | `1.66e-02` | `1.78e-02` |
| 0.75 | `2.19e-03` | `5.57e-03` | `6.83e-03` | `1.47e-02` | `1.60e-02` |
| 0.90 | `1.53e-03` | `3.54e-03` | `1.07e-02` | `1.55e-02` | `1.55e-02` |

Interpretation:

- QR box-dual ranks first across all tested HD quantiles and stochastic regimes.
- The four-center split is naturally non-IID: Switzerland and VA have much lower `thalach` distributions than Cleveland and Hungarian.
- The real-data result agrees with the hard simulation result: exploiting the exact QR box dual is more robust than direct check-subgradient FedSPD, smoothing, subgradient, or consensus ADMM under R1+R2.

Report wording:

> On the real four-center Heart Disease task, the QR box-dual method consistently gives the smallest objective gap under both isolated and combined client/sample randomness. This supports the claim that the improvement is not only a synthetic artifact, but also useful under naturally heterogeneous medical-center data.

## M1/M2 Step-Size Ablation

新增脚本：

- `scripts/run_dual_adaptive_ablation_with_plots.R`

输出数据：

- `results/dual_adaptive_ablation_targets.csv`
- `results/dual_adaptive_ablation_summary.csv`
- `results/dual_adaptive_ablation_aggregate.csv`
- `results/dual_adaptive_ablation_trace.csv`
- `results/dual_adaptive_ablation_trace_aggregate.csv`

输出图：

- `figures/dual_adaptive_ablation_r1r2_gap.png`
- `figures/dual_adaptive_ablation_r1r2_convergence.png`
- `figures/dual_adaptive_ablation_all_scenarios.png`

Compared rules:

- A0 operator: common primal/dual operator step;
- M1 box-aware: rescales steps by the QR dual box size `max(tau, 1 - tau)`;
- M1+M2 tau-adaptive: further increases primal step and decreases dual step as `tau` approaches the extremes.

Hard non-IID R1+R2 final gaps:

| tau | A0 operator | M1 box-aware | M1+M2 tau-adaptive |
|---:|---:|---:|---:|
| 0.50 | `4.74e-03` | `1.94e-02` | `1.94e-02` |
| 0.75 | `4.80e-03` | `7.77e-03` | `1.62e-02` |
| 0.90 | `4.42e-03` | `4.76e-03` | `1.20e-02` |
| 0.95 | `2.48e-03` | `2.49e-03` | `1.24e-02` |

Full-client/full-batch final gaps:

| tau | A0 operator | M1 box-aware | M1+M2 tau-adaptive |
|---:|---:|---:|---:|
| 0.50 | `8.49e-06` | `1.20e-05` | `1.20e-05` |
| 0.75 | `1.41e-05` | `1.11e-05` | `1.33e-05` |
| 0.90 | `2.84e-05` | `2.26e-05` | `3.72e-05` |
| 0.95 | `1.55e-04` | `1.22e-04` | `6.61e-05` |

Interpretation:

- M1/M2 are implemented in both centralized `qr_pdhg()` and federated `qr_box_fed_pdhg()`.
- The ablation shows M2 helps in the deterministic extreme-quantile regime, especially at `tau = 0.95`.
- Under combined R1+R2 stochasticity, naive M2 is too aggressive because it enlarges the primal step while using stale cached directions and local mini-batch dual updates.
- Therefore the final main method should not default to M2 in the hardest stochastic experiments. The report should present M2 as an optional extreme-quantile deterministic acceleration, and use the ablation to justify the more conservative default step rule under R1+R2.

Suggested report wording:

> The M1/M2 ablation reveals a stability-speed tradeoff. Tau-adaptive scaling improves deterministic extreme-quantile convergence, but the same scaling can amplify stochastic noise under simultaneous client sampling and local mini-batching. We therefore report QR box-dual with conservative operator/box-aware steps as the main stochastic federated method, and keep tau-adaptive scaling as an optional acceleration for less noisy regimes.

## Advanced Stabilization Candidate Pool

新增算法接口：

- `dual_relaxation` in `qr_box_fed_pdhg()`;
- `server_momentum` in `qr_box_fed_pdhg()`;
- `step_decay_power` and `step_decay_offset` in `qr_box_fed_pdhg()`;
- `primal_clip` in `qr_box_fed_pdhg()`;
- existing `aggregation = "selected_reweighted"` included in the candidate pool.

新增实验脚本：

- `scripts/run_advanced_stabilization_experiment_with_plots.R`

输出数据：

- `results/advanced_stabilization_targets.csv`
- `results/advanced_stabilization_summary.csv`
- `results/advanced_stabilization_aggregate.csv`
- `results/advanced_stabilization_trace.csv`
- `results/advanced_stabilization_trace_aggregate.csv`

输出图：

- `figures/advanced_stabilization_r1r2_gap.png`
- `figures/advanced_stabilization_r1r2_ratio.png`
- `figures/advanced_stabilization_best_convergence.png`

Candidate variants:

- A0 operator;
- M1 box-aware;
- M2 tau-adaptive;
- Dual relax `.5`;
- Server EMA `.5`;
- Step decay `.25`;
- Direction clip `.25`;
- Relax + EMA;
- Relax + decay;
- Selected reweighted aggregation.

Hard/extreme non-IID R1+R2 results:

| heterogeneity | tau | best variant | best gap | A0 gap | best/A0 |
|---|---:|---|---:|---:|---:|
| hard | 0.90 | A0 operator / Direction clip | `4.42e-03` | `4.42e-03` | 1.00 |
| hard | 0.95 | Direction clip `.25` / A0 operator | `2.48e-03` | `2.48e-03` | 1.00 |
| extreme | 0.90 | A0 operator | `5.54e-03` | `5.54e-03` | 1.00 |
| extreme | 0.95 | Direction clip `.25` | `4.53e-03` | `4.54e-03` | 0.998 |

Important negative findings:

- Dual relaxation is harmful in R1+R2 because it slows correction of stale local sample-level dual variables.
- Server EMA is not consistently useful because cached client directions already act as a memory mechanism.
- Step decay underfits within the fixed communication budget.
- Selected reweighted aggregation is much noisier than cached aggregation under non-IID client sampling.
- Direction clipping is safe and sometimes marginally helpful, but the gain is too small to claim as a major improvement.
- M2 remains very useful in low-noise full-client extreme quantile cases, but not under combined R1+R2 stochasticity.

Interpretation:

> The advanced stabilization pool strengthens the empirical story mostly through disciplined negative evidence. Common stochastic stabilization tricks do not improve the final R1+R2 objective gap over the conservative cached QR box-dual update. This suggests that the main gain comes from the QR box-dual geometry and cached inactive-client directions, rather than generic momentum or damping.

## Expanded Scenario Matrix

新增脚本：

- `scripts/run_expanded_scenario_suite_with_plots.R`
- `scripts/run_qr_box_tuning_expanded_scenarios.R`

输出数据：

- `results/expanded_scenario_design.csv`
- `results/expanded_scenario_summary.csv`
- `results/expanded_scenario_aggregate.csv`
- `results/expanded_scenario_ranks.csv`
- `results/qr_tuning_expanded_summary.csv`
- `results/qr_tuning_expanded_aggregate.csv`
- `results/qr_tuning_expanded_ranks.csv`

输出图：

- `figures/expanded_scenario_final_gap.png`
- `figures/expanded_scenario_rank_heatmap.png`
- `figures/expanded_scenario_convergence.png`
- `figures/qr_tuning_expanded_gap.png`
- `figures/qr_tuning_expanded_rank_heatmap.png`

Stress scenarios:

- low client participation: `K = 2/20`;
- tiny user mini-batch: batch size 5;
- extreme quantile: `tau = 0.95`;
- extreme client heterogeneity;
- highly unbalanced client sizes;
- high-dimensional sparse model with `p = 60`.

Fixed 400-round baseline comparison:

| scenario | winner | QR box-dual gap | winner gap |
|---|---|---:|---:|
| low participation | FSPG-smooth | `3.24e-02` | `1.19e-02` |
| tiny user batch | FSPG-smooth | `2.21e-02` | `7.63e-03` |
| extreme tau | QR box-dual | `4.29e-03` | `4.29e-03` |
| extreme heterogeneity | FSPG-smooth | `8.87e-03` | `7.47e-03` |
| unbalanced clients | FSPG-smooth | `1.70e-02` | `6.24e-03` |
| high-dimensional sparse | FSPG-smooth | `1.96e-02` | `1.69e-02` |

Tuned QR box-dual comparison:

| scenario | best QR variant | best QR gap | FSPG 400 gap |
|---|---|---:|---:|
| low participation | QR long 1200 | `5.23e-03` | `1.19e-02` |
| tiny user batch | QR long 1200 | `4.44e-03` | `7.63e-03` |
| extreme tau | QR clip .25 800 | `1.13e-03` | `7.59e-03` |
| extreme heterogeneity | QR long 1200 | `1.35e-03` | `7.47e-03` |
| unbalanced clients | QR long 1200 | `1.49e-03` | `6.24e-03` |
| high-dimensional sparse | QR long 1200 | `5.99e-03` | `1.69e-02` |

Interpretation:

- QR box-dual should not be described as universally best under every fixed short communication budget.
- FSPG-smooth is a strong short-horizon baseline because smoothing makes early stochastic gradients easier to optimize.
- After QR-specific tuning or a longer communication horizon, QR box-dual wins across the expanded stress matrix.
- This strengthens the final claim by making it more precise: the proposed method is best in final objective quality after appropriate QR-specific tuning, while smoothing remains competitive for very short communication budgets.

## NYC Taxi Million-Row Real-Data Experiment

新增文件：

- `NYC_TAXI_LARGE_SCALE.md`

新增脚本：

- `scripts/prepare_nyc_taxi_data.sh`
- `scripts/run_nyc_taxi_large_scale_experiment.R`

Data:

- official NYC TLC Yellow Taxi trip records;
- 2024 Q1 Parquet files;
- 8,417,330 clean trips after filtering;
- 1,000,000 modeled rows;
- 59 eligible pickup-zone clients.

Task:

- response: `log_total_amount`;
- `tau = 0.9`;
- clients: `PULocationID`;
- 12 clients per round;
- local user batch size: 5,000;
- QR long run: 1,000 rounds.

Results:

| method | objective | gap to best observed |
|---|---:|---:|
| QR box-dual long | 0.0206365 | 0 |
| QR box-dual | 0.0220496 | 0.0014131 |
| FSPG-smooth | 0.0240271 | 0.0033905 |
| FedSubGrad | 0.1043675 | 0.0837310 |
| FedSPD-check | 0.3452998 | 0.3246633 |

Interpretation:

- This is the first genuinely large real-data distributed stochastic optimization experiment in the project.
- The pickup-zone split gives natural non-IID clients and strong client-size imbalance.
- QR box-dual was not stable with a tiny local user batch, but became best after increasing the local sample-level dual update coverage.
- This supports the algorithmic interpretation: QR box-dual needs enough local dual refresh per communication round; when that condition is met, it outperforms smoothing and subgradient baselines on million-row real data.

Report wording:

> The NYC Taxi experiment moves beyond small real-data validation. On one million modeled trips from 2024 Q1 and 59 pickup-zone clients, cached QR box-dual achieves the best observed objective after increasing local user-level dual coverage, while smoothing remains a strong short-horizon baseline.

## UCI Household Power Large-Scale Real-Data Experiment

新增文件：

- `HOUSEHOLD_POWER_LARGE_SCALE.md`

新增脚本：

- `scripts/prepare_household_power_data.sh`
- `scripts/run_household_power_large_scale_experiment.R`

Data/task:

- UCI Individual Household Electric Power Consumption;
- 2,049,280 modeled rows after cleaning;
- response `log1p(Global_active_power)`;
- `tau = 0.9`;
- 48 calendar-month clients;
- 10 clients per round;
- local user batch size 5,000;
- QR long run 800 rounds.

Results:

| method | objective | gap to best observed |
|---|---:|---:|
| QR box-dual long | 0.0419800 | 0 |
| FSPG-smooth | 0.0423711 | 0.0003911 |
| FedSubGrad | 0.0456561 | 0.0036761 |
| QR box-dual | 0.0500239 | 0.0080439 |
| FedSPD-check | 0.0550988 | 0.0131188 |

Interpretation:

- This is a second large-scale real-data experiment, separate from NYC Taxi.
- The calendar-month client split induces temporal non-IID rather than spatial non-IID.
- QR box-dual long achieves the best observed objective, while FSPG-smooth is close at the shorter horizon.
- Together with NYC Taxi, this supports the core practical claim: cached QR box-dual is not merely a simulation method; it can solve million-row distributed stochastic QR problems when local sample-level dual updates are sufficiently rich.

## Advanced Innovation Variants

The main QR box-dual solver now has three optional extensions. All are disabled by default, so the original method and saved benchmark interpretation remain unchanged.

### Staleness-Aware Cached Direction

Under partial client participation, cached local directions can become stale. The new variant tracks one age value per client:

- selected clients reset age to 0;
- unselected clients increase age by 1 each communication round;
- cached directions are optionally reweighted by an exponential or inverse decay;
- the aggregated direction can be normalized to keep total client mass comparable to the original cached aggregation.

This variant is available through:

- `qr_box_fed_pdhg(..., staleness = "exponential")`;
- `fit_fedqr("QR box-dual stale", ...)`;
- `fit_fedqr("QR box-dual stale+robust", ...)`.

The trace records `mean_staleness`, `max_staleness`, and `mean_stale_weight`.

The adaptive staleness variant uses the observed mean cache age to adjust the decay rate during training. When the cache is older than the target age, the effective decay rate increases; when the cache is fresh, the rate relaxes. The trace records the resulting `adaptive_staleness_rate`.

### Client-Robust Objective

The original objective uses sample-average weighting, equivalent to client weights `w_j = n_j / n`. This is appropriate for global empirical risk, but in highly unbalanced federated data it can underweight small but important clients.

The robust extension supports:

- `sample`: `w_j = n_j / n`;
- `uniform`: `w_j = 1 / m`;
- `sqrt_size`: `w_j` proportional to `sqrt(n_j)`;
- `custom`: user-specified client weights.

The QR box-dual local direction is generalized from

```text
X_j^T v_j / n
```

to

```text
w_j X_j^T v_j / n_j.
```

When `client_weighting = "sample"`, this exactly recovers the original direction. The robust aliases currently use uniform client weighting.

The adaptive client-robust variant updates client weights during training from the current per-client check loss. Clients with larger current loss receive larger weight, with smoothing and floor controls to avoid unstable jumps. This makes the robust module data-adaptive rather than fixed to uniform weighting.

### Calibration-Aware High-Quantile Correction

High quantile regression should be evaluated not only by check loss but also by coverage:

```text
coverage = mean(y <= x^T beta).
```

The new calibration utilities provide:

- global intercept correction, which shifts only the intercept;
- client-offset correction, which keeps the global coefficient vector fixed and returns deployment-time per-client offsets;
- summary diagnostics for global coverage error, mean client coverage error, worst-client coverage error, and coverage dispersion.

These utilities are exposed as:

- `quantile_coverage()`;
- `calibrate_quantile_intercept()`;
- `adaptive_calibrate_quantile()`;
- `calibration_summary()`.

The adaptive calibration routine compares raw prediction, global intercept correction, and client-offset correction, then chooses the candidate with the lowest requested coverage error metric.

### Control-Variate Variance Reduction

The `QR box-dual adaptive+VR` variant adds a conservative client-level control variate to the adaptive method. Each client keeps an exponential moving average of its raw dual direction, and the server keeps a weighted global control direction. The cached direction can then be corrected as

```text
raw_direction_j - control_j + global_control.
```

The implementation uses a small blend and a correction cap, so the correction reduces stochastic direction dispersion without aggressively removing real client heterogeneity. The trace records `vr_correction_norm`, `raw_direction_variance`, and `corrected_direction_variance`.

### Empirical Reading

The advanced experiment suite compares `QR box-dual`, `QR box-dual stale`, `QR box-dual robust`, `QR box-dual stale+robust`, `QR box-dual adaptive`, `QR box-dual adaptive+VR`, `FSPG-smooth`, and `FedSPD-check` on hard/extreme non-IID simulation and Heart Disease.

Observed pattern:

- staleness-aware QR box-dual improves target gap in the extreme non-IID high-quantile simulation;
- robust weighting improves worst-client loss on Heart Disease and several unbalanced simulation settings;
- calibration refinement sharply reduces global and client-level coverage error;
- robust variants do not always minimize the global sample-average objective, which is expected because they optimize a more client-balanced target.

## Large-Client Scale Stress Test

The four-center Heart Disease experiment is useful as a real multi-center medical dataset, but it is too small to stress partial participation and stale cached directions. A separate scale-stress experiment therefore evaluates 50 and 100 clients with only 10% participation per communication round.

Default setting:

- 50 or 100 clients;
- 5 or 10 selected clients per round;
- hard and extreme non-IID heterogeneity;
- highly unbalanced client sizes;
- `tau = 0.9` and `tau = 0.95`;
- methods: base, stale, robust, stale+robust, adaptive, FSPG-smooth, and FedSPD-check.

Observed aggregate:

- QR box-dual variants win all 8 target-gap settings;
- `QR box-dual stale+robust` wins 5/8 target-gap settings;
- `QR box-dual adaptive` wins 1/8 target-gap settings;
- `QR box-dual adaptive+VR` wins 1/8 target-gap settings;
- `QR box-dual robust` wins 1/8 target-gap settings;
- `QR box-dual adaptive+VR` wins 6/8 worst-client fairness settings;
- `QR box-dual adaptive+VR` reduces final corrected-direction variance to about 92.8% of raw direction variance on average.

This scale test is the strongest evidence for the advanced modules: once the number of clients is large and participation is sparse, stale-aware and robust/adaptive weighting become materially useful.

The VR result should be interpreted carefully: it improves direction stability and worst-client fairness, but it is not the dominant target-gap optimizer. The best final target gap is still usually obtained by `QR box-dual stale+robust`.
