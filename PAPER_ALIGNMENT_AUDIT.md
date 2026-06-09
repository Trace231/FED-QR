# Paper Alignment Audit

## Reference Checked

The project summary cites "FedSPD-DP '22" as the primal-dual federated reference. The matching paper is:

- Yiwei Li, Weijian Wang, Tianqing Zhu, Jingpeng Li, Bo Chang, and Shui Yu.
  "Federated Stochastic Primal-dual Learning with Differential Privacy."
  arXiv:2204.12284.

Local copies:

- `references/fedspd_dp_2204.12284.pdf`
- `references/fedspd_src/FedSPD-DP_arxiv.tex`

## Important Finding

The current `rfedqr` implementation is **not** a faithful implementation of FedSPD-DP Algorithm 1.

It is a QR-specialized Fenchel/PDHG prototype:

- sample-level dual variables for the quantile check loss;
- box projection `clip[tau - 1, tau]`;
- a global primal vector `beta`;
- client sampling and mini-batching layered onto that update.

FedSPD-DP Algorithm 1 has a different state and update structure:

- client-local primal variables `x_i`;
- a server global consensus variable `x_0`;
- client-level dual variables `lambda_i` for the consensus constraints `x_i = x_0`;
- local proximal SGD inner loop with `Q` steps;
- inactive clients keep cached `x_i`, `lambda_i`, and uploaded model `y_i`;
- server aggregates cached uploaded models from all clients;
- optional DP Gaussian perturbation on uploaded local model.

Therefore, the current implementation should not be described as "a faithful reproduction of FedSPD-DP".

## What Current Code Matches

Current code matches the QR-specific saddle-point idea:

```text
min_beta mean_i rho_tau(y_i - x_i^T beta) + P_lambda(beta)
```

Using the Fenchel conjugate of the check loss:

```text
rho_tau(u) = sup_{v in [tau - 1, tau]} v u
```

This yields a primal-dual update with:

```text
v_i <- clip(v_i + sigma * (y_i - x_i^T beta_bar), tau - 1, tau)
beta <- prox_{eta P}(beta + eta * X^T v / n)
```

That part is mathematically meaningful for quantile regression, and the centralized version has been checked against `quantreg::rq()`.

## What Current Code Does Not Match

Compared with FedSPD-DP Algorithm 1, current `fed_qr_spd()` lacks:

1. Explicit local primal models `x_i`.
2. Explicit server consensus model `x_0`.
3. Consensus dual variables `lambda_i`.
4. Local `Q`-step proximal SGD inner loop.
5. FedSPD-DP aggregation:

   ```text
   x_0^t = (1/N) sum_i ytilde_i^t
   ytilde_i^t = x_i^t - lambda_i^t / rho + noise
   ```

6. Cached inactive-client uploads.
7. FedSPD-DP time-varying inverse step `gamma_i^t`.
8. DP sensitivity/noise accounting.

## Why This Matters

The current implementation is useful as a first QR-PDHG sanity check, but it is too loose to support a claim such as:

> "We faithfully implement the FedSPD-DP paper for federated quantile regression."

A faithful project should instead say one of the following:

1. **FedSPD-DP-faithful route**:
   We implement Algorithm 1 from FedSPD-DP and adapt the local loss to quantile regression through a smoothing or proximal surrogate.

2. **QR-PDHG route**:
   We implement a QR-specific Fenchel/PDHG algorithm, using FedSPD-DP only as motivation for partial client participation and local stochasticity.

3. **Hybrid route**:
   We implement both:
   - a FedSPD-DP-faithful baseline;
   - a QR-specific box-PDHG method.

The hybrid route is the strongest for a course project because it prevents the main method from being only a toy and makes the algorithmic comparison cleaner.

## Recommended Correction

To make the project non-toy and paper-grounded, the next implementation should add a separate function:

```text
fedspd_consensus_qr()
```

It should follow FedSPD-DP Algorithm 1:

1. Maintain `x0`, `x_i`, `lambda_i`, and cached uploaded `ytilde_i`.
2. At each communication round:
   - server computes `x0 = mean_i ytilde_i`;
   - server samples `K` clients without replacement;
   - each active client runs `Q` local proximal stochastic updates;
   - active client updates `lambda_i`;
   - active client uploads `ytilde_i = x_i - lambda_i / rho` without DP noise initially;
   - inactive clients keep cached states.
3. For quantile regression:
   - start with smoothed quantile loss to satisfy the paper's differentiable `f_i` assumption;
   - then optionally add a nonsmooth QR-specific variant using subgradient or Fenchel splitting.
4. Compare this paper-faithful implementation against the current QR box-PDHG implementation.

Status: implemented as `fedspd_dp_qr()` in `R/fedspd_dp_qr.R`, with the first reproduction script in `scripts/run_fedspd_dp_reproduction.R`.

## Wording to Use in the Report

Until the consensus implementation is added, use:

> We first implemented a QR-specific Fenchel-PDHG prototype and validated it against `quantreg`. After auditing FedSPD-DP Algorithm 1, we separated this prototype from the paper-faithful consensus primal-dual implementation, which maintains local primal variables, consensus dual variables, partial client participation, and cached inactive-client uploads.

Avoid:

> We faithfully reproduce FedSPD-DP.
