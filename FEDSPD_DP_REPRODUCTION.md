# FedSPD-DP 论文复现记录

## 目标论文

当前复现对象：

- Yiwei Li, Weijian Wang, Tianqing Zhu, Jingpeng Li, Bo Chang, and Shui Yu.
  "Federated Stochastic Primal-dual Learning with Differential Privacy."
  arXiv:2204.12284.

本地文件：

- `references/fedspd_dp_2204.12284.pdf`
- `references/fedspd_src/FedSPD-DP_arxiv.tex`

## 当前实现

核心实现：

- `R/fedspd_dp_qr.R`

复现实验脚本：

- `scripts/run_fedspd_dp_reproduction.R`

输出结果：

- `results/fedspd_dp_reproduction_summary.csv`
- `results/fedspd_dp_reproduction_trace.csv`

## Algorithm 1 对齐情况

已实现论文 Algorithm 1 的关键状态：

- server consensus model: `x0`
- client local primal models: `x_i`
- client consensus dual variables: `lambda_i`
- cached uploaded local models: `y_tilde_i`
- partial client participation: each round samples `K` clients without replacement
- inactive clients: keep cached `x_i`, `lambda_i`, and `y_tilde_i`
- local `Q`-step proximal stochastic update
- server aggregation using all cached uploads

已实现论文核心更新：

```text
x0^t = (1/N) sum_i ytilde_i
```

```text
x_i^{t,r} =
prox_{R_i / (gamma_i^t + rho)}
(
  (gamma_i^t x_i^{t,r-1}
   + rho x0^t
   + lambda_i^{t-1}
   - grad f_i(x_i^{t,r-1}; B_i^{t,r}))
  / (gamma_i^t + rho)
)
```

```text
x_i^t = average of Q local inner iterates
```

```text
lambda_i^t = lambda_i^{t-1} + rho (x0^t - x_i^t)
```

```text
ytilde_i^t = x_i^t - lambda_i^t / rho + noise
```

## QR 适配

FedSPD-DP 论文假设本地 `f_i` 可微，而原始 check loss 不可微。为了忠实保留论文算法假设，当前复现版使用 smoothed quantile loss：

```text
rho_{tau,mu}(u) = tau * u + mu * log(1 + exp(-u / mu))
```

其对 `beta` 的梯度为：

```text
grad_beta = - X^T [tau - sigmoid(-residual / mu)] / batch_size
```

这样做的含义：

- `fedspd_dp_qr()` 是 FedSPD-DP Algorithm 1 的 QR-smoothed adaptation；
- 之前的 `fed_qr_spd()` 是 QR box-geometry Fenchel/PDHG 方法；
- 二者应该分开叙述和比较。

当前也已实现原始 check loss 的 subgradient 模式：

```text
loss = "check"
```

对应局部 stochastic subgradient：

```text
psi = tau - 1{y - x beta < 0}
grad_beta = - X^T psi / batch_size
```

因此 `fedspd_dp_qr()` 现在支持两类 QR 适配：

- `loss = "smooth"`：忠实满足 FedSPD-DP 论文中 `f_i` 可微的假设；
- `loss = "check"`：直接优化原始非光滑 QR check loss，属于 FedSPD-DP 框架下的 nonsmooth QR adaptation。

## DP 与步长

已实现 DP 噪声开关：

- `dp = FALSE`：不加 DP 噪声，用于先验证优化结构；
- `dp = TRUE`：按论文 sensitivity 公式加入 Gaussian noise。

已实现 DP noise scale：

```text
Q = 1:
sigma_{i,t} = 4G sqrt(2 log(1.25 / delta)) / [epsilon (rho + gamma_i^t)]

Q > 1:
sigma_{i,t} =
4QG sqrt(2 log(1.25 / delta)) /
[(Q - 1) epsilon (rho + gamma_i^t)]
```

已实现三种 `gamma`：

- `gamma_rule = "paper"`：使用论文 Theorem 2 的 `gamma_i^t`；
- `gamma_rule = "sqrt"`：实验用简化 `gamma0 * sqrt(t)`；
- `gamma_rule = "constant"`：固定 `gamma0`。

论文 Theorem 2 的步长公式已实现为 `fedspd_gamma()`。

## 第一轮复现结果

设置：

- 模拟 QR 数据；
- `n = 800, p = 10`;
- 20 clients;
- `tau = 0.5`;
- smoothed QR with `mu = 0.05`;
- `rho = 20`;
- `rounds = 200`;
- batch size = 30。

结果：

| label | K | Q | DP | gamma | final QR objective | consensus gap | beta L2 error |
|---|---:|---:|---|---|---:|---:|---:|
| full clients | 20 | 5 | no | sqrt | 0.331485 | 0.00246 | 0.1167 |
| half clients | 10 | 5 | no | sqrt | 0.382602 | 0.00582 | 0.4715 |
| quarter clients | 5 | 5 | no | sqrt | 0.499605 | 0.00863 | 0.9078 |
| half clients | 10 | 1 | no | sqrt | 0.485047 | 0.00540 | 0.8594 |
| half clients | 10 | 10 | no | sqrt | 0.358359 | 0.00688 | 0.3477 |
| half clients | 10 | 5 | yes | sqrt | 0.409988 | 0.04890 | 0.5970 |
| full clients | 20 | 5 | no | paper | 0.341009 | 0.00248 | 0.2275 |
| central QR-PDHG | NA | NA | no | NA | 0.329151 | NA | 0.1112 |

Interpretation:

- Full-client FedSPD-DP is close to centralized QR-PDHG, which validates the consensus implementation.
- Smaller `K` worsens performance, matching the paper's PCP tradeoff.
- Larger `Q` improves communication-round efficiency in this setting, matching the paper's local-SGD motivation.
- DP noise worsens accuracy and consensus, matching the paper's privacy-utility tradeoff.
- The paper gamma schedule is conservative on this QR-smoothed problem, but still converges in the right direction.

## Remaining Work for a Strong Reproduction

The current reproduction implements Algorithm 1 faithfully, but it is not yet a full empirical reproduction of every paper figure.

Remaining work:

1. Add paper-style comparison baselines:
   - FedAvg / DP-FedAvg;
   - DP-ADMM-like baseline if time allows.
2. Reproduce K-sweep curves:
   - objective vs communication round for different `K`.
3. Reproduce Q-sweep curves:
   - objective vs communication round for different local steps `Q`.
4. Reproduce DP tradeoff curves:
   - different `epsilon` or total privacy budgets.
5. Add plots from the CSV traces.
6. Then compare this paper-faithful FedSPD-DP baseline with QR box-PDHG.

## Nonsmooth QR Adaptation Result

新增脚本：

- `scripts/run_qr_adaptation_comparison.R`

输出：

- `results/qr_adaptation_comparison_summary.csv`
- `results/qr_adaptation_comparison_trace.csv`

设置：

- `tau = 0.5, 0.75, 0.9`;
- 20 clients;
- full clients vs half clients;
- compare `loss = "smooth"` and `loss = "check"`;
- no DP noise.

结果摘要：

| tau | loss | clients | final QR objective | consensus gap | beta L2 error |
|---:|---|---|---:|---:|---:|
| 0.50 | smooth | full | 0.329956 | 0.00217 | 0.0936 |
| 0.50 | check | full | 0.329961 | 0.00224 | 0.1005 |
| 0.50 | central | central | 0.329151 | NA | 0.1112 |
| 0.75 | smooth | full | 0.326003 | 0.00215 | 0.2120 |
| 0.75 | check | full | 0.326317 | 0.00220 | 0.2186 |
| 0.75 | central | central | 0.321028 | NA | 0.2361 |
| 0.90 | smooth | full | 0.214195 | 0.00215 | 0.4862 |
| 0.90 | check | full | 0.214369 | 0.00212 | 0.4914 |
| 0.90 | central | central | 0.200561 | NA | 0.3062 |

Interpretation:

- The nonsmooth `check` mode is stable inside the paper-faithful FedSPD-DP consensus framework.
- `smooth` and `check` behave similarly under the current step settings.
- Full participation remains much closer to centralized QR than partial participation.
- The gap at `tau = 0.9` leaves room for the QR-specific box-dual method and step-size refinements.

Next:

- Implement FedSPD-DP-box QR or a clean hybrid comparison where:
  - `fedspd_dp_qr(loss = "smooth")` is the paper-assumption baseline;
  - `fedspd_dp_qr(loss = "check")` is the direct nonsmooth QR adaptation;
  - `fed_qr_spd()` / next box version is the QR geometry-specialized method.
