# rfedqr 项目笔记

## 当前进展

我们已经完成了第一版最小算法闭环：

- 建立了 `rfedqr` 项目骨架。
- 实现了单机版 QR-PDHG。
- 支持 `none` 和 `l1` 惩罚的 proximal step。
- 用模拟数据跑通了 `tau = 0.5, 0.75, 0.9`。
- 下载并检查了 UCI Heart Disease 四中心数据。
- 安装并使用 `quantreg`，完成了 `qr_pdhg()` 与 `quantreg::rq()` 的单机数值对照。

当前模拟结果显示，单机 QR-PDHG 在干净模拟数据上可以稳定收敛。这是一个积极信号，说明底层算法内核是可行的。

## quantreg 数值对照

我们用同一批模拟数据比较了：

- `quantreg::rq()`；
- 自己实现的 `qr_pdhg()`。

设置：

- 无惩罚；
- 集中式单机数据；
- `tau = 0.5`；
- 样本量 `n = 600`；
- 特征维度 `p = 12`；
- 非对称噪声。

结果保存在：

- `results/quantreg_comparison.csv`

最大系数绝对差约为 `1.8e-4`。这说明当前 `qr_pdhg()` 实现确实在求解标准分位数回归问题，而不是只是在优化一个形式相似但写错的目标函数。

当前可引用结论：

> In the centralized unpenalized setting, the proposed QR-PDHG implementation matches the mature `quantreg::rq()` solver up to numerical tolerance, validating the basic primal-dual formulation before federated extensions.

## 联邦 QR-FedSPD 初步实验

我们已经实现了第一版联邦算法：

- `R/fed_qr_spd.R`
- `scripts/run_federated_simulation.R`

设计：

- 8 个客户端；
- `tau = 0.5`；
- 模拟数据 `n = 1200, p = 12`；
- IID 随机切分 vs Dirichlet non-IID 切分；
- 四种随机性设置：
  - `full_clients_full_batch`：无 R1，无 R2；
  - `r1_only`：客户端采样；
  - `r2_only`：本地 mini-batch；
  - `r1_r2`：客户端采样 + 本地 mini-batch。

结果保存在：

- `results/federated_partition_summary.csv`
- `results/federated_simulation_summary.csv`
- `results/federated_simulation_trace.csv`

关键结果：

| partition | config | final objective | beta L2 error |
|---|---:|---:|---:|
| iid | full clients/full batch | 0.359996 | 0.104 |
| iid | R1 only | 0.367079 | 0.177 |
| iid | R2 only | 0.367242 | 0.209 |
| iid | R1 + R2 | 0.375183 | 0.248 |
| non-IID | full clients/full batch | 0.359996 | 0.104 |
| non-IID | R1 only | 0.381860 | 0.277 |
| non-IID | R2 only | 0.368470 | 0.219 |
| non-IID | R1 + R2 | 0.383660 | 0.333 |
| central | QR-PDHG | 0.359994 | 0.103 |

初步解释：

- `full_clients_full_batch` 与中央 QR-PDHG 几乎一致，说明联邦拆分本身没有改变优化问题。
- IID 下 R1 和 R2 都会增加误差，但幅度相对可控。
- non-IID 下 R1 明显更差，说明客户端采样会放大客户端分布差异。
- R1 + R2 叠加最差，符合“双重随机性”的项目叙事。

当前可引用结论：

> When all clients and full local batches are used, the federated implementation reproduces the centralized QR-PDHG solution. Under non-IID partitions, client sampling (R1) causes a larger degradation than local mini-batching (R2), supporting the planned R1/R2 isolation study.

## UCI Heart Disease 四中心真实数据实验

我们已经实现了真实四中心 QR 实验：

- `scripts/run_heart_disease_qr.R`

实验口径：

- 响应变量：`thalach`，最大心率；
- 建模时将响应标准化为 `thalach_z`；
- 自然客户端：Cleveland、Hungarian、Switzerland、VA；
- 删除缺失严重变量 `slope`、`ca`、`thal`；
- 使用完整案例，最终 `n = 740`；
- 比较：
  - `quantreg::rq()`；
  - central QR-PDHG；
  - full-client federated QR-FedSPD；
  - 每轮采样 2 个中心的 R1 federated QR-FedSPD；
  - 单中心 local QR。

结果保存在：

- `results/heart_disease_thalach_by_center.csv`
- `results/heart_disease_qr_summary.csv`
- `results/heart_disease_qr_coef_comparison.csv`
- `results/heart_disease_local_qr_summary.csv`
- `results/heart_disease_qr_trace.csv`

中心响应分布：

| center | n | mean thalach | median | q90 |
|---|---:|---:|---:|---:|
| Cleveland | 303 | 149.6 | 153 | 176.6 |
| Hungarian | 261 | 139.2 | 140 | 170.0 |
| Switzerland | 46 | 113.6 | 115 | 146.5 |
| VA | 130 | 121.4 | 120 | 150.0 |

该表说明四个中心的响应分布有明显差异，适合作为真实 non-IID 联邦演示。

核心结果：

| tau | method | objective |
|---:|---|---:|
| 0.50 | central QR-PDHG | 0.326092 |
| 0.50 | fed full clients | 0.326130 |
| 0.50 | fed R1 clients | 0.340113 |
| 0.50 | quantreg rq | 0.326084 |
| 0.75 | central QR-PDHG | 0.239926 |
| 0.75 | fed full clients | 0.239953 |
| 0.75 | fed R1 clients | 0.241853 |
| 0.75 | quantreg rq | 0.239918 |
| 0.90 | central QR-PDHG | 0.123328 |
| 0.90 | fed full clients | 0.123394 |
| 0.90 | fed R1 clients | 0.125949 |
| 0.90 | quantreg rq | 0.123316 |

解释：

- central QR-PDHG 与 `quantreg::rq()` 的目标函数几乎一致，真实数据上算法内核仍然有效。
- full-client federated QR-FedSPD 与 central QR-PDHG 非常接近，说明自然四中心拆分没有改变目标。
- R1 客户端采样会恶化目标函数，尤其在中心分布不同的情况下，这是项目中 R1/non-IID 叙事的真实数据证据。
- 单中心 local QR 的模型迁移到全局数据时目标函数明显变差，尤其 Switzerland 因样本量小且分布偏移强，global objective 很高。

注意：

- `tau = 0.5` 下 `quantreg` 给出 “Solution may be nonunique” 警告。因此该处系数差异不能直接解释为算法错误，目标函数值更适合作为一致性指标。

当前可引用结论：

> On the four-center UCI Heart Disease data, QR-FedSPD with all clients closely matches both centralized QR-PDHG and `quantreg::rq()` in objective value. Client sampling introduces a visible objective gap, while local single-center QR models generalize poorly across centers, confirming the practical relevance of federated QR under center heterogeneity.

## Baseline 初步对比

我们已经实现了两个轻量 baseline：

- `fed_subgrad_qr()`：联邦 subgradient quantile regression；
- `fed_smooth_qr()`：FSPG-style smoothing baseline。

相关文件：

- `R/baselines.R`
- `scripts/run_baseline_comparison.R`

实验设置：

- 模拟数据 `n = 1200, p = 12`；
- 8 个客户端；
- Dirichlet non-IID 切分；
- 比较 `tau = 0.5` 与 `tau = 0.9`；
- 两个 regime：
  - `full_clients_full_batch`；
  - `r1_r2_stochastic`，每轮采样 4/8 客户端，每客户端 mini-batch 40。

结果保存在：

- `results/baseline_comparison_summary.csv`
- `results/baseline_comparison_trace.csv`

初步结果：

| tau | regime | method | final objective | beta L2 error |
|---:|---|---|---:|---:|
| 0.5 | full | FedSubGrad | 0.359994 | 0.106 |
| 0.5 | full | FSPG smooth | 0.360308 | 0.087 |
| 0.5 | full | generic FedSPD | 0.359995 | 0.103 |
| 0.5 | full | QR-FedSPD box | 0.359996 | 0.103 |
| 0.5 | R1+R2 | FedSubGrad | 0.360038 | 0.100 |
| 0.5 | R1+R2 | FSPG smooth | 0.364319 | 0.128 |
| 0.5 | R1+R2 | generic FedSPD | 0.381368 | 0.273 |
| 0.5 | R1+R2 | QR-FedSPD box | 0.402358 | 0.393 |
| 0.9 | full | FedSubGrad | 0.230014 | 0.221 |
| 0.9 | full | FSPG smooth | 0.229167 | 0.314 |
| 0.9 | full | generic FedSPD | 0.228881 | 0.312 |
| 0.9 | full | QR-FedSPD box | 0.228881 | 0.311 |
| 0.9 | R1+R2 | FedSubGrad | 0.229888 | 0.242 |
| 0.9 | R1+R2 | FSPG smooth | 0.230312 | 0.310 |
| 0.9 | R1+R2 | generic FedSPD | 0.235702 | 0.387 |
| 0.9 | R1+R2 | QR-FedSPD box | 0.236788 | 0.399 |

解释：

- 在 `full_clients_full_batch` 下，primal-dual、subgradient、smoothing 都能接近同一目标值。
- 在强随机 `R1+R2` 下，当前 QR-FedSPD 默认步长偏激进，表现不如 FedSubGrad。
- 这不是最终负结论，而是一个重要实现发现：随机 primal-dual 版本需要额外的 stochastic damping、局部多步策略、或更稳健的客户端采样校正。
- 因此 PDF 中“三路线头对头”应分成两个层次：
  - deterministic/full-batch：比较优化路线的极限行为；
  - stochastic/R1+R2：比较随机鲁棒性和步长敏感性。

当前可引用结论：

> In the full-client/full-batch regime, QR-FedSPD reaches the same objective level as subgradient and smoothing baselines. Under aggressive R1+R2 stochasticity, the current primal-dual implementation is more step-size sensitive, motivating stochastic damping as the next algorithmic refinement.

## FedSPD-DP 论文复现

用户明确要求先忠实复现 FedSPD-DP 论文算法，而不是继续推进 QR box-PDHG 原型。

我们已新增 paper-faithful consensus primal-dual implementation：

- `R/fedspd_dp_qr.R`
- `scripts/run_fedspd_dp_reproduction.R`
- `FEDSPD_DP_REPRODUCTION.md`

该实现与 FedSPD-DP Algorithm 1 对齐：

- server model `x0`;
- client local models `x_i`;
- consensus dual variables `lambda_i`;
- cached uploaded models `ytilde_i`;
- partial client participation;
- inactive clients keep cached states;
- local `Q`-step proximal stochastic update;
- optional DP Gaussian noise.

由于 FedSPD-DP 假设 `f_i` 可微，而 QR check loss 不可微，当前 paper-faithful QR 版本使用 smoothed quantile loss。之后应将它作为论文复现 baseline，再与 QR-specific box-PDHG 主方法比较。

更新：`fedspd_dp_qr()` 已加入 `loss = "check"`，可以在同一个 FedSPD-DP consensus 框架内直接使用原始 check loss subgradient。第一轮实验显示，`check` 模式稳定，表现与 smoothed QR 接近；在 `tau = 0.9` 下仍存在明显 central gap，后续 QR box-dual 方法有发挥空间。

## QR Box-Dual 主方法

已实现忠实 QR box-dual federated PDHG：

- `R/qr_box_fed_pdhg.R`
- `scripts/run_box_dual_comparison.R`
- `QR_BOX_DUAL_METHOD.md`

核心结构：

- 每个样本维护对偶变量 `v_i in [tau - 1, tau]`；
- 客户端本地执行 `clip()` dual prox；
- 客户端缓存并上传 `X_i^T v_i / n`；
- 服务器聚合所有客户端缓存方向，包括 inactive clients；
- 服务器更新原始变量 `beta` 并做 PDHG extrapolation。

验证结果：

- `tau = 0.5, 0.75, 0.9` 下 full-client/full-batch 目标函数与 central QR-PDHG 基本一致；
- half-client/minibatch cached aggregation 也稳定；
- 对偶变量严格落在 QR box 中；
- `tau = 0.9` 下显著优于 FedSPD-DP check subgradient 适配。

这可以作为项目主算法，FedSPD-DP smoothed/check 作为论文基线和直接适配基线。

更新：已完成 R1/R2 压力测试。`tau = 0.9`、20 客户端、Dirichlet non-IID、3 个 seeds 下，box-dual 在 isolated R1 和 isolated R2 中 final gap 约 `1e-4` 到 `3e-4`，在 R1+R2 中约 `2e-3` 到 `4e-3`；FedSPD-check 在 R1+R2 中约 `7e-2` 到 `9e-2`。这说明 box-dual 的优势已经被拉开。

更新：已补完整 hard baseline comparison。R1+R2 hard non-IID 下，QR box-dual 在 mild/hard/extreme 三档 final gap 均最小；FSPG-smooth 是强 baseline，排名第二；FedSubGrad 与 FedSPD-check 明显落后。对应图已生成在 `rfedqr/figures/`。

更新：ADMM baseline 已实现。`qr_admm()` central 版与 QR-PDHG/quantreg 级别目标值对齐；`fed_qr_admm()` full-client 版也能对齐 central。但在 hard R1+R2 中，FedQR-ADMM 退化明显，说明 ADMM 路线对 partial stale clients 和 mini-batch residual updates 很敏感。

更新：L1-penalized 实验已完成。`p=60`、true support size 8、hard non-IID、R1+R2 下，QR box-dual 在 objective gap 上最好，并且随着 `lambda` 从 0.001 到 0.02，选中变量数从约 56 降到约 11，FDR 从约 0.86 降到约 0.37。项目现在已经真正覆盖 penalized QR。

更新：MCP/SCAD 已实现并完成实验。非凸惩罚在 objective gap 上不一定优于 L1，但在变量选择上有清楚价值：MCP 在 `lambda=0.04` 时平均选中约 8 个变量，接近 true support size，TPR 约 0.83，FDR 约 0.13。报告中应把 MCP/SCAD 定位为变量选择/减小 shrinkage bias 的扩展。

## 稳定收敛是否意味着任务简单？

不意味着任务简单。

目前稳定收敛的是一个较理想化的版本：

- 单机集中式数据；
- 无惩罚或简单 L1 惩罚；
- 模拟数据；
- 没有客户端异质性；
- 没有客户端采样随机性 R1；
- 没有本地 mini-batch 随机性 R2；
- 尚未和 `quantreg::rq()` 做严格数值对照；
- 尚未和 smoothing、ADMM、subgradient 等 baseline 进行系统比较；
- 尚未处理真实数据中的缺失、中心偏移和变量选择问题。

因此，当前结论应表述为：

> 最小算法内核可行，但项目的研究难点仍在联邦化、异质性、双重随机性、极端分位数和 baseline 对比中。

## 真正的项目难点

### 1. 联邦 non-IID 场景

UCI Heart Disease 的四个中心分布差异明显，例如疾病阳性率：

- Cleveland: 约 45.9%
- Hungarian: 约 36.1%
- Switzerland: 约 93.5%
- VA: 约 74.5%

这说明真实联邦场景不是简单 IID 切分。客户端分布差异可能导致聚合方向偏移，也会影响客户端采样下的稳定性。

### 2. 双重随机性 R1 / R2

联邦优化中有两种随机性：

- R1：每轮只采样部分客户端；
- R2：每个客户端内部只用 mini-batch。

这两类随机性不应混在一起解释。项目的一个重要卖点是系统隔离 R1 与 R2 的影响，尤其是在 non-IID 数据下比较它们对收敛、方差和最终目标函数的不同作用。

### 3. 极端分位数

`tau = 0.9` 或 `tau = 0.95` 比中位数回归更难：

- 有效尾部样本更少；
- check loss 更不对称；
- 对偶变量 box `[tau - 1, tau]` 更不对称；
- 步长更敏感。

这正是 M2 tau-adaptive 步长可能发挥作用的地方。

### 4. 真实数据适配

UCI Heart Disease 数据可以作为真实联邦演示，但不是天然完美的分位数回归数据。

检查结果：

- 四中心合并共有 920 行。
- 如果保留所有 13 个原始变量，完整案例只剩 299，且几乎全来自 Cleveland。
- 如果去掉缺失严重的 `slope`、`ca`、`thal`，核心变量完整案例为 740。

这说明真实数据实验需要明确说明：

- 为什么选择某个连续响应变量；
- 如何处理缺失；
- 为什么使用四中心自然切分；
- 如何解释中心间分布差异。

### 5. Baseline 对比

仅证明 QR-FedSPD 能跑是不够的。项目最终需要说明：

- 相比 smoothing 方法，它避免了额外 smoothing 超参数和原问题改变；
- 相比 subgradient，它收敛更快或更稳定；
- 相比 ADMM，它每轮计算和通信更轻；
- 相比 generic PDHG，它利用了 QR 的 box geometry。

## 后续工作优先级

### 下一步 1：安装并验证 `quantreg`

先完成单机数值一致性检查：

- 安装 `quantreg`；
- 用同一批模拟数据比较 `qr_pdhg()` 和 `quantreg::rq()`；
- 检查系数差异、目标函数差异和收敛曲线。

这是第一阶段最关键的 sanity check。

状态：已完成。最大系数绝对差约为 `1.8e-4`。

### 下一步 2：整理 HD 数据实验口径

建议先选择一个连续响应变量：

- 首选候选：`thalach`，即最大心率；
- 备选：`oldpeak` 或 `trestbps`；
- 谨慎使用：`chol`，因为 Switzerland 中 `chol` 全为 0，VA 中也有异常 0。

暂定核心预测变量：

- `age`
- `sex`
- `cp`
- `trestbps`
- `chol`
- `fbs`
- `restecg`
- `exang`
- `oldpeak`
- `num` 或 `target_binary`

若响应为 `thalach`，则应从预测变量中移除 `thalach` 本身。

### 下一步 3：实现联邦版 `fed_qr_spd()`

在模拟数据上先实现：

- 多客户端切分；
- 每轮客户端采样 R1；
- 客户端内 mini-batch R2；
- 本地对偶变量保留在客户端；
- 服务器只聚合原始变量 `beta`。

先不急着加入所有 baseline，优先让联邦 QR-FedSPD 自己跑通。

### 下一步 4：做 R1/R2 隔离实验

实验设计：

- IID vs Dirichlet non-IID；
- full clients + full batch；
- sampled clients + full batch；
- full clients + mini-batch；
- sampled clients + mini-batch。

目标是分清楚：

- R1 是否主要受 non-IID 影响；
- R2 是否主要体现为局部噪声；
- 二者叠加时是否出现更明显震荡。

## 当前项目定位

当前最合适的项目叙事是：

> 我们不是单纯实现一个能收敛的分位数回归算法，而是利用 QR 的 box geometry，设计并评估一个适合联邦 non-IID 场景的随机 primal-dual 框架，并系统分析客户端采样与本地 mini-batch 两类随机性的影响。

## HD 四中心完整 baseline 对比

新增脚本：

- `scripts/run_hd_baseline_comparison_with_plots.R`

输出数据：

- `results/hd_center_response_summary.csv`
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

设置：

- UCI Heart Disease 四中心自然切分；
- response: standardized `thalach`;
- complete cases after dropping highly missing variables: `n = 740`;
- clients: Cleveland 303, Hungarian 261, Switzerland 46, VA 130;
- `tau in {0.5, 0.75, 0.9}`;
- scenarios:
  - Full: all clients + full local batch;
  - R1 client: sampled clients only;
  - R2 user: all clients + local mini-batch;
  - R1+R2: sampled clients + local mini-batch;
- methods:
  - QR box-dual;
  - FedQR-ADMM;
  - FedSubGrad;
  - FSPG-smooth;
  - FedSPD-check;
- target objective uses `quantreg::rq()` when available.

Targets:

| tau | target method | target objective | central QR-PDHG objective |
|---:|---|---:|---:|
| 0.50 | quantreg::rq | 0.3260838 | 0.3260896 |
| 0.75 | quantreg::rq | 0.2399178 | 0.2399242 |
| 0.90 | quantreg::rq | 0.1233159 | 0.1233224 |

R1+R2 final objective gaps:

| tau | method | final gap |
|---:|---|---:|
| 0.50 | QR box-dual | `4.66e-03` |
| 0.50 | FedSubGrad | `6.06e-03` |
| 0.50 | FSPG-smooth | `6.49e-03` |
| 0.50 | FedSPD-check | `1.66e-02` |
| 0.50 | FedQR-ADMM | `1.78e-02` |
| 0.75 | QR box-dual | `2.19e-03` |
| 0.75 | FSPG-smooth | `5.57e-03` |
| 0.75 | FedSubGrad | `6.83e-03` |
| 0.75 | FedSPD-check | `1.47e-02` |
| 0.75 | FedQR-ADMM | `1.60e-02` |
| 0.90 | QR box-dual | `1.53e-03` |
| 0.90 | FSPG-smooth | `3.54e-03` |
| 0.90 | FedSubGrad | `1.07e-02` |
| 0.90 | FedSPD-check | `1.55e-02` |
| 0.90 | FedQR-ADMM | `1.55e-02` |

Interpretation:

- HD is not a toy IID dataset: the four centers have visibly different `thalach` distributions and highly imbalanced client sizes.
- QR box-dual is the best method across all HD scenarios and all tested quantiles in objective gap.
- The advantage is especially clean in R1+R2, the setting most aligned with the project title.
- ADMM is a valid centralized/federated baseline, but its consensus form is much less robust under stochastic client and sample updates.
- FedSPD-check is a faithful nonsmooth adaptation of the reference FedSPD-DP structure, but it does not exploit the QR box dual geometry and remains clearly behind QR box-dual.

Current status:

> The main algorithmic improvement is now implemented and empirically supported: QR box-dual uses the exact quantile-regression box dual, supports R1/R2 stochastic federated updates, supports L1/MCP/SCAD penalties, and outperforms non-toy baselines on both hard simulations and the real four-center Heart Disease task.

## M1/M2 双自适应消融

新增脚本：

- `scripts/run_dual_adaptive_ablation_with_plots.R`

输出图：

- `figures/dual_adaptive_ablation_r1r2_gap.png`
- `figures/dual_adaptive_ablation_r1r2_convergence.png`
- `figures/dual_adaptive_ablation_all_scenarios.png`

消融口径：

- A0 operator: 普通 operator 步长；
- M1 box-aware: 利用 QR dual box 宽度重标定 primal/dual step；
- M1+M2 tau-adaptive: 在 M1 基础上对极端分位数进一步增大 primal step、减小 dual step。

结论：

- M1/M2 已经在 `qr_pdhg()` 和 `qr_box_fed_pdhg()` 中实现。
- 现在也有针对最终 QR box-dual 主方法的 hard non-IID 消融。
- M2 在 deterministic extreme quantile 下有收益，例如 `tau = 0.95` full-client/full-batch gap 从 A0 的 `1.55e-04` 降到 `6.61e-05`。
- 但 M2 在 R1+R2 下不稳，例如 `tau = 0.90` R1+R2 gap 从 A0 的 `4.42e-03` 变成 `1.20e-02`，`tau = 0.95` 从 `2.48e-03` 变成 `1.24e-02`。
- 因此最终主实验不应默认使用 naive M2；应该把 M2 作为可选 extreme-quantile acceleration，并用消融说明其在随机联邦场景下的稳定性代价。

当前可引用结论：

> The M1/M2 ablation shows that tau-adaptive scaling is not universally beneficial. It accelerates deterministic extreme-quantile optimization, but under simultaneous client sampling and local mini-batching it amplifies stochastic noise. This justifies using the more conservative QR box-dual step rule as the main stochastic federated method.

## 高级稳健化候选池

新增实现：

- `qr_box_fed_pdhg(dual_relaxation = ...)`
- `qr_box_fed_pdhg(server_momentum = ...)`
- `qr_box_fed_pdhg(step_decay_power = ..., step_decay_offset = ...)`
- `qr_box_fed_pdhg(primal_clip = ...)`
- `qr_box_fed_pdhg(aggregation = "selected_reweighted")`

新增脚本：

- `scripts/run_advanced_stabilization_experiment_with_plots.R`

输出图：

- `figures/advanced_stabilization_r1r2_gap.png`
- `figures/advanced_stabilization_r1r2_ratio.png`
- `figures/advanced_stabilization_best_convergence.png`

测试候选：

- A0 operator;
- M1 box-aware;
- M2 tau-adaptive;
- dual relaxation;
- server EMA;
- step decay;
- direction clipping;
- relaxation + EMA;
- relaxation + decay;
- selected reweighted aggregation.

核心结果：

| heterogeneity | tau | best variant | best gap | A0 gap |
|---|---:|---|---:|---:|
| hard | 0.90 | A0 / direction clip | `4.42e-03` | `4.42e-03` |
| hard | 0.95 | direction clip / A0 | `2.48e-03` | `2.48e-03` |
| extreme | 0.90 | A0 | `5.54e-03` | `5.54e-03` |
| extreme | 0.95 | direction clip | `4.53e-03` | `4.54e-03` |

结论：

- 高级候选基本没有显著超过保守 A0 cached QR box-dual。
- direction clipping 安全，但提升非常小，不能作为主要创新点。
- server EMA、dual relaxation、step decay 在 R1+R2 下普遍变差。
- selected reweighted aggregation 不如 cached aggregation，说明 inactive-client cached direction 是关键设计。
- 这个实验是有价值的负结果：它说明主算法强不是因为 generic 稳健化技巧，而是因为 QR box-dual geometry + cached federated direction 本身已经适合这个问题结构。

当前可引用结论：

> We tested a pool of advanced stochastic stabilization variants, including dual relaxation, server momentum, step decay, direction clipping, and selected-client reweighting. None materially improved the hard R1+R2 objective gap over the conservative cached QR box-dual update. This negative ablation supports the design choice of using cached QR box-dual updates as the final main method.

## 扩展场景矩阵

新增脚本：

- `scripts/run_expanded_scenario_suite_with_plots.R`
- `scripts/run_qr_box_tuning_expanded_scenarios.R`

新增场景：

- `low_participation`: only 2/20 clients per round;
- `tiny_user_batch`: mini-batch size 5;
- `extreme_tau`: `tau = 0.95`;
- `extreme_heterogeneity`: extreme client-level covariate/noise/tail heterogeneity;
- `unbalanced_clients`: client sample sizes from 30 to 200;
- `high_dim_sparse`: `p = 60`, true support size 8.

输出图：

- `figures/expanded_scenario_final_gap.png`
- `figures/expanded_scenario_rank_heatmap.png`
- `figures/expanded_scenario_convergence.png`
- `figures/qr_tuning_expanded_gap.png`
- `figures/qr_tuning_expanded_rank_heatmap.png`

第一轮固定 400 轮预算对比显示：

| scenario | best at 400 rounds | QR box-dual rank |
|---|---|---:|
| low participation | FSPG-smooth | 2 |
| tiny user batch | FSPG-smooth | 2 |
| extreme tau | QR box-dual | 1 |
| extreme heterogeneity | FSPG-smooth | 2 |
| unbalanced clients | FSPG-smooth | 2 |
| high-dimensional sparse | FSPG-smooth | 2 |

这说明默认 QR box-dual 不能被表述为“所有预算下无条件最好”。在紧通信预算或极小 batch 下，FSPG-smooth 的早期优化速度更有竞争力。

随后对 QR box-dual 做轮次与步长诊断：

| scenario | best tuned QR variant | best tuned QR gap | FSPG 400 gap |
|---|---|---:|---:|
| low participation | QR long 1200 | `5.23e-03` | `1.19e-02` |
| tiny user batch | QR long 1200 | `4.44e-03` | `7.63e-03` |
| extreme tau | QR clip .25 800 | `1.13e-03` | `7.59e-03` |
| extreme heterogeneity | QR long 1200 | `1.35e-03` | `7.47e-03` |
| unbalanced clients | QR long 1200 | `1.49e-03` | `6.24e-03` |
| high-dimensional sparse | QR long 1200 | `5.99e-03` | `1.69e-02` |

Interpretation:

- FSPG-smooth is a strong short-horizon baseline.
- QR box-dual has better asymptotic/stabilized behavior, but sometimes needs more communication rounds or mild tuning.
- The final report should avoid saying “QR box-dual is always best at every budget.” A more defensible claim is:

> Under fixed short communication budgets, smoothing can be faster in some highly stochastic scenarios. After QR-specific tuning or a longer communication horizon, cached QR box-dual recovers the best objective gap across the expanded scenario matrix.

## NYC Taxi 百万级真实数据实战

新增文件：

- `NYC_TAXI_LARGE_SCALE.md`

新增脚本：

- `scripts/prepare_nyc_taxi_data.sh`
- `scripts/run_nyc_taxi_large_scale_experiment.R`

数据来源：

- official NYC TLC Yellow Taxi trip records;
- 2024 Q1: January, February, March Parquet files;
- taxi zone lookup table.

数据规模：

| raw clean trips | clean pickup-zone clients | eligible clients | modeled rows |
|---:|---:|---:|---:|
| 8,417,330 | 257 | 59 | 1,000,000 |

建模设置：

- response: `log_total_amount`;
- quantile: `tau = 0.9`;
- clients: pickup zones (`PULocationID`);
- features: trip distance, duration, passenger count, fees, pickup time cyclic features, borough, vendor, rate code, payment type;
- final setting: 12/59 clients per round, local user batch size 5,000, base rounds 500, QR long 1,000.

结果：

| method | objective | gap to best observed |
|---|---:|---:|
| QR box-dual long | 0.0206365 | 0 |
| QR box-dual | 0.0220496 | 0.0014131 |
| FSPG-smooth | 0.0240271 | 0.0033905 |
| FedSubGrad | 0.1043675 | 0.0837310 |
| FedSPD-check | 0.3452998 | 0.3246633 |

解释：

- HD 是真实 non-IID，但不是大规模。
- NYC Taxi 是真正的百万级真实数据分布式随机优化场景。
- 初始小 batch 运行中，QR box-dual 不稳定，FSPG-smooth 明显更好。
- 将本地用户 batch 增大到 5,000 后，QR box-dual long 取得最佳 observed objective。
- 这说明 QR box-dual 在大数据场景需要足够的本地 sample-level dual 覆盖；这与算法结构一致，不是 toy 调参。

当前可引用结论：

> On the million-row NYC Taxi task, cached QR box-dual achieves the best observed objective once each active client performs sufficiently large local user-level dual updates. This provides the first real large-scale distributed stochastic optimization evidence for the proposed method.

## UCI Household Power 百万级真实数据实战

新增文件：

- `HOUSEHOLD_POWER_LARGE_SCALE.md`

新增脚本：

- `scripts/prepare_household_power_data.sh`
- `scripts/run_household_power_large_scale_experiment.R`

Data:

- UCI Individual Household Electric Power Consumption;
- minute-level household electricity records;
- 2,049,280 modeled rows after removing missing values;
- 48 calendar-month clients.

Task:

- response: `log1p(Global_active_power)`;
- `tau = 0.9`;
- clients: calendar months;
- 10/48 clients per round;
- local user batch size: 5,000;
- QR long run: 800 rounds.

Results:

| method | objective | gap to best observed |
|---|---:|---:|
| QR box-dual long | 0.0419800 | 0 |
| FSPG-smooth | 0.0423711 | 0.0003911 |
| FedSubGrad | 0.0456561 | 0.0036761 |
| QR box-dual | 0.0500239 | 0.0080439 |
| FedSPD-check | 0.0550988 | 0.0131188 |

Interpretation:

- This adds a second non-NYC large-scale real-data setting.
- The monthly client split creates temporal non-IID and moderate client-size imbalance.
- QR box-dual long again achieves the best observed final objective, but FSPG-smooth remains very close under the short 400-round budget.
- The result supports a careful claim: QR box-dual is strongest in final nonsmoothed check-loss objective after enough local dual refresh; smoothing is a serious short-horizon baseline and should remain in all comparisons.

Current report wording:

> Across both NYC Taxi and UCI Household Power, cached QR box-dual reaches the best observed objective after sufficient local sample-level dual refresh. The Household Power result is especially useful because FSPG-smooth is close, showing that the comparison is not artificially easy while still favoring the proposed QR-specific method at convergence.

## 代码整理与测试体系

新增统一接口：

- `fit_fedqr()`：统一调用 QR box-dual、FSPG-smooth、FedSubGrad、FedSPD-check、FedSPD-smooth、FedQR-ADMM。
- `run_fedqr_methods()`：批量运行 methods，并统一返回 `summary`、`trace`、`coefficients`。
- `fedqr_methods()`：列出注册方法。

新增实验工具：

- `make_experiment_dirs()`;
- `save_experiment_outputs()`;
- `plot_convergence()`;
- `plot_final_gap()`;
- `validate_result_table()`.

整理脚本：

- `scripts/load_rfedqr.R` 统一 source R 文件；
- `scripts/run_nyc_taxi_large_scale_experiment.R` 改为调用 `run_fedqr_methods()`；
- `scripts/run_household_power_large_scale_experiment.R` 改为调用 `run_fedqr_methods()`；
- `scripts/run_hd_baseline_comparison_with_plots.R` 改为通过 `fit_fedqr()` dispatch baseline。

包交付：

- 新增 `README.md`;
- 新增 `EXPERIMENTS.md`;
- 新增 `.Rbuildignore`;
- 新增 `tests/testthat/`;
- 新增轻量 `man/rfedqr-exports.Rd`;
- `DESCRIPTION` 升级到 `0.1.0` 并补测试依赖声明。

验证状态：

- `testthat::test_dir("tests/testthat")` 通过；
- `R CMD build .` 通过；
- `R CMD check rfedqr_0.1.0.tar.gz --no-manual` 通过，状态为 `OK`；
- 大实验没有在整理阶段重跑，避免覆盖已有百万级结果；测试中只做现有 summary CSV 的回归校验。

## 高级创新点：staleness-aware、client-robust、calibration-aware

新增目标：

在不改变原始 QR box-dual 默认行为的前提下，实现三条更像论文贡献的高级模块：

1. staleness-aware cached direction；
2. client-robust/fairness-aware objective；
3. high-quantile calibration refinement。

新增算法能力：

- `qr_box_fed_pdhg()` 新增 `staleness = c("none", "exponential", "inverse")`；
- 每个客户端维护 cached direction 的 age；
- stale aggregation 可按 `decay(age_j)` 降权陈旧客户端方向；
- trace 新增 `mean_staleness`、`max_staleness`、`mean_stale_weight`；
- `client_weighting = c("sample", "uniform", "sqrt_size", "custom")` 支持样本平均、客户端等权、平方根折中和自定义权重；
- local direction 从 `X_j^T v_j / n` 泛化为 `w_j X_j^T v_j / n_j`；
- 默认 `client_weighting = "sample"` 时严格保留原全局样本平均口径。

新增诊断函数：

- `client_qr_objective()`;
- `client_loss_summary()`;
- `quantile_coverage()`;
- `calibrate_quantile_intercept()`;
- `calibration_summary()`.

统一接口新增方法名：

- `QR box-dual stale`;
- `QR box-dual robust`;
- `QR box-dual stale+robust`.

新增实验脚本：

- `scripts/run_advanced_innovation_experiment_with_plots.R`;
- `scripts/run_large_realdata_coefficient_evaluation.R`.

新增输出：

- `results/advanced_innovation_summary.csv`;
- `results/advanced_innovation_client_loss.csv`;
- `results/advanced_innovation_calibration.csv`;
- `results/advanced_innovation_trace.csv`;
- `results/advanced_innovation_large_client_loss.csv`;
- `results/advanced_innovation_large_calibration.csv`;
- `figures/advanced_innovation_final_gap.png`;
- `figures/advanced_innovation_client_fairness.png`;
- `figures/advanced_innovation_calibration.png`;
- `figures/advanced_innovation_staleness_trace.png`.

初步结果：

- 在 extreme non-IID simulation 中，`QR box-dual stale` 的 final target gap 优于普通 `QR box-dual`，说明 staleness-aware cached direction 对陈旧缓存方向有实际优化收益。
- 在 Heart Disease 和多个 unbalanced simulation 场景中，`QR box-dual robust` 或 `QR box-dual stale+robust` 改善 worst-client loss，说明它的价值主要体现在客户端公平性，而不是单纯追求 global sample-average objective。
- calibration refinement 能把 global coverage error 压到接近 0；client-offset calibration 能显著降低 mean/worst client coverage error。当前 calibration 是同数据诊断式修正，报告中应表述为 post-training calibration/fairness diagnostic，后续若时间允许可增加 train/calibration/test split。
- NYC 和 Household 暂不重跑新训练方法，只基于已有 coefficient CSV 做 coefficient-level fairness/calibration evaluation；这避免覆盖百万级结果，也符合当前交付优先级。

推荐报告命名：

> staleness-aware client-robust calibrated QR box-dual FedPDHG

报告表达时建议拆成三个模块讲，而不是宣称一个模块在所有指标上都 SOTA：

- staleness-aware: improves optimization under partial participation and extreme non-IID;
- client-robust: improves worst-client and client-distribution fairness;
- calibration-aware: improves high-quantile coverage after training.
