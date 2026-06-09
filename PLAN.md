# rfedqr 项目实施计划

## 目标

本项目要实现一个面向联邦数据场景的惩罚分位数回归工具包 `rfedqr`。第一阶段先把单机 QR-PDHG 算法跑通，并用模拟数据验证它和经典集中式分位数回归的一致性；第二阶段再扩展到联邦随机 primal-dual 算法，系统研究客户端采样与本地 mini-batch 两类随机性。

## 阶段 0：项目骨架与最小闭环

交付物：

- `R/qr_pdhg.R`：单机 QR-PDHG。
- `R/prox.R`：惩罚项 proximal operator，第一版支持 `none` 和 `l1`。
- `R/objective.R`：check loss 与目标函数。
- `R/simulate.R`：模拟数据生成与客户端切分工具。
- `scripts/run_simulation.R`：第一批模拟实验。
- `scripts/inspect_heart_disease.R`：下载并概览 UCI Heart Disease 数据。

优先验证：

1. 无惩罚 QR-PDHG 是否能收敛。
2. `tau = 0.5, 0.75, 0.9` 下目标函数是否下降。
3. 如果本地安装了 `quantreg`，系数是否与 `quantreg::rq()` 接近。
4. UCI Heart Disease 四个中心的数据规模、缺失、响应变量分布是否适合做真实联邦演示。

## 阶段 1：M1 / M2 步长消融

M1：Box-Aware 步长。

- 利用 QR 对偶变量 `v_i in [tau - 1, tau]` 的 box 约束。
- 与通用 PDHG 步长比较收敛速度和稳定性。

M2：tau-Adaptive 步长。

- 极端分位数下减小对偶步长、增大原始步长。
- 比较 `tau = 0.1, 0.5, 0.9`。

指标：

- 目标函数值 vs iteration。
- 目标函数值 vs wall-clock。
- 系数误差。
- 调参敏感性。

## 阶段 2：联邦模拟

实现：

- 按 IID 与 Dirichlet non-IID 生成客户端。
- R1：每轮客户端采样。
- R2：客户端内 mini-batch。
- 只通信原始变量 `beta`，对偶变量 `v_i` 保留在本地。

核心问题：

- non-IID 下 R1 是否显著增加方差。
- R2 是否主要影响局部噪声但对异质性不敏感。
- M1/M2 是否在联邦环境仍然稳定。

## 阶段 3：baseline 对比

优先级：

1. FedSubGrad：实现简单，作为慢但稳的下界。
2. FSPG-style smoothing：将 check loss 平滑化后做随机梯度。
3. FedSPD-DP-style generic PDHG：不使用 QR box 专门步长。
4. FS-QRADMM-style baseline：如果时间允许再做，因实现和调参成本较高。

比较维度：

- 通信轮次。
- 每轮计算成本。
- 总 wall-clock。
- 最终目标函数值。
- 系数估计误差。

## 阶段 4：真实数据

### UCI Heart Disease

用途：小规模真实联邦演示。

- 四个自然中心：Cleveland、Hungarian、Switzerland、VA。
- 样本量小，适合端到端 sanity check。
- 重点检查缺失值和中心间分布差异。

### NYC Taxi

用途：scalability 与异质性测试。

- 按 borough 或 pickup zone 切分客户端。
- 数据量从 1M 到 10M 逐步扩大。
- 放在后期，避免前期被数据工程拖慢。

## 当前建议顺序

1. 运行 `scripts/run_simulation.R`，确认单机 QR-PDHG 能收敛。
2. 运行 `scripts/inspect_heart_disease.R`，了解 HD 数据结构和缺失情况。
3. 根据模拟结果微调默认步长。
4. 加入联邦版 `fed_qr_spd()`。
5. 做 R1/R2 隔离实验。
6. 再处理真实数据实验与报告图表。

