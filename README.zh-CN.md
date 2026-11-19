# TabulaStromal

本次更新为 **2026-11-19** 经审计结果快照，包含 [F1 图谱](F1_StromalAtlas/)、[F2 模型表现](F2_ModelPerformance/)、[F3 下游机制关联](F3_ds_analysis/)、[Rebenchmark](TabulaStromal_Rebenchmark_2026/) 与 [CellAge](CellAge_Exploration/) 五个模块。

[英文首页及模型图](README.md) · [40 页下游 PPT](F3_ds_analysis/TabulaStromal_F3_Downstream_Results_20261119.pptx) · [三段英文 Results 与 Word/PDF](docs/results/) · [审计清单](audit/README.md)

![当前 Figure 2](F2_ModelPerformance/figures/Fig2_StromalClock_ModelPerformance.png)

当前 compact EN 保留 750/350 个 Celltype_Gene 特征、678/337 个基因、各 27 个状态。开发集 R²≈0.950 是训练内拟合；ALL100 外部 R²=0.377/0.391 的模型选择使用了外部结局，不能视为独立验证；REUSED70 也参与过选择。Figure 2 B–D 是旧 fold-4 benchmark，E/F 是当前开发集拟合。MAD 的两种实际定义已分别说明。

F3 包含稳定入选特征、功能富集、LISA2、RNA prior network、成人 raw-count ECM 和当前选定状态的 CellAge 结果。图和汇总表均已加入网页。共有基因、富集和调控关联不能直接写成共同因果衰老机制；相对 hUSI 状态差异不是随年龄增加或 age acceleration。

明确排除 `VenEC_CLU_eye_audit` 与 `F3_features/Age10_20_audit_20261117` 审计子树及全部后代。合法结果中的 Celltype 名称不等于被排除的审计目录。

完整元数据发现 1,469,578 条（含环境/历史文件），实际读取并复核哈希的核心文件 432 个。公开 229 个原始图像/PDF、201 个汇总/证据表和经路径泛化的原始 PPTX。完整 cache 路径清单只留在本地 ignored outputs；不上传表达矩阵、模型或单样本预测/评分行。

运行方法、依赖、可复现性边界和授权信息见[英文首页](README.md)、[数据规范](docs/data.md)、[复现说明](docs/reproducing.md)及 [NOTICE](NOTICE)。本仓库结果草稿不表示论文已经发表。
