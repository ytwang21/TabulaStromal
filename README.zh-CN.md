# TabulaStromal

本仓库整理人类成纤维细胞与内皮细胞图谱、年龄预测、稳定特征及 CellAge 关联分析的项目代码。英文入口和命令见 [README.md](README.md)。

包含的主要工作流：

- `workflows/atlas/`：图谱审计、标志物、轨迹及 Figure 1 相关代码。
- `workflows/benchmark/`：13 模型基准测试、分组交叉验证、外部评估、稳定性及最终 ElasticNet 分支。
- `workflows/cellage/`：已有 hUSI、SenCID recommendation 输出的审计和关联分析。
- `workflows/model_performance/`：历史 Figure 2 构图代码。
- `legacy/`：历史注释、整合及探索笔记本；已移除运行输出。
- `results/`：真实结果的聚合指标和稳定特征表，不含个人层面的预测和元数据。

建议使用 Python 3.11，从仓库根目录执行：

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements/benchmark.txt
python scripts/validate_repository.py
python scripts/run_workflow.py benchmark --dry-run
```

实际数据未随代码发布，默认配置中的输入不存在属于预期情况。输入格式、项目原目录映射及在原服务器生成本地配置的方法见 [数据说明](docs/data.md) 和 [复现说明](docs/reproducing.md)。最终 EN 管线需要原始验证过的 checkpoint 及历史图形来源文件，不是无数据即可一键运行的示例；独立推理接口可使用自行训练或合法取得的可信 bundle。

请保留科学范围的区别：样本标签不等于已经独立核实的供体；最终 EN 不等于全部模型中的内部最优模型；FT30 ALL100 包含适配样本，不能替代独立留出结果；hUSI 不等于以年为单位的生物年龄，recSID 不等于已经核实的衰老强度。负的外部 R² 在仓库中如实报告。

新增四个结果目录的图件展示、文件与 Figure 清单、证据登记表，以及 [三个英文 Results 与可编辑 Word](docs/results/)。完整逐文件元数据清单保存在本地忽略目录，公开审计摘要覆盖全部扫描记录。

原始矩阵、个人层面的预测/元数据、模型和运行缓存未复制到公开仓库。项目尚未选择整体开源许可证，详见 [NOTICE](NOTICE)。
