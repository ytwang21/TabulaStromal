#!/usr/bin/env python
"""Build the complete meta-cell Monocle/STREAM/CellRank analysis notebook."""
from pathlib import Path
import nbformat as nbf
HERE=Path(__file__).resolve().parent
md=nbf.v4.new_markdown_cell;code=nbf.v4.new_code_cell
nb=nbf.v4.new_notebook(metadata={"kernelspec":{"display_name":"Python 3 (StromalAtlas)","language":"python","name":"python3"},"language_info":{"name":"python","version":"3.11"}})
nb.cells=[
md("""# StromalAtlas：全细胞 atlas + meta-cell 三工具轨迹比较

UMAP、覆盖、marker/signature 使用原来的全细胞 atlas。轨迹使用用户指定的 Fibro 47,646 / Endo 19,011 meta-cell 对象，**不使用 PAGA**。三种方法均在相同、固定的 6,000 个代表性 meta-cells 上比较（完整 ID 记录），不将单细胞替代成 meta-cell；可在新运行配置 max_fit_metacells=0 使用全部 meta-cells。

实际比较：Monocle 2 DDRTree 的六个特征设置、STREAM 两种表示的真实弹性主树、CellRank 真实转移矩阵及命运概率。树形清晰、marker 趋势或概率流向均不等于发育方向的独立验证。"""),
code("""import os,sys,json,subprocess
from pathlib import Path
import pandas as pd
from IPython.display import display,Image,Markdown
candidates=[Path.cwd(),Path.cwd()/"StromalAtlas"/"code",Path.cwd()/"code"]
CODE=next(p.resolve() for p in candidates if (p/"stromal_atlas_pipeline.py").exists())
sys.path.insert(0,str(CODE))
import stromal_atlas_pipeline
RUN_DIR = CODE.parent / "results" / "runs" / "atlas_metacell_20261031"
CONFIG = CODE / "config.json"
print("Python:",sys.executable,"Run:",RUN_DIR)
display(pd.Series(json.loads(CONFIG.read_text()),name="value").to_frame())"""),
code("""def run_phase(stage):
    command=[sys.executable,"-u",str(CODE/"stromal_atlas_pipeline.py"),"--config",str(CONFIG),"--run-dir",str(RUN_DIR),"--stage",stage]
    if (RUN_DIR/"run_manifest.json").exists():command.append("--resume")
    subprocess.run(command,check=True,cwd=CODE)
def show(relative,width=1100):
    display(Image(filename=str(RUN_DIR/relative),width=width))
def table(relative):return pd.read_csv(RUN_DIR/relative)
"""),
md("## Phase 0：全细胞 atlas 审计。已有成功阶段 --resume 会跳过；本次复用前次已验证的 atlas 结果，轨迹重新执行。"),
code('run_phase("audit")'),
code("""for lin in ["fibro","endo"]:
    d=json.loads((RUN_DIR/"shared"/f"data_audit_{lin}.json").read_text())
    display(pd.Series({k:d[k] for k in ["n_cells","n_genes","n_tissues","n_samples","n_donors","fields","obsm","has_raw","layers"]},name=lin).to_frame())
    print("obs.columns:",d["obs_columns"])
    display(pd.Series(d["categories"]["general"],name="cells").to_frame())"""),
md("## Phase 1：UMAP 与覆盖。虚线为自动密度轮廓，分离岛保持同色；Unknown 年龄保留，样本不是供体。"),
code('run_phase("atlas")'),
code("""for lin in ["fibro","endo"]:
    show(f"{lin}/{lin}_umap_subtype_dashed.png")
    show(f"{lin}/{lin}_dataset_summary.png",850)
    display(table(f"{lin}/summary_tissue.csv"))
    display(table(f"{lin}/summary_age_group.csv"))"""),
md("## Phase 2：全细胞 marker/signature。marker 组内字母序，缺失基因记录；固定表达匹配控制的 module score 是描述性注释支持，不是独立验证。"),
code('run_phase("markers")'),
code("""for lin in ["fibro","endo"]:
    display(json.loads((RUN_DIR/lin/"marker_availability.json").read_text()))
    show(f"{lin}/{lin}_marker_dotplot.png")
    show(f"{lin}/{lin}_signature_heatmap.png",650)
    show(f"{lin}/{lin}_signature_sample_boxplots.png")"""),
md("""## Phase 3：meta-cell 轨迹（三工具）

Monocle：HVG750/1000/1500 及 DEG-informed50/150/350；平均 log 表达使用 uninormal，不伪装成 counts。DEG-informed 特征依赖注释，不能独立验证预期分支。

STREAM：复现七中心 seed 后，实际执行 ElPiGraph 弹性树优化，比较 PCA10 与历史 UMAP 表示；不会把 seed MST 当成最终优化结果。

CellRank：0.8 PseudotimeKernel + 0.2 ConnectivityKernel。200-landmark Markov lumping 上进行 GPCCA，随后在完整共同拟合集合的稀疏转移矩阵上求解命运概率。流向是 T@UMAP−UMAP，不是 DPT 梯度。DPT 根方向仍是先验，不是 RNA velocity。

三个工具使用相同的 meta-cell ID。少数 meta-cell 含数千至数万原始细胞，是历史 leftover 聚合产生的长尾，记录并提供诊断。"""),
code('run_phase("trajectory")'),
code("""for lin in ["fibro","endo"]:
    display(json.loads((RUN_DIR/lin/"trajectory/metacell_audit.json").read_text()))
    display(json.loads((RUN_DIR/lin/"trajectory/prepared/cohort.json").read_text()))
    display(table(f"{lin}/trajectory/three_method_comparison.csv"))
    display(table(f"{lin}/trajectory/pseudotime_concordance.csv"))
    show(f"{lin}/trajectory/{lin}_monocle_feature_comparison.png",1300)
    show(f"{lin}/trajectory/{lin}_three_method_comparison.png",1300)
    show(f"{lin}/trajectory/{lin}_marker_trends.png",1300)
    display(Markdown((RUN_DIR/lin/"trajectory/trajectory_report.md").read_text()))"""),
md("## Phase 4：重新排版。A–H 保留全细胞 atlas 结果，I–J 使用三工具比较后选出的 meta-cell 原生轨迹。另给出两谱系×三工具并排大图。不强行画预设的 subtype 间分化箭头。"),
code('run_phase("figure")'),
code("""show("integrated_figure/Fig_StromalAtlas_Overview.png",1500)
show("integrated_figure/Fig_StromalAtlas_Trajectory_Comparison.png",1500)
print((RUN_DIR/"integrated_figure/figure_caption.txt").read_text())
display(json.loads((RUN_DIR/"validation_report.json").read_text()))
display(pd.DataFrame(json.loads((RUN_DIR/"run_manifest.json").read_text())["stages"]).T)
"""),
md("""## 后台运行与分步骤

参见 README.md。`bash run_nohup.sh` 新建完整运行；`--reuse-atlas-run ../results/runs/atlas_20261031` 可复用 atlas 部分并重新比较 meta-cell 轨迹。逐工具使用 `metacell_trajectory.py --method prepare|monocle|stream|cellrank|compare --lineage fibro|endo --run-dir PATH`；所有调用使用同一目录。每个已成功的原生拟合均有单独检查点，失败后不会重跑全部方法。

当前不存在可证明“正确”发育树的 lineage tracing/真实时间/velocity 证据；展示推荐是有明确假设的描述性方法比较。""")]
nbf.write(nb,HERE/"stromal_atlas_pipeline.ipynb")
print(HERE/"stromal_atlas_pipeline.ipynb")
