#!/usr/bin/env python
"""Produce a concise phase/output/scientific findings report from actual saved artifacts."""
import argparse
import json
from pathlib import Path
import pandas as pd


def report(root):
    root=Path(root)
    manifest=json.loads((root/'run_manifest.json').read_text())
    lines=['# 实际运行交付记录',f'\n运行目录：`{root.resolve()}`', '\n## 阶段状态\n']
    for key,status in manifest['stages'].items():
        stage,lin=key.split(':')
        out=root/lin if lin!='shared' else root/'integrated_figure'
        lines.append(f'- **{key}**：{status["status"]}；耗时 {status.get("seconds","—")} 秒；输出 `{out}`；问题：{status.get("error","无执行错误")}。')
    for lin in ['fibro','endo']:
        p=root/lin;a=root/'shared'/f'data_audit_{lin}.json'
        if not a.exists():continue
        d=json.loads(a.read_text())
        lines.extend([f'\n## {lin}\n',f'- 审计：{d["n_cells"]:,} cells，{d["n_genes"]:,} genes，{d["n_tissues"]} tissues，{d["n_samples"]} samples；donor 未提供。',
                      f'- 三层注释：general {len(d["categories"]["general"])} / subtype {len(d["categories"]["subtype"])} / Celltype {len(d["categories"]["celltype"])}。',
                      f'- Age_group Unknown：{d["missing"]["age_group"]/d["n_cells"]:.2%}。'])
        tissue=p/'summary_tissue.csv'
        if tissue.exists():
            t=pd.read_csv(tissue,index_col=0)
            lines.append(f'- 覆盖：细胞贡献最多组织 {t.n_cells.idxmax()}（{t.n_cells.max():,} cells）；样本最多组织 {t.n_samples.idxmax()}（{t.n_samples.max():,} samples）。')
        age=p/'summary_age_group.csv'
        if age.exists():
            t=pd.read_csv(age,index_col=0).drop(index='Unknown',errors='ignore')
            if len(t):lines.append(f'- 已知年龄分组内细胞最多：{t.n_cells.idxmax()}（{t.n_cells.max():,}）；样本最多：{t.n_samples.idxmax()}（{t.n_samples.max():,}）。A/H 含历史文本映射，不能全部视为精确岁数。')
        availability=p/'marker_availability.json'
        if availability.exists():
            data=json.loads(availability.read_text())
            lines.append(f'- 缺失 marker：{json.dumps(data["missing"],ensure_ascii=False)}。不补造缺失表达。')
        concord=p/'signature_concordance.json'
        if concord.exists():
            value=json.loads(concord.read_text())['top_score_agreement']
            lines.append(f'- 最高 signature 与注释一致：{value:.2%}（描述性，marker 与注释非独立，不是分类准确率）。')
        metrics=p/'trajectory'/'trajectory_method_metrics.csv'
        if metrics.exists():
            t=pd.read_csv(metrics);r=t.iloc[0]
            lines.append(f'- 轨迹 reference `{r.method}`：根分量有限 DPT {r.finite_dpt_fraction:.2%}；替代根 Spearman {r.root_sensitivity_spearman:.3f}；组织 eta² {r.pseudotime_tissue_eta2:.3f}。')
            lines.append(f'- 完整方法比较、稳定性判据和限制：`{p / "trajectory" / "trajectory_report.md"}`。不能据此证实分化或转分化。')
        comparison=p/'trajectory'/'three_method_comparison.csv'
        if comparison.exists():
            audit=json.loads((p/'trajectory/metacell_audit.json').read_text())
            selected=json.loads((p/'trajectory/trajectory_selection.json').read_text())
            lines.append(f'- 轨迹输入：{audit["n_meta_cells"]:,} 个原始 meta-cells；共同拟合集合 {audit["modeling_cohort"]["n_fit_metacells"]:,} 个；单个 meta-cell 的原始细胞数中位数 {audit["n_cells_in_meta_quantiles"]["median"]:g}。')
            lines.append(f'- 已实际比较 Monocle 2、STREAM、CellRank；主图展示 `{selected["model"]}`。不使用 PAGA；该展示排序不是正确分化方向的验证。')
            lines.append(f'- 比较表：`{comparison}`；原始轨迹坐标、特征集、转移矩阵/真实命运概率及解释：`{p / "trajectory" / "trajectory_report.md"}`。')
        lines.append(f'- 单图文件数：PNG {len(list(p.rglob("*.png")))} / PDF {len(list(p.rglob("*.pdf")))}。')
    remaining=[f'{s}:{l}' for s in ['audit','atlas','markers','trajectory'] for l in ['fibro','endo'] if manifest['stages'].get(f'{s}:{l}',{}).get('status')!='complete']
    if manifest['stages'].get('figure:shared',{}).get('status')!='complete':remaining.append('figure:shared')
    lines.extend(['\n## 剩余工作\n', '未完成阶段：'+(', '.join(remaining) if remaining else '无；本次全流程已完成。'),
                  '\n科研解释待办：主图使用前检查 trajectory_report.md 的不稳定/混杂提示，并核实缺失 marker、注释和年龄元数据。此为结果解释限制，不是未完成的计算阶段。'])
    path=root/'delivery_report.md';path.write_text('\n'.join(lines)+'\n')
    return path


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('run_dir',type=Path)
    print(report(parser.parse_args().run_dir))
