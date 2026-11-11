#!/usr/bin/env python
"""Independently audit metrics and require identical ALL100 sample pairs.

No trained model or original benchmark file is modified. The renderer uses this
same guard and the independently recomputed metrics, rather than mixing tables.
"""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import argparse
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
METRICS = ('R2', 'MAE', 'MAD', 'Age_MAD', 'MedAE', 'RMSE', 'PCC')
KEYS = ['model', 'K', 'fold']


def independent_metrics(age, prediction):
    """Direct formulas; independent of tabulastromal.metrics/sklearn scoring."""
    y = np.asarray(age, dtype=float)
    p = np.asarray(prediction, dtype=float)
    if y.ndim != 1 or y.shape != p.shape or not len(y):
        raise ValueError('Expected nonempty paired one-dimensional ages/predictions')
    if not np.isfinite(y).all() or not np.isfinite(p).all():
        raise ValueError('Nonfinite values: sample filtering is forbidden')
    e = y - p
    yc, pc = y - y.mean(), p - p.mean()
    sse, sst, pss = float(e @ e), float(yc @ yc), float(pc @ pc)
    age_mad = float(np.abs(yc).mean())
    return dict(n=len(y), R2=1-sse/sst if len(y) > 1 and sst > 0 else np.nan,
                MAE=float(np.abs(e).mean()), MAD=age_mad, Age_MAD=age_mad,
                MedAE=float(np.median(np.abs(e))), RMSE=float(np.sqrt(sse/len(y))),
                PCC=float((yc @ pc)/np.sqrt(sst*pss)) if len(y) > 1 and sst > 0 and pss > 0 else np.nan,
                Age_mean=float(y.mean()), Age_SD=float(np.sqrt(sst/len(y))),
                Age_min=float(y.min()), Age_max=float(y.max()), SSE=sse, SST=sst,
                normalized_MSE=sse/sst if sst > 0 else np.nan,
                max_absolute_error=float(np.abs(e).max()))


def verified_external_metrics(metrics, predictions, split):
    """Fail on any scope/table mismatch; return freshly computed ALL100 scores."""
    if split.empty or split.sample_id.duplicated().any():
        raise ValueError('External split must identify each sample exactly once')
    ids = set(split.sample_id.astype(str))
    if metrics.empty or metrics.duplicated(KEYS).any():
        raise ValueError('Metric rows must be unique per model/K/fold')
    if not metrics.source.eq('external_all_post_adaptation').all():
        raise ValueError('Primary metric source must be external_all_post_adaptation')
    if not predictions.source.eq('external_all_post_adaptation').all():
        raise ValueError('Primary predictions must be external_all_post_adaptation')
    mkeys = set(map(tuple, metrics[KEYS].to_numpy()))
    pkeys = set(map(tuple, predictions[KEYS].drop_duplicates().to_numpy()))
    if mkeys != pkeys:
        raise ValueError('Prediction/metric model-K-fold keys differ')
    result = metrics.copy()
    index = {tuple(r): i for i, r in zip(result.index, result[KEYS].to_numpy())}
    evidence = []
    for key, part in predictions.groupby(KEYS, sort=True):
        sample_ids = part.sample_id.astype(str)
        if sample_ids.duplicated().any() or set(sample_ids) != ids:
            raise ValueError(f'{key}: metrics require ALL100 exact unique sample IDs')
        if 'Age' in split:
            expected = split.assign(sample_id=split.sample_id.astype(str)).set_index('sample_id')
            np.testing.assert_allclose(part.Age, expected.loc[sample_ids, 'Age'], rtol=0, atol=1e-10)
        rowindex = index[key]
        row = result.loc[rowindex]
        computed = independent_metrics(part.Age, part.prediction)
        discrepancy = 0.
        for name in METRICS:
            np.testing.assert_allclose(computed[name], row[name], rtol=1e-8, atol=1e-8,
                                       equal_nan=True, err_msg=f'{key}: {name} not ALL100')
            if np.isfinite(computed[name]) and np.isfinite(row[name]):
                discrepancy = max(discrepancy, abs(computed[name]-float(row[name])))
            result.loc[rowindex, name] = computed[name]
        for name in ['n', 'n_evaluation']:
            if name in row and int(row[name]) != len(ids):
                raise ValueError(f'{key}: {name} is not the ALL100 sample count')
        if 'baseline_prediction' in part:
            before = independent_metrics(part.Age, part.baseline_prediction)
            for name in METRICS:
                if 'baseline_'+name in row:
                    np.testing.assert_allclose(before[name], row['baseline_'+name], rtol=1e-8, atol=1e-8,
                                               equal_nan=True, err_msg=f'{key}: baseline {name} not ALL100')
                    result.loc[rowindex, 'baseline_'+name] = before[name]
                if 'delta_'+name in row:
                    result.loc[rowindex, 'delta_'+name] = computed[name]-before[name]
        digest = hashlib.sha256(json.dumps(sorted(ids), ensure_ascii=False).encode()).hexdigest()
        evidence.append(dict(zip(KEYS, key), source='external_all_post_adaptation', **computed,
                             n_for_R2=len(ids), n_for_MAE=len(ids), n_for_MAD=len(ids), n_for_PCC=len(ids),
                             same_pairs_for_all_metrics=True, sample_ids_sha256=digest,
                             max_saved_metric_difference=discrepancy))
    return result, pd.DataFrame(evidence)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'Stromal_Combined_Figure/outputs')
    args = parser.parse_args()
    checks, diagnostics, summary = [], [], {}
    for label, root in [('Fibroblast', ROOT/'Fibro_Accelerated_20261102/results/Fibro'), ('Endothelial', ROOT/'Endo')]:
        split = pd.read_csv(root/'external/external_ft30_split.csv', low_memory=False)
        p = pd.read_csv(root/'external/external_ft30_all_folds_predictions.csv', low_memory=False)
        m = pd.read_csv(root/'external/external_ft30_all_folds_metrics.csv', low_memory=False)
        _, audit = verified_external_metrics(m, p, split)
        audit.insert(0, 'lineage', label)
        checks.append(audit)
        # Also verify the unnumbered best-fold figure inputs independently.
        _, original = verified_external_metrics(
            pd.read_csv(root/'external/external_ft30_all_sample_metrics.csv', low_memory=False),
            pd.read_csv(root/'external/external_ft30_all_sample_predictions.csv', low_memory=False), split)
        cv = pd.read_csv(root/'benchmark/cv_fold_metrics.csv', low_memory=False)
        cp = pd.read_csv(root/'benchmark/cv_predictions.csv', low_memory=False)
        selected = cv[(cv.split=='Valid') & (cv.model=='ElasticNet')].sort_values(
            ['R2', 'MAE', 'K'], ascending=[False, True, True]).groupby('fold').head(1)
        for row in selected.itertuples():
            en = p[(p.model=='ElasticNet') & (p.K==row.K) & (p.fold==row.fold)]
            internal = cp[(cp.model=='ElasticNet') & (cp.K==row.K) & (cp.fold==row.fold)]
            groups = [(name, internal[internal.split==name]) for name in ['Train', 'Valid']]
            groups += [('External_all100', en), ('External_heldout70', en[en.is_heldout])]
            for name, part in groups:
                scores = independent_metrics(part.Age, part.prediction)
                internal_export_delta = 0.
                if name in ['Train', 'Valid']:
                    saved = cv[(cv.model=='ElasticNet') & (cv.K==row.K) & (cv.fold==row.fold) & (cv.split==name)].iloc[0]
                    for metric in METRICS:
                        # Archived internal predictions were exported from float32.
                        # Observed max EN discrepancy is 2.84e-6 years (MedAE);
                        # the external ALL100 guard above remains at 1e-8.
                        np.testing.assert_allclose(scores[metric], saved[metric], rtol=1e-8, atol=5e-6, equal_nan=True)
                        internal_export_delta=max(internal_export_delta,abs(scores[metric]-float(saved[metric])))
                diagnostics.append(dict(lineage=label, fold=row.fold, K=row.K, split=name, **scores,
                    internal_export_max_metric_difference=internal_export_delta,
                    K_selection='same-fold internal Valid R2; external not used'))
        summary[label] = dict(all_fold_groups=len(audit), best_fold_groups=len(original), n_all100=len(split),
                              n_adaptation=int((split.set=='adaptation').sum()), n_heldout=int((split.set=='holdout').sum()),
                              max_metric_difference=float(audit.max_saved_metric_difference.max()), status='passed')
        audit.to_csv(root/'reports/external_metric_scope_audit.csv', index=False)
    scope = pd.concat(checks, ignore_index=True)
    diag = pd.DataFrame(diagnostics).sort_values(['fold', 'lineage', 'split'])
    args.output.mkdir(parents=True, exist_ok=True)
    scope.to_csv(args.output/'external_metric_scope_audit.csv', index=False)
    diag.to_csv(args.output/'EN_all_folds_metric_diagnostics.csv', index=False)
    report = dict(status='passed', checked_at=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),
                  metrics=['R2','MAE','MAD','PCC','RMSE','MedAE'], external_scope='ALL100 for every primary metric',
                  independent_direct_formulas=True, mixed_70_100_detected=False, lineages=summary)
    (args.output/'external_metric_scope_audit.json').write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
    lines = ['# External 全量指标核验', '', '全部520组及旧版104组best-fold均通过。每组所有指标使用同一份100%外验年龄/预测配对，无70%/100%混用。', '',
             '`R² = 1 − Σ(y−预测)² / Σ(y−平均年龄)² = 1 − RMSE² / Age_SD²`；Age_SD使用ddof=0。',
             'MAE是绝对误差；R²还受年龄分布宽窄和平方误差影响。不同队列间MAE下降并不保证R²上升。R²不是PCC²。',
             'MAD = mean(abs(Age − mean(Age)))，是年龄标签离散度，不是模型误差；MedAE另报。', '',
             '每折E/F的ElasticNet K只由该折内部Valid R²选择。全量含拟合过的30%；70%留出指标独立保存，不进入主图。', '',
             '|折|谱系|K|分组|n|MAE|R²|PCC|MAD|RMSE|年龄SD|', '|---:|---|---:|---|---:|---:|---:|---:|---:|---:|---:|']
    for row in diag[diag.split.isin(['Valid','External_all100'])].itertuples():
        lines.append(f'|{row.fold}|{row.lineage}|{row.K}|{row.split}|{row.n}|{row.MAE:.5f}|{row.R2:.5f}|{row.PCC:.5f}|{row.MAD:.5f}|{row.RMSE:.5f}|{row.Age_SD:.5f}|')
    lines += ['', '逐组样本数/ID哈希及所有指标：external_metric_scope_audit.csv。',
              'EN训练/验证/全量外验/70%诊断：EN_all_folds_metric_diagnostics.csv。',
              '绘图入口强制核验全部样本集合并从同一配对向量重新计算所有External指标。原始模型、预测与内部CV保留。']
    (ROOT/'external_metric_scope_audit.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
