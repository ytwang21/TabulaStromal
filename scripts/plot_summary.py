#!/usr/bin/env python3
"""Render a scope-labelled figure directly from the released aggregate CSV."""
import csv
from pathlib import Path
import os

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / 'outputs/runtime/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    with (ROOT / 'results/final_clock_results_summary.csv').open() as handle:
        rows = list(csv.DictReader(handle))
    scopes = [('internal_original', 'Internal CV\n(selected K)'),
              ('noadapt', 'External\nno adaptation'),
              ('FT30_heldout70', 'External held-out 70%\nafter 30% adaptation')]
    colors = ['#cf635e', '#3981bd']
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for ax, metric, title in zip(axes, ['R2', 'MAE'], ['R²', 'MAE (years)']):
        for j, row in enumerate(rows):
            x = [i + (j - .5) * .35 for i in range(3)]
            values = [float(row[prefix + '_' + metric]) for prefix, _ in scopes]
            ax.bar(x, values, width=.32, color=colors[j], label=row['lineage'])
            for position, value in zip(x, values):
                label = f'{value:.4f}' if metric == 'R2' and 0 < abs(value) < .005 else f'{value:.2f}'
                ax.annotate(label, (position, value), ha='center', fontsize=9,
                            xytext=(0, 4 if value >= 0 else -12), textcoords='offset points')
        ax.set_xticks(range(3), [label for _, label in scopes], fontsize=9)
        ax.set_ylabel(title)
        ax.axhline(0, linewidth=.7, color='#445566')
        ax.spines[['top', 'right']].set_visible(False)
        ax.margins(y=.23)
    axes[0].legend(frameon=False, ncol=2, loc='lower left', fontsize=9)
    fig.suptitle('Final ElasticNet clocks: performance depends on evaluation scope', fontsize=13)
    fig.text(.5, .025, 'Source: released final_clock_results_summary.csv. Adaptation samples excluded from the held-out panel.',
             ha='center', fontsize=8, color='#445566')
    fig.tight_layout(rect=(0, .07, 1, .94))
    dest = ROOT / 'docs/assets/clock_performance.png'
    dest.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest, dpi=170)
    plt.close(fig)
    print(dest.relative_to(ROOT))


if __name__ == '__main__':
    main()
