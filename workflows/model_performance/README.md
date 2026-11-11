# F2_ModelPerformance workflow

Historical Figure 2 builder. External source-fold results use FT30 ALL100, including fitted samples.

[Original result figures and tables](../../F2_ModelPerformance/) · [Reproduction instructions](../../docs/reproducing.md)

![Related results](../../F2_ModelPerformance/figures/TabulaStromal_Workflow_V6.png)

# Historical Figure 2

`make_fig2.py` preserves the original fold-4 figure builder. It requires archived scatter predictions, metric diagnostics, source manifests, and an overview PDF; these inputs are not included in this source release. Its historical `SOURCE` path must point to the restored original figure archive. Do not infer a runnable asset bundle from the presence of this script.

Released aggregate data are under `results/historical_figure2/`. The external values use FT30 ALL100 (including adaptation samples). Fold-4 K=20,000 EN scatter plots do not represent the final full-development endothelial K=10,000 model. Historical “biological age” axis wording is not evidence of biological-age validation. README displays the supplied workflow overview followed by the original B–D heatmap panel exports. Its caption identifies the historical scope; the final-clock metrics are reported separately in the results table.
