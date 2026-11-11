# AUDIT ONLY: extracted reference source; do not run this file. Blood data paths are provenance, not project inputs.
# REFERENCE CELL 19
age_col="Age"
id_col="Individual"
prop_suffix="_ratio"
test_size=0.2
random_state=42
vt_threshold=1e-5
k_best=8000
bootstrap_n=300
stable_threshold=0.7
l1_ratio_grid=[0.5, 0.7, 0.9, 0.95, 0.99]
final_l1_ratio_grid=[0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1]
cv_splits=5
output_dir="fig3/iter_model/inte"
model_name="en_model_inte_bootstrap.pkl"

# REFERENCE CELL 22
# =========================
# 2. 划分训练集 / 验证集
# =========================
print("Splitting train and test sets...")
X_prop_train, X_prop_test, X_expr_train, X_expr_test, y_train, y_test = train_test_split(
    X_prop,
    X_expr,
    y,
    test_size=test_size,
    random_state=random_state
)

# REFERENCE CELL 23
# =========================
# 3. 训练集上做方差过滤 + SelectKBest
# =========================
print("Performing feature selection on training set...")
vt = VarianceThreshold(threshold=vt_threshold)
X_expr_train_vt = vt.fit_transform(X_expr_train)
X_expr_test_vt = vt.transform(X_expr_test)

k_real = min(k_best, X_expr_train_vt.shape[1])

skb = SelectKBest(score_func=f_regression, k=k_real)
X_expr_train_selected = skb.fit_transform(X_expr_train_vt, y_train)
X_expr_test_selected = skb.transform(X_expr_test_vt)

selected_expr_cols = X_expr_train.columns[vt.get_support()][skb.get_support()]

# REFERENCE CELL 24
# =========================
# 4. 标准化
# =========================
print("Standardizing features...")
scaler_prop = StandardScaler()
scaler_expr = StandardScaler()

X_prop_train_scaled = scaler_prop.fit_transform(X_prop_train)
X_prop_test_scaled = scaler_prop.transform(X_prop_test)

X_expr_train_scaled = scaler_expr.fit_transform(X_expr_train_selected)
X_expr_test_scaled = scaler_expr.transform(X_expr_test_selected)

all_cols = list(prop_cols) + list(selected_expr_cols)

X_train_final = np.hstack([X_prop_train_scaled, X_expr_train_scaled])
X_test_final = np.hstack([X_prop_test_scaled, X_expr_test_scaled])

X_train_final_df = pd.DataFrame(X_train_final, columns=all_cols, index=X_prop_train.index)
X_test_final_df = pd.DataFrame(X_test_final, columns=all_cols, index=X_prop_test.index)

# REFERENCE CELL 25
# =========================
# 5. 初步 ElasticNetCV 选择参数
# =========================
print("Performing initial ElasticNetCV for parameter tuning...")
base_model = ElasticNetCV(
    l1_ratio=l1_ratio_grid,
    cv=cv_splits,
    max_iter=5000,
    n_jobs=-1,
    random_state=random_state
)

base_model.fit(X_train_final_df, y_train)

best_alpha = base_model.alpha_
best_l1 = base_model.l1_ratio_

print("Best alpha:", best_alpha)
print("Best l1_ratio:", best_l1)

# REFERENCE CELL 26
# =========================
# 6. Bootstrap 特征稳定性筛选
# =========================
print("Performing bootstrap for feature stability selection...")
feature_counter = Counter()

for i in range(bootstrap_n):
    X_boot, y_boot = resample(
        X_train_final_df,
        y_train,
        replace=True,
        n_samples=len(y_train),
        random_state=i
    )

    model = ElasticNet(
        alpha=best_alpha,
        l1_ratio=best_l1,
        max_iter=5000,
        random_state=random_state
    )

    model.fit(X_boot, y_boot)

    selected = X_train_final_df.columns[model.coef_ != 0]

    for f in selected:
        feature_counter[f] += 1

print("Bootstrap finished.")

selection_freq = pd.DataFrame({
    "Feature": list(feature_counter.keys()),
    "Frequency": [feature_counter[f] / bootstrap_n for f in feature_counter]
}).sort_values("Frequency", ascending=False)

selection_freq.to_csv(
    os.path.join(output_dir, "bootstrap_selection_frequency.csv"),
    index=False
)

# REFERENCE CELL 27
# =========================
# 8. 获取稳定特征
# =========================
stable_features = selection_freq[
    selection_freq["Frequency"] >= stable_threshold
].copy()

selected_features = stable_features["Feature"].tolist()

if len(selected_features) == 0:
    raise ValueError(
        f"没有特征达到 stable_threshold={stable_threshold}，建议降低阈值，比如 0.5 或 0.6。"
    )

X_train_selected = X_train_final_df[selected_features]
X_test_selected = X_test_final_df[selected_features]

stable_features.to_csv(
    os.path.join(output_dir, "stable_features.csv"),
    index=False
)

print(f"Stable features number: {len(selected_features)}")


# REFERENCE CELL 28
# =========================
# 9. 用稳定特征重新训练最终模型
# =========================
print("Training final model with stable features...")
cv = KFold(
    n_splits=cv_splits,
    shuffle=True,
    random_state=random_state
)

model_reduced = ElasticNetCV(
    l1_ratio=final_l1_ratio_grid,
    alphas=None,
    cv=cv,
    max_iter=5000,
    n_jobs=-1,
    random_state=random_state
)

scores_reduced = cross_val_score(
    model_reduced,
    X_train_selected,
    y_train,
    cv=cv,
    scoring="r2"
)

print(
    f"训练集交叉验证 R²："
    f"{scores_reduced.mean():.3f} ± {scores_reduced.std():.3f}"
)

model_reduced.fit(X_train_selected, y_train)

# REFERENCE CELL 29
# =========================
# 10. 训练集和验证集预测
# =========================
print("Evaluating final model on train and test sets...")
y_train_pred = model_reduced.predict(X_train_selected)
y_test_pred = model_reduced.predict(X_test_selected)

train_r2 = r2_score(y_train, y_train_pred)
test_r2 = r2_score(y_test, y_test_pred)

train_mae = mean_absolute_error(y_train, y_train_pred)
test_mae = mean_absolute_error(y_test, y_test_pred)

train_r, train_p = pearsonr(y_train, y_train_pred)
test_r, test_p = pearsonr(y_test, y_test_pred)

metrics = {
    "best_alpha_initial": best_alpha,
    "best_l1_ratio_initial": best_l1,
    "final_alpha": model_reduced.alpha_,
    "final_l1_ratio": model_reduced.l1_ratio_,
    "n_stable_features": len(selected_features),
    "cv_r2_mean_train": scores_reduced.mean(),
    "cv_r2_std_train": scores_reduced.std(),
    "train_r2": train_r2,
    "train_mae": train_mae,
    "train_pearson_r": train_r,
    "train_pearson_p": train_p,
    "test_r2": test_r2,
    "test_mae": test_mae,
    "test_pearson_r": test_r,
    "test_pearson_p": test_p
}

metrics_df = pd.DataFrame([metrics])
metrics_df.to_csv(os.path.join(output_dir, "model_metrics.csv"), index=False)

print("Train R²:", round(train_r2, 3))
print("Test R²:", round(test_r2, 3))
print("Train MAE:", round(train_mae, 3))
print("Test MAE:", round(test_mae, 3))
print("Test Pearson r:", round(test_r, 3))

# REFERENCE CELL 32
# =========================
# 13. 保存最终模型和预处理器
# =========================
print("Saving final model and preprocessors...")
joblib.dump(
    {
        "vt": vt,
        "skb": skb,
        "scaler_prop": scaler_prop,
        "scaler_expr": scaler_expr,
        "model": model_reduced,
        "prop_cols": prop_cols,
        "expr_cols": expr_cols,
        "selected_expr_cols": list(selected_expr_cols),
        "selected_features": selected_features,
        "stable_features": stable_features,
        "metrics": metrics
    },
    os.path.join(output_dir, model_name)
)

# REFERENCE CELL 34
import os
import joblib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from collections import Counter
from sklearn.model_selection import train_test_split, KFold, cross_val_score
from sklearn.feature_selection import VarianceThreshold, SelectKBest, f_regression
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNetCV, ElasticNet
from sklearn.utils import resample
from sklearn.metrics import r2_score, mean_absolute_error
from scipy.stats import pearsonr


def build_aging_clock_elasticnet_bootstrap(
    df_all_features,
    age_col="Age",
    id_col="Individual",
    prop_suffix="_ratio",
    test_size=0.2,
    random_state=42,
    vt_threshold=1e-5,
    k_best=8000,
    bootstrap_n=300,
    stable_threshold=0.7,
    l1_ratio_grid=[0.5, 0.7, 0.9, 0.95, 0.99],
    final_l1_ratio_grid=[0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1],
    cv_splits=5,
    output_dir="external_inputs/inte_AIDA",
    model_name="en_model_inte_AIDA_bootstrap.pkl"
):
    os.makedirs(output_dir, exist_ok=True)

    # =========================
    # 1. 准备 X / y
    # =========================
    print("Preparing X and y...")
    y = df_all_features[age_col].values
    X = df_all_features.drop(columns=[age_col])

    prop_cols = [col for col in X.columns if col.endswith(prop_suffix)]
    expr_cols = [col for col in X.columns if col not in prop_cols and col != id_col]

    X_prop = X[prop_cols].copy()
    X_expr = X[expr_cols].copy()

    # =========================
    # 2. 划分训练集 / 验证集
    # =========================
    print("Splitting train and test sets...")
    X_prop_train, X_prop_test, X_expr_train, X_expr_test, y_train, y_test = train_test_split(
        X_prop,
        X_expr,
        y,
        test_size=test_size,
        random_state=random_state
    )

    # =========================
    # 3. 训练集上做方差过滤 + SelectKBest
    # =========================
    print("Performing feature selection on training set...")
    vt = VarianceThreshold(threshold=vt_threshold)
    X_expr_train_vt = vt.fit_transform(X_expr_train)
    X_expr_test_vt = vt.transform(X_expr_test)

    k_real = min(k_best, X_expr_train_vt.shape[1])

    skb = SelectKBest(score_func=f_regression, k=k_real)
    X_expr_train_selected = skb.fit_transform(X_expr_train_vt, y_train)
    X_expr_test_selected = skb.transform(X_expr_test_vt)

    selected_expr_cols = X_expr_train.columns[vt.get_support()][skb.get_support()]

    # =========================
    # 4. 标准化
    # =========================
    print("Standardizing features...")
    scaler_prop = StandardScaler()
    scaler_expr = StandardScaler()

    X_prop_train_scaled = scaler_prop.fit_transform(X_prop_train)
    X_prop_test_scaled = scaler_prop.transform(X_prop_test)

    X_expr_train_scaled = scaler_expr.fit_transform(X_expr_train_selected)
    X_expr_test_scaled = scaler_expr.transform(X_expr_test_selected)

    all_cols = list(prop_cols) + list(selected_expr_cols)

    X_train_final = np.hstack([X_prop_train_scaled, X_expr_train_scaled])
    X_test_final = np.hstack([X_prop_test_scaled, X_expr_test_scaled])

    X_train_final_df = pd.DataFrame(X_train_final, columns=all_cols, index=X_prop_train.index)
    X_test_final_df = pd.DataFrame(X_test_final, columns=all_cols, index=X_prop_test.index)

    # =========================
    # 5. 初步 ElasticNetCV 选择参数
    # =========================
    print("Performing initial ElasticNetCV for parameter tuning...")
    base_model = ElasticNetCV(
        l1_ratio=l1_ratio_grid,
        cv=cv_splits,
        max_iter=5000,
        n_jobs=-1,
        random_state=random_state
    )

    base_model.fit(X_train_final_df, y_train)

    best_alpha = base_model.alpha_
    best_l1 = base_model.l1_ratio_

    print("Best alpha:", best_alpha)
    print("Best l1_ratio:", best_l1)

    # =========================
    # 6. Bootstrap 特征稳定性筛选
    # =========================
    print("Performing bootstrap for feature stability selection...")
    feature_counter = Counter()

    for i in range(bootstrap_n):
        X_boot, y_boot = resample(
            X_train_final_df,
            y_train,
            replace=True,
            n_samples=len(y_train),
            random_state=i
        )

        model = ElasticNet(
            alpha=best_alpha,
            l1_ratio=best_l1,
            max_iter=5000,
            random_state=random_state
        )

        model.fit(X_boot, y_boot)

        selected = X_train_final_df.columns[model.coef_ != 0]

        for f in selected:
            feature_counter[f] += 1

    print("Bootstrap finished.")

    selection_freq = pd.DataFrame({
        "Feature": list(feature_counter.keys()),
        "Frequency": [feature_counter[f] / bootstrap_n for f in feature_counter]
    }).sort_values("Frequency", ascending=False)

    selection_freq.to_csv(
        os.path.join(output_dir, "bootstrap_selection_frequency.csv"),
        index=False
    )

    # =========================
    # 8. 获取稳定特征
    # =========================
    stable_features = selection_freq[
        selection_freq["Frequency"] >= stable_threshold
    ].copy()

    selected_features = stable_features["Feature"].tolist()

    if len(selected_features) == 0:
        raise ValueError(
            f"没有特征达到 stable_threshold={stable_threshold}，建议降低阈值，比如 0.5 或 0.6。"
        )

    X_train_selected = X_train_final_df[selected_features]
    X_test_selected = X_test_final_df[selected_features]

    stable_features.to_csv(
        os.path.join(output_dir, "stable_features.csv"),
        index=False
    )

    print(f"Stable features number: {len(selected_features)}")

    # =========================
    # 9. 用稳定特征重新训练最终模型
    # =========================
    print("Training final model with stable features...")
    cv = KFold(
        n_splits=cv_splits,
        shuffle=True,
        random_state=random_state
    )

    model_reduced = ElasticNetCV(
        l1_ratio=final_l1_ratio_grid,
        alphas=None,
        cv=cv,
        max_iter=5000,
        n_jobs=-1,
        random_state=random_state
    )

    scores_reduced = cross_val_score(
        model_reduced,
        X_train_selected,
        y_train,
        cv=cv,
        scoring="r2"
    )

    print(
        f"训练集交叉验证 R²："
        f"{scores_reduced.mean():.3f} ± {scores_reduced.std():.3f}"
    )

    model_reduced.fit(X_train_selected, y_train)

    # =========================
    # 10. 训练集和验证集预测
    # =========================
    print("Evaluating final model on train and test sets...")
    y_train_pred = model_reduced.predict(X_train_selected)
    y_test_pred = model_reduced.predict(X_test_selected)

    train_r2 = r2_score(y_train, y_train_pred)
    test_r2 = r2_score(y_test, y_test_pred)

    train_mae = mean_absolute_error(y_train, y_train_pred)
    test_mae = mean_absolute_error(y_test, y_test_pred)

    train_r, train_p = pearsonr(y_train, y_train_pred)
    test_r, test_p = pearsonr(y_test, y_test_pred)

    metrics = {
        "best_alpha_initial": best_alpha,
        "best_l1_ratio_initial": best_l1,
        "final_alpha": model_reduced.alpha_,
        "final_l1_ratio": model_reduced.l1_ratio_,
        "n_stable_features": len(selected_features),
        "cv_r2_mean_train": scores_reduced.mean(),
        "cv_r2_std_train": scores_reduced.std(),
        "train_r2": train_r2,
        "train_mae": train_mae,
        "train_pearson_r": train_r,
        "train_pearson_p": train_p,
        "test_r2": test_r2,
        "test_mae": test_mae,
        "test_pearson_r": test_r,
        "test_pearson_p": test_p
    }

    metrics_df = pd.DataFrame([metrics])
    metrics_df.to_csv(os.path.join(output_dir, "model_metrics.csv"), index=False)

    print("Train R²:", round(train_r2, 3))
    print("Test R²:", round(test_r2, 3))
    print("Train MAE:", round(train_mae, 3))
    print("Test MAE:", round(test_mae, 3))
    print("Test Pearson r:", round(test_r, 3))


    # =========================
    # 12. 训练集 + 验证集画到同一张图
    # =========================
    print("Plotting train and test predictions...")

    def plot_prediction_train_test(
        y_train,
        y_train_pred,
        y_test,
        y_test_pred,
        title=None,
        test_label="Validation",
        save_path=None,
        show=True
    ):
        import numpy as np
        import matplotlib.pyplot as plt

        from scipy.stats import pearsonr
        from sklearn.metrics import r2_score, mean_absolute_error

        # =========================
        # 1. 数据整理
        # =========================
        y_train = np.asarray(y_train).flatten()
        y_train_pred = np.asarray(y_train_pred).flatten()
        y_test = np.asarray(y_test).flatten()
        y_test_pred = np.asarray(y_test_pred).flatten()

        # 去除无效值
        train_mask = np.isfinite(y_train) & np.isfinite(y_train_pred)
        test_mask = np.isfinite(y_test) & np.isfinite(y_test_pred)

        y_train = y_train[train_mask]
        y_train_pred = y_train_pred[train_mask]

        y_test = y_test[test_mask]
        y_test_pred = y_test_pred[test_mask]

        # =========================
        # 2. 计算评价指标
        # =========================
        train_r2 = r2_score(y_train, y_train_pred)
        train_mae = mean_absolute_error(y_train, y_train_pred)
        train_mad = np.median(np.abs(y_train - y_train_pred))
        train_pcc = pearsonr(y_train, y_train_pred)[0]

        test_r2 = r2_score(y_test, y_test_pred)
        test_mae = mean_absolute_error(y_test, y_test_pred)
        test_mad = np.median(np.abs(y_test - y_test_pred))
        test_pcc = pearsonr(y_test, y_test_pred)[0]

        # =========================
        # 3. 论文风格基础设置
        # =========================
        train_color = "#4C78A8"      # 蓝色
        test_color = "#E76F51"       # 橙红
        ideal_color = "#8C8C8C"      # 灰色 y=x
        regression_color = "#CC4C4C" # 深红回归线

        with plt.rc_context({
            "font.family": "Arial",
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 8.5,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "xtick.major.size": 3.5,
            "ytick.major.size": 3.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42
        }):

            # 适合论文多面板拼图的比例
            fig, ax = plt.subplots(figsize=(4.5, 4.2))

            # =========================
            # 4. 固定坐标范围（论文统一）
            # =========================
            xmin = 0
            xmax = 120

            ax.set_xlim(xmin, xmax)
            ax.set_ylim(xmin, xmax)
            ax.set_aspect("equal", adjustable="box")

            # 每20岁一个刻度
            ticks = np.arange(0, 121, 20)
            ax.set_xticks(ticks)
            ax.set_yticks(ticks)

            # =========================
            # 5. 绘制散点
            # =========================
            ax.scatter(
                y_train,
                y_train_pred,
                s=22,
                color=train_color,
                alpha=0.55,
                edgecolors="white",
                linewidths=0.25,
                label=f"Train, n={len(y_train)}",
                zorder=2,
                rasterized=True
            )

            ax.scatter(
                y_test,
                y_test_pred,
                s=27,
                color=test_color,
                alpha=0.80,
                edgecolors="white",
                linewidths=0.3,
                label=f"{test_label}, n={len(y_test)}",
                zorder=3,
                rasterized=True
            )

            # =========================
            # 6. 理想线 y=x
            # =========================
            ax.plot(
                [xmin, xmax],
                [xmin, xmax],
                linestyle="--",
                color=ideal_color,
                linewidth=1.2,
                dashes=(4, 3),
                label=r"Ideal line, $y=x$",
                zorder=1
            )

            # =========================
            # 7. Validation/Test回归线
            # =========================
            slope, intercept = np.polyfit(
                y_test,
                y_test_pred,
                deg=1
            )

            x_line = np.linspace(xmin, xmax, 300)
            y_line = slope * x_line + intercept

            ax.plot(
                x_line,
                y_line,
                color=regression_color,
                linewidth=1.8,
                label=f"{test_label} regression",
                zorder=4
            )

            # =========================
            # 8. 指标文本
            # =========================
            metric_text = (
                r"$\bf{Train}$" + "\n"
                + rf"$R^2$ = {train_r2:.2f}" + "\n"
                + f"MAE = {train_mae:.2f}" + "\n"
                + f"MAD = {train_mad:.2f}" + "\n"
                + f"PCC = {train_pcc:.2f}" + "\n\n"
                + rf"$\bf{{{test_label}}}$" + "\n"
                + rf"$R^2$ = {test_r2:.2f}" + "\n"
                + f"MAE = {test_mae:.2f}" + "\n"
                + f"MAD = {test_mad:.2f}" + "\n"
                + f"PCC = {test_pcc:.2f}"
            )

            ax.text(
                0.04,
                0.96,
                metric_text,
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=8.5,
                linespacing=1.15,
                bbox=dict(
                    boxstyle="round,pad=0.35",
                    facecolor="white",
                    edgecolor="#B5B5B5",
                    linewidth=0.6,
                    alpha=0.92
                ),
                zorder=10
            )

            # =========================
            # 9. 坐标轴与标题
            # =========================
            ax.set_xlabel("Chronological age (years)")
            ax.set_ylabel("Predicted age (years)")

            if title is not None:
                ax.set_title(
                    title,
                    fontsize=12,
                    fontweight="bold",
                    pad=8
                )

            # 去掉顶部和右侧边框
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            ax.spines["left"].set_color("#333333")
            ax.spines["bottom"].set_color("#333333")

            ax.tick_params(
                axis="both",
                direction="out",
                colors="#333333"
            )

            # 使用很弱的辅助网格
            ax.grid(
                True,
                linestyle="-",
                linewidth=0.4,
                alpha=0.16,
                color="#808080"
            )

            ax.set_axisbelow(True)

            # =========================
            # 10. 图例
            # =========================
            legend = ax.legend(
                loc="lower right",
                frameon=False,
                handlelength=2.3,
                handletextpad=0.6,
                borderaxespad=0.6,
                labelspacing=0.5
            )

            plt.tight_layout(pad=0.6)

            # =========================
            # 11. 保存
            # =========================
            if save_path is not None:
                fig.savefig(
                    os.path.join(save_path, "ElasticNet_train_validation.png"),
                    dpi=600,
                    bbox_inches="tight",
                    facecolor="white"
                )

                fig.savefig(
                    os.path.join(save_path, "ElasticNet_train_validation.pdf"),
                    dpi=600,
                    bbox_inches="tight",
                    facecolor="white"
                )

            if show:
                plt.show()
            else:
                plt.close(fig)

        # 返回评价指标，方便后续整理结果表
        metrics = {
            "Train_R2": train_r2,
            "Train_MAE": train_mae,
            "Train_MAD": train_mad,
            "Train_PCC": train_pcc,
            f"{test_label}_R2": test_r2,
            f"{test_label}_MAE": test_mae,
            f"{test_label}_MAD": test_mad,
            f"{test_label}_PCC": test_pcc,
            f"{test_label}_Slope": slope,
            f"{test_label}_Intercept": intercept
        }

        return fig, ax, metrics

    fig, ax, prediction_metrics = plot_prediction_train_test(
        y_train=y_train,
        y_train_pred=y_train_pred,
        y_test=y_test,
        y_test_pred=y_test_pred,
        title=None,
        test_label="Validation",
        save_path=output_dir
)

    # =========================
    # 13. 保存最终模型和预处理器
    # =========================
    print("Saving final model and preprocessors...")
    joblib.dump(
        {
            "vt": vt,
            "skb": skb,
            "scaler_prop": scaler_prop,
            "scaler_expr": scaler_expr,
            "model": model_reduced,
            "prop_cols": prop_cols,
            "expr_cols": expr_cols,
            "selected_expr_cols": list(selected_expr_cols),
            "selected_features": selected_features,
            "stable_features": stable_features,
            "metrics": metrics
        },
        os.path.join(output_dir, model_name)
    )

    print(f"Model saved to: {os.path.join(output_dir, model_name)}")

    return {
        "model": model_reduced,
        "metrics": metrics_df,
        "selection_freq": selection_freq,
        "stable_features": stable_features,
        "selected_features": selected_features,
        "X_train_selected": X_train_selected,
        "X_test_selected": X_test_selected,
        "y_train": y_train,
        "y_test": y_test,
        "y_train_pred": y_train_pred,
        "y_test_pred": y_test_pred
    }

# REFERENCE CELL 36
result = build_aging_clock_elasticnet_bootstrap(
    df_all_features=df_all_features,
    age_col="Age",
    id_col="Individual",
    prop_suffix="_ratio",
    test_size=0.2,
    bootstrap_n=300,
    stable_threshold=0.7,
    output_dir="external_inputs/inte_AIDA"
)

# REFERENCE CELL 39
import os
import joblib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from collections import Counter
from sklearn.model_selection import train_test_split, KFold, cross_val_score
from sklearn.feature_selection import VarianceThreshold, SelectKBest, f_regression
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNetCV, ElasticNet
from sklearn.utils import resample
from sklearn.metrics import r2_score, mean_absolute_error
from scipy.stats import pearsonr


def build_aging_clock_elasticnet_bootstrap(
    df_all_features,
    age_col="Age",
    id_col="Individual",
    prop_suffix="_ratio",
    test_size=0.2,
    random_state=42,
    vt_threshold=1e-5,
    k_best=8000,
    bootstrap_n=300,
    stable_threshold=0.7,
    l1_ratio_grid=[0.5, 0.7, 0.9, 0.95, 0.99],
    final_l1_ratio_grid=[0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1],
    cv_splits=5,
    output_dir="external_inputs/inte_AIDA_CIMA",
    model_name="en_model_inte_AIDA_CIMA_bootstrap.pkl"
):
    os.makedirs(output_dir, exist_ok=True)

    # =========================
    # 1. 准备 X / y
    # =========================
    print("Preparing X and y...")
    y = df_all_features[age_col].values
    X = df_all_features.drop(columns=[age_col])

    prop_cols = [col for col in X.columns if col.endswith(prop_suffix)]
    expr_cols = [col for col in X.columns if col not in prop_cols and col != id_col]

    X_prop = X[prop_cols].copy()
    X_expr = X[expr_cols].copy()

    # =========================
    # 2. 划分训练集 / 验证集
    # =========================
    print("Splitting train and test sets...")
    X_prop_train, X_prop_test, X_expr_train, X_expr_test, y_train, y_test = train_test_split(
        X_prop,
        X_expr,
        y,
        test_size=test_size,
        random_state=random_state
    )

    # =========================
    # 3. 训练集上做方差过滤 + SelectKBest
    # =========================
    print("Performing feature selection on training set...")
    vt = VarianceThreshold(threshold=vt_threshold)
    X_expr_train_vt = vt.fit_transform(X_expr_train)
    X_expr_test_vt = vt.transform(X_expr_test)

    k_real = min(k_best, X_expr_train_vt.shape[1])

    skb = SelectKBest(score_func=f_regression, k=k_real)
    X_expr_train_selected = skb.fit_transform(X_expr_train_vt, y_train)
    X_expr_test_selected = skb.transform(X_expr_test_vt)

    selected_expr_cols = X_expr_train.columns[vt.get_support()][skb.get_support()]

    # =========================
    # 4. 标准化
    # =========================
    print("Standardizing features...")
    scaler_prop = StandardScaler()
    scaler_expr = StandardScaler()

    X_prop_train_scaled = scaler_prop.fit_transform(X_prop_train)
    X_prop_test_scaled = scaler_prop.transform(X_prop_test)

    X_expr_train_scaled = scaler_expr.fit_transform(X_expr_train_selected)
    X_expr_test_scaled = scaler_expr.transform(X_expr_test_selected)

    all_cols = list(prop_cols) + list(selected_expr_cols)

    X_train_final = np.hstack([X_prop_train_scaled, X_expr_train_scaled])
    X_test_final = np.hstack([X_prop_test_scaled, X_expr_test_scaled])

    X_train_final_df = pd.DataFrame(X_train_final, columns=all_cols, index=X_prop_train.index)
    X_test_final_df = pd.DataFrame(X_test_final, columns=all_cols, index=X_prop_test.index)

    # =========================
    # 5. 初步 ElasticNetCV 选择参数
    # =========================
    print("Performing initial ElasticNetCV for parameter tuning...")
    base_model = ElasticNetCV(
        l1_ratio=l1_ratio_grid,
        cv=cv_splits,
        max_iter=5000,
        n_jobs=-1,
        random_state=random_state
    )

    base_model.fit(X_train_final_df, y_train)

    best_alpha = base_model.alpha_
    best_l1 = base_model.l1_ratio_

    print("Best alpha:", best_alpha)
    print("Best l1_ratio:", best_l1)

    # =========================
    # 6. Bootstrap 特征稳定性筛选
    # =========================
    print("Performing bootstrap for feature stability selection...")
    feature_counter = Counter()

    for i in range(bootstrap_n):
        X_boot, y_boot = resample(
            X_train_final_df,
            y_train,
            replace=True,
            n_samples=len(y_train),
            random_state=i
        )

        model = ElasticNet(
            alpha=best_alpha,
            l1_ratio=best_l1,
            max_iter=5000,
            random_state=random_state
        )

        model.fit(X_boot, y_boot)

        selected = X_train_final_df.columns[model.coef_ != 0]

        for f in selected:
            feature_counter[f] += 1

    print("Bootstrap finished.")

    selection_freq = pd.DataFrame({
        "Feature": list(feature_counter.keys()),
        "Frequency": [feature_counter[f] / bootstrap_n for f in feature_counter]
    }).sort_values("Frequency", ascending=False)

    selection_freq.to_csv(
        os.path.join(output_dir, "bootstrap_selection_frequency.csv"),
        index=False
    )

    # =========================
    # 8. 获取稳定特征
    # =========================
    stable_features = selection_freq[
        selection_freq["Frequency"] >= stable_threshold
    ].copy()

    selected_features = stable_features["Feature"].tolist()

    if len(selected_features) == 0:
        raise ValueError(
            f"没有特征达到 stable_threshold={stable_threshold}，建议降低阈值，比如 0.5 或 0.6。"
        )

    X_train_selected = X_train_final_df[selected_features]
    X_test_selected = X_test_final_df[selected_features]

    stable_features.to_csv(
        os.path.join(output_dir, "stable_features.csv"),
        index=False
    )

    print(f"Stable features number: {len(selected_features)}")

    # =========================
    # 9. 用稳定特征重新训练最终模型
    # =========================
    print("Training final model with stable features...")
    cv = KFold(
        n_splits=cv_splits,
        shuffle=True,
        random_state=random_state
    )

    model_reduced = ElasticNetCV(
        l1_ratio=final_l1_ratio_grid,
        alphas=None,
        cv=cv,
        max_iter=5000,
        n_jobs=-1,
        random_state=random_state
    )

    scores_reduced = cross_val_score(
        model_reduced,
        X_train_selected,
        y_train,
        cv=cv,
        scoring="r2"
    )

    print(
        f"训练集交叉验证 R²："
        f"{scores_reduced.mean():.3f} ± {scores_reduced.std():.3f}"
    )

    model_reduced.fit(X_train_selected, y_train)

    # =========================
    # 10. 训练集和验证集预测
    # =========================
    print("Evaluating final model on train and test sets...")
    y_train_pred = model_reduced.predict(X_train_selected)
    y_test_pred = model_reduced.predict(X_test_selected)

    train_r2 = r2_score(y_train, y_train_pred)
    test_r2 = r2_score(y_test, y_test_pred)

    train_mae = mean_absolute_error(y_train, y_train_pred)
    test_mae = mean_absolute_error(y_test, y_test_pred)

    train_r, train_p = pearsonr(y_train, y_train_pred)
    test_r, test_p = pearsonr(y_test, y_test_pred)

    metrics = {
        "best_alpha_initial": best_alpha,
        "best_l1_ratio_initial": best_l1,
        "final_alpha": model_reduced.alpha_,
        "final_l1_ratio": model_reduced.l1_ratio_,
        "n_stable_features": len(selected_features),
        "cv_r2_mean_train": scores_reduced.mean(),
        "cv_r2_std_train": scores_reduced.std(),
        "train_r2": train_r2,
        "train_mae": train_mae,
        "train_pearson_r": train_r,
        "train_pearson_p": train_p,
        "test_r2": test_r2,
        "test_mae": test_mae,
        "test_pearson_r": test_r,
        "test_pearson_p": test_p
    }

    metrics_df = pd.DataFrame([metrics])
    metrics_df.to_csv(os.path.join(output_dir, "model_metrics.csv"), index=False)

    print("Train R²:", round(train_r2, 3))
    print("Test R²:", round(test_r2, 3))
    print("Train MAE:", round(train_mae, 3))
    print("Test MAE:", round(test_mae, 3))
    print("Test Pearson r:", round(test_r, 3))


    # =========================
    # 12. 训练集 + 验证集画到同一张图
    # =========================
    print("Plotting train and test predictions...")

    def plot_prediction_train_test(
        y_train,
        y_train_pred,
        y_test,
        y_test_pred,
        title=None,
        test_label="Validation",
        save_path=None,
        show=True
    ):
        import numpy as np
        import matplotlib.pyplot as plt

        from scipy.stats import pearsonr
        from sklearn.metrics import r2_score, mean_absolute_error

        # =========================
        # 1. 数据整理
        # =========================
        y_train = np.asarray(y_train).flatten()
        y_train_pred = np.asarray(y_train_pred).flatten()
        y_test = np.asarray(y_test).flatten()
        y_test_pred = np.asarray(y_test_pred).flatten()

        # 去除无效值
        train_mask = np.isfinite(y_train) & np.isfinite(y_train_pred)
        test_mask = np.isfinite(y_test) & np.isfinite(y_test_pred)

        y_train = y_train[train_mask]
        y_train_pred = y_train_pred[train_mask]

        y_test = y_test[test_mask]
        y_test_pred = y_test_pred[test_mask]

        # =========================
        # 2. 计算评价指标
        # =========================
        train_r2 = r2_score(y_train, y_train_pred)
        train_mae = mean_absolute_error(y_train, y_train_pred)
        train_mad = np.median(np.abs(y_train - y_train_pred))
        train_pcc = pearsonr(y_train, y_train_pred)[0]

        test_r2 = r2_score(y_test, y_test_pred)
        test_mae = mean_absolute_error(y_test, y_test_pred)
        test_mad = np.median(np.abs(y_test - y_test_pred))
        test_pcc = pearsonr(y_test, y_test_pred)[0]

        # =========================
        # 3. 论文风格基础设置
        # =========================
        train_color = "#4C78A8"      # 蓝色
        test_color = "#E76F51"       # 橙红
        ideal_color = "#8C8C8C"      # 灰色 y=x
        regression_color = "#CC4C4C" # 深红回归线

        with plt.rc_context({
            "font.family": "Arial",
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 8.5,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "xtick.major.size": 3.5,
            "ytick.major.size": 3.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42
        }):

            # 适合论文多面板拼图的比例
            fig, ax = plt.subplots(figsize=(4.5, 4.2))

            # =========================
            # 4. 固定坐标范围（论文统一）
            # =========================
            xmin = 0
            xmax = 120

            ax.set_xlim(xmin, xmax)
            ax.set_ylim(xmin, xmax)
            ax.set_aspect("equal", adjustable="box")

            # 每20岁一个刻度
            ticks = np.arange(0, 121, 20)
            ax.set_xticks(ticks)
            ax.set_yticks(ticks)

            # =========================
            # 5. 绘制散点
            # =========================
            ax.scatter(
                y_train,
                y_train_pred,
                s=22,
                color=train_color,
                alpha=0.55,
                edgecolors="white",
                linewidths=0.25,
                label=f"Train, n={len(y_train)}",
                zorder=2,
                rasterized=True
            )

            ax.scatter(
                y_test,
                y_test_pred,
                s=27,
                color=test_color,
                alpha=0.80,
                edgecolors="white",
                linewidths=0.3,
                label=f"{test_label}, n={len(y_test)}",
                zorder=3,
                rasterized=True
            )

            # =========================
            # 6. 理想线 y=x
            # =========================
            ax.plot(
                [xmin, xmax],
                [xmin, xmax],
                linestyle="--",
                color=ideal_color,
                linewidth=1.2,
                dashes=(4, 3),
                label=r"Ideal line, $y=x$",
                zorder=1
            )

            # =========================
            # 7. Validation/Test回归线
            # =========================
            slope, intercept = np.polyfit(
                y_test,
                y_test_pred,
                deg=1
            )

            x_line = np.linspace(xmin, xmax, 300)
            y_line = slope * x_line + intercept

            ax.plot(
                x_line,
                y_line,
                color=regression_color,
                linewidth=1.8,
                label=f"{test_label} regression",
                zorder=4
            )

            # =========================
            # 8. 指标文本
            # =========================
            metric_text = (
                r"$\bf{Train}$" + "\n"
                + rf"$R^2$ = {train_r2:.2f}" + "\n"
                + f"MAE = {train_mae:.2f}" + "\n"
                + f"MAD = {train_mad:.2f}" + "\n"
                + f"PCC = {train_pcc:.2f}" + "\n\n"
                + rf"$\bf{{{test_label}}}$" + "\n"
                + rf"$R^2$ = {test_r2:.2f}" + "\n"
                + f"MAE = {test_mae:.2f}" + "\n"
                + f"MAD = {test_mad:.2f}" + "\n"
                + f"PCC = {test_pcc:.2f}"
            )

            ax.text(
                0.04,
                0.96,
                metric_text,
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=8.5,
                linespacing=1.15,
                bbox=dict(
                    boxstyle="round,pad=0.35",
                    facecolor="white",
                    edgecolor="#B5B5B5",
                    linewidth=0.6,
                    alpha=0.92
                ),
                zorder=10
            )

            # =========================
            # 9. 坐标轴与标题
            # =========================
            ax.set_xlabel("Chronological age (years)")
            ax.set_ylabel("Predicted age (years)")

            if title is not None:
                ax.set_title(
                    title,
                    fontsize=12,
                    fontweight="bold",
                    pad=8
                )

            # 去掉顶部和右侧边框
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            ax.spines["left"].set_color("#333333")
            ax.spines["bottom"].set_color("#333333")

            ax.tick_params(
                axis="both",
                direction="out",
                colors="#333333"
            )

            # 使用很弱的辅助网格
            ax.grid(
                True,
                linestyle="-",
                linewidth=0.4,
                alpha=0.16,
                color="#808080"
            )

            ax.set_axisbelow(True)

            # =========================
            # 10. 图例
            # =========================
            legend = ax.legend(
                loc="lower right",
                frameon=False,
                handlelength=2.3,
                handletextpad=0.6,
                borderaxespad=0.6,
                labelspacing=0.5
            )

            plt.tight_layout(pad=0.6)

            # =========================
            # 11. 保存
            # =========================
            if save_path is not None:
                fig.savefig(
                    os.path.join(save_path, "ElasticNet_train_validation.png"),
                    dpi=600,
                    bbox_inches="tight",
                    facecolor="white"
                )

                fig.savefig(
                    os.path.join(save_path, "ElasticNet_train_validation.pdf"),
                    dpi=600,
                    bbox_inches="tight",
                    facecolor="white"
                )

            if show:
                plt.show()
            else:
                plt.close(fig)

        # 返回评价指标，方便后续整理结果表
        metrics = {
            "Train_R2": train_r2,
            "Train_MAE": train_mae,
            "Train_MAD": train_mad,
            "Train_PCC": train_pcc,
            f"{test_label}_R2": test_r2,
            f"{test_label}_MAE": test_mae,
            f"{test_label}_MAD": test_mad,
            f"{test_label}_PCC": test_pcc,
            f"{test_label}_Slope": slope,
            f"{test_label}_Intercept": intercept
        }

        return fig, ax, metrics

    fig, ax, prediction_metrics = plot_prediction_train_test(
        y_train=y_train,
        y_train_pred=y_train_pred,
        y_test=y_test,
        y_test_pred=y_test_pred,
        title=None,
        test_label="Validation",
        save_path=output_dir
)

    # =========================
    # 13. 保存最终模型和预处理器
    # =========================
    print("Saving final model and preprocessors...")
    joblib.dump(
        {
            "vt": vt,
            "skb": skb,
            "scaler_prop": scaler_prop,
            "scaler_expr": scaler_expr,
            "model": model_reduced,
            "prop_cols": prop_cols,
            "expr_cols": expr_cols,
            "selected_expr_cols": list(selected_expr_cols),
            "selected_features": selected_features,
            "stable_features": stable_features,
            "metrics": metrics
        },
        os.path.join(output_dir, model_name)
    )

    print(f"Model saved to: {os.path.join(output_dir, model_name)}")

    return {
        "model": model_reduced,
        "metrics": metrics_df,
        "selection_freq": selection_freq,
        "stable_features": stable_features,
        "selected_features": selected_features,
        "X_train_selected": X_train_selected,
        "X_test_selected": X_test_selected,
        "y_train": y_train,
        "y_test": y_test,
        "y_train_pred": y_train_pred,
        "y_test_pred": y_test_pred
    }

# REFERENCE CELL 40
result = build_aging_clock_elasticnet_bootstrap(
    df_all_features=df_all_features,
    age_col="Age",
    id_col="Individual",
    prop_suffix="_ratio",
    test_size=0.2,
    bootstrap_n=300,
    stable_threshold=0.7,
    output_dir="external_inputs/inte_AIDA_CIMA"
)

# REFERENCE CELL 51
# 目标变量 y 和特征矩阵 X
y = df_all_features['Age'].values
X = df_all_features.drop(columns='Age')

# 区分细胞比例特征（假设以 "_prop" 结尾）和表达特征（其余）
prop_cols = [col for col in X.columns if col.endswith('_ratio')]
expr_cols = [col for col in X.columns if col not in prop_cols and col != 'Individual']
expr_cols

X_prop = X[prop_cols].copy()
X_expr = X[expr_cols].copy()

# 去除低方差表达特征
vt = VarianceThreshold(threshold=1e-5)
X_expr_vt = vt.fit_transform(X_expr)

# 选择最相关的前 k 个特征
skb = SelectKBest(score_func=f_regression, k=8000)   # k=8000
X_expr_selected = skb.fit_transform(X_expr_vt, y)

# 获取保留下来的列名
selected_expr_cols = X_expr.columns[vt.get_support()][skb.get_support()]

# 合并，标准化
scaler_prop = StandardScaler()
scaler_expr = StandardScaler()

X_prop_scaled = scaler_prop.fit_transform(X_prop)
X_expr_scaled = scaler_expr.fit_transform(X_expr_selected)

# 合并，带列名
all_cols = list(X_prop.columns) + list(selected_expr_cols)
X_final = np.hstack([X_prop_scaled, X_expr_scaled])
X_final_df = pd.DataFrame(X_final, columns=all_cols)
X_final_df

# REFERENCE CELL 52
# 选择最好的参数
from sklearn.linear_model import ElasticNetCV

base_model = ElasticNetCV(
    l1_ratio=[0.5, 0.7, 0.9, 0.95, 0.99],
    cv=5,
    max_iter=5000,
    n_jobs=-1
)

base_model.fit(X_final_df, y)

best_alpha = base_model.alpha_
best_l1 = base_model.l1_ratio_

print("Best alpha:", best_alpha)
print("Best l1_ratio:", best_l1)

# REFERENCE CELL 53
# bootstrap 计算频率
from sklearn.linear_model import ElasticNet
from sklearn.utils import resample
from collections import Counter
import pandas as pd
import numpy as np

n_boot = 300   # 200 已经足够稳定
feature_counter = Counter()

for i in range(n_boot):

    X_boot, y_boot = resample(
        X_final_df,
        y,
        replace=True,
        n_samples=len(y),
        random_state=i
    )

    model = ElasticNet(
        alpha=best_alpha,
        l1_ratio=best_l1,
        max_iter=5000
    )

    model.fit(X_boot, y_boot)

    selected = X_final_df.columns[model.coef_ != 0]

    for f in selected:
        feature_counter[f] += 1

print("Bootstrap finished.")
# 计算选择频率
selection_freq = pd.DataFrame({
    "Feature": list(feature_counter.keys()),
    "Frequency": [feature_counter[f]/n_boot for f in feature_counter]
})

selection_freq = selection_freq.sort_values("Frequency", ascending=False)

# REFERENCE CELL 55
# 上面是为了画图又做了一遍bootstrap，下面直接用已经筛选好的特征
stable_features = pd.read_csv("external_inputs/bootstrap_selected_stable_features.csv")
selected_features = stable_features['Feature'].tolist()
X_selected = X_final_df[selected_features]

# 3. 定义交叉验证
cv = KFold(n_splits=5, shuffle=True, random_state=42)

# 4. 训练模型并评估
model_reduced = ElasticNetCV(
    l1_ratio=[.1, .5, .7, .9, .95, .99, 1],
    alphas=None,
    cv=cv,
    max_iter=5000,
    n_jobs=-1
)

# 交叉验证性能
scores_reduced = cross_val_score(model_reduced, X_selected, y, cv=cv, scoring='r2')
print(f"使用stable的173特征的交叉验证 R² 得分（平均±标准差）：{scores_reduced.mean():.3f} ± {scores_reduced.std():.3f}")
