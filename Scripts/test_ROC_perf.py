#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import auc, roc_curve
from statsmodels.nonparametric.smoothers_lowess import lowess

def smooth_and_export(
    x, y,
    frac=0.2,
    out_prefix="smoothed_curve",
    highlight_x=None   # <-- value to highlight
):
    x = np.asarray(x)
    y = np.asarray(y)

    # sort
    idx = np.argsort(x)
    x_sorted = x[idx]
    y_sorted = y[idx]

    # LOWESS smoothing
    smoothed = lowess(y_sorted, x_sorted, frac=frac)
    x_smooth = smoothed[:, 0]
    y_smooth = smoothed[:, 1]

    # ---- STORE ----
    df = pd.DataFrame({
        "x": x_smooth,
        "y_smooth": y_smooth
    })

    # ---- OPTIONAL HIGHLIGHT ----
    highlight_point = None
    if highlight_x is not None:
        # find closest x
        idx_closest = np.argmin(np.abs(x_smooth - highlight_x))
        hx = x_smooth[idx_closest]
        hy = y_smooth[idx_closest]

        highlight_point = (hx, hy)

        # add to dataframe (single-row export)
        df_highlight = pd.DataFrame({"x": [hx], "y_smooth": [hy]})
        df_highlight.to_csv(f"{out_prefix}_highlight.csv", index=False)

    # ---- EXPORT CSV ----
    csv_path = f"{out_prefix}.csv"
    df.to_csv(csv_path, index=False)

    # ---- PLOT ----
    plt.figure()
    plt.scatter(x, y, alpha=0.4, label="raw")
    plt.plot(x_smooth, y_smooth, label="smoothed")

    # highlight point on plot
    if highlight_point is not None:
        hx, hy = highlight_point
        plt.scatter([hx], [hy], s=80, zorder=3, label="highlight")
        plt.annotate(f"({hx:.3g}, {hy:.3g})",
                     (hx, hy),
                     textcoords="offset points",
                     xytext=(5,5))

    plt.legend()

    # ---- EXPORT IMAGE ----
    img_path = f"{out_prefix}.png"
    plt.savefig(img_path, dpi=300, bbox_inches="tight")
    plt.close()

    return df, highlight_point, csv_path, img_path

def empirical_pvalues(neg_values, values):
    neg_sorted = np.sort(neg_values)
    n = len(neg_sorted)

    def pval(x):
        left = np.searchsorted(neg_sorted, x, side='left') / n
        right = 1 - np.searchsorted(neg_sorted, x, side='right') / n
        return 2 * min(left, right)

    return np.array([pval(v) for v in values])


def evaluate_quantile_thresholds(df, lfc_col):
    df = df.copy()
    df[lfc_col] = pd.to_numeric(df[lfc_col])

    df = df[df['Controls'].isin(['negative_control', 'positive_control'])]

    df['y_true'] = df['Controls'] == 'positive_control'

    neg = df.loc[df['Controls'] == 'negative_control', lfc_col].dropna()

    X_values = np.linspace(0.0001, 0.01, 200)

    results = []

    for X in X_values:
        line_min = np.quantile(neg, X)
        line_max = np.quantile(neg, 1 - X)
        data=df.copy()
        data = data.loc[(((df[lfc_col] > line_min) & (df[lfc_col] < line_max)) | df['y_true']),:]
        ranks=data["|".join(['neg', 'rank'])]
        inverted = max(ranks) + 1 - ranks ### Ranks are highest (1) to lowest (max) Rauc want scores with highest (max) lowest (min)
        y_pred = inverted
        y_true = data['y_true']
        y_pred, y_true = np.asarray(y_pred), np.asarray(y_true)
        fpr, tpr, _ = roc_curve(y_true, y_pred)
        print(fpr)
        roc_auc = auc(fpr, tpr)

        results.append({
            'X': X,
            'line_min': line_min,
            'line_max': line_max,
            'TPR': tpr,
            'FPR': fpr,
            'AUC': roc_auc
        })

    perf_df = pd.DataFrame(results)
    print(perf_df)
    # best threshold (Youden’s J)
    perf_df['J'] = perf_df['TPR'] - perf_df['FPR']

    #best_J = perf_df.loc[perf_df['J'].idxmax()]
    best_J=None
    best_AUC = perf_df.loc[perf_df['AUC'].idxmax()]

    return perf_df, roc_auc, best_J, best_AUC




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True,
                        help='xlsx file')
    parser.add_argument('--lfc_col', default='lfc')
    parser.add_argument('-o', '--output', default='quantile_analysis')

    args = parser.parse_args()
    all_sheets = pd.read_excel(args.input, sheet_name=None)

    labels = []
    datasets=[]
    for label, data_full in all_sheets.items():
        data, roc, best_j, best_auc =evaluate_quantile_thresholds(data_full, args.lfc_col)
        #smooth_and_export(data['X'],data['J'],out_prefix=label +'ROC_J_PERF', highlight_x= best_j['X'] )
        smooth_and_export(data['X'],data['AUC'],out_prefix=label +'ROC_auc_PERF', highlight_x= best_auc['X'] )



if __name__ == "__main__":
    main()