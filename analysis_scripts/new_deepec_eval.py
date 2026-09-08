#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Avalia a taxa de Falsos Positivos (FPR) por categoria de decoy para o DeepEC,
usando CSVs resumo já agregados (um por categoria), considerando apenas o
EC Level 4 (nível mais rigoroso — corresponde ao EC completo).

Cada CSV de entrada deve ter as colunas:
  Top-K, EC Level, FPR (%), N FP, N TN, N Com Pred, N Sem Pred, N Total

Para a linha com EC Level == 4, calcula três categorias:
  - Predicts original EC (FP)  = FPR (%) da linha
  - Predicts wrong EC (FP)     = (N Com Pred / N Total * 100) - FPR (%)
  - Predicts non-enzyme (TN)   = N Sem Pred / N Total * 100

Uso:
  python eval_deepec_by_category.py \
    Results_DeepEC/Cataliticals/deepec_decoy_summary.csv \
    Results_DeepEC/5A/deepec_decoy_summary.csv \
    Results_DeepEC/10A/deepec_decoy_summary.csv \
    Results_DeepEC/15A/deepec_decoy_summary.csv \
    --labels Catalytic 5Å 10Å 15Å \
    --top_k Top-1 \
    --out_dir resultados_deepec
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

EC_LEVEL = 4


# ── Lógica de Avaliação ──────────────────────────────────────────────────────

def evaluate_deepec_category(csv_path: str, label: str, top_k: str = None) -> dict:
    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"Erro ao carregar {csv_path}: {e}")
        return None

    if top_k:
        df = df[df["Top-K"] == top_k]

    row = df[df["EC Level"] == EC_LEVEL]
    if row.empty:
        print(f"Aviso: {csv_path} não tem linha para EC Level {EC_LEVEL}.")
        return None
    row = row.iloc[0]

    total      = row["N Total"]
    fpr        = row["FPR (%)"]
    n_com_pred = row["N Com Pred"]
    n_sem_pred = row["N Sem Pred"]

    fp_correct = fpr
    fp_wrong   = (n_com_pred / total * 100) - fpr
    tn         = n_sem_pred / total * 100

    print(f"  {label}: original EC={fp_correct:.1f}% | "
          f"wrong EC={fp_wrong:.1f}% | non-enzyme={tn:.1f}%")

    return {
        "Category":                    label,
        "Predicts original EC (FP)":   fp_correct,
        "Predicts wrong EC (FP)":      fp_wrong,
        "Predicts non-enzyme (TN)":    tn,
    }


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="FPR DeepEC por Categoria (EC Level 4, a partir de CSVs resumo)")
    parser.add_argument("csvs",    nargs="+", help="CSVs resumo das categorias")
    parser.add_argument("--labels", nargs="+", help="Labels das categorias")
    parser.add_argument("--top_k", default=None, help="Filtrar por um valor de Top-K (ex: Top-1). Se omitido, usa a(s) linha(s) de EC Level 4 direto.")
    parser.add_argument("--out_dir", default="resultados_deepec", help="Diretório de saída")

    args = parser.parse_args()

    if args.labels and len(args.labels) != len(args.csvs):
        print("Erro: Número de labels deve ser igual ao de CSVs.")
        exit(1)

    labels  = args.labels if args.labels else [f"Cat {i+1}" for i in range(len(args.csvs))]
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Processando {len(args.csvs)} categorias (EC Level {EC_LEVEL})...")

    results = []
    for path, label in zip(args.csvs, labels):
        r = evaluate_deepec_category(path, label, args.top_k)
        if r:
            results.append(r)

    summary_df = pd.DataFrame(results)
    summary_df.to_csv(out_dir / "summary_deepec.csv", index=False)

    # ── Plot ─────────────────────────────────────────────────────────────────

    sns.set_style("whitegrid")

    categories  = summary_df["Category"].tolist()
    col_orig    = "Predicts original EC (FP)"
    col_wrong   = "Predicts wrong EC (FP)"
    col_nonenz  = "Predicts non-enzyme (TN)"

    vals_orig   = summary_df[col_orig].values
    vals_wrong  = summary_df[col_wrong].values
    vals_nonenz = summary_df[col_nonenz].values

    x      = np.arange(len(categories))
    width  = 0.25

    colors = {
        col_orig:   "#C0392B",   # vermelho escuro — FP correto (pior caso)
        col_wrong:  "#E67E22",   # laranja         — FP errado
        col_nonenz: "#2471A3",   # azul escuro     — TN (correto)
    }

    fig, ax = plt.subplots(figsize=(13, 7))

    bars1 = ax.bar(x - width, vals_orig,   width, label=col_orig,   color=colors[col_orig],   edgecolor='white', linewidth=0.8)
    bars2 = ax.bar(x,          vals_wrong,  width, label=col_wrong,  color=colors[col_wrong],  edgecolor='white', linewidth=0.8)
    bars3 = ax.bar(x + width,  vals_nonenz, width, label=col_nonenz, color=colors[col_nonenz], edgecolor='white', linewidth=0.8)

    for bars in [bars1, bars2, bars3]:
        for bar in bars:
            h = bar.get_height()
            if h > 0.5:
                ax.annotate(
                    f'{h:.1f}%',
                    xy=(bar.get_x() + bar.get_width() / 2, h),
                    xytext=(0, 6),
                    textcoords="offset points",
                    ha='center', va='bottom',
                    fontsize=14, fontweight='bold', color='#2C3E50'
                )

    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=14)
    ax.set_xlabel("Decoy category", fontsize=15, labelpad=10)
    ax.set_ylabel("Sequences (%)", fontsize=15, labelpad=10)
    ax.tick_params(axis='y', labelsize=13)
    ax.set_ylim(0, 105)

    ax.yaxis.grid(True, linestyle='--', alpha=0.5, color='gray')
    ax.set_axisbelow(True)
    sns.despine(left=False, bottom=False)

    ax.legend(
        title="Prediction outcome",
        fontsize=13,
        title_fontsize=14,
        frameon=True,
        framealpha=0.9,
        edgecolor='gray',
        loc='upper left',
        bbox_to_anchor=(1.01, 1.0),
        handlelength=1.2,
        handleheight=0.8,
        borderpad=0.6,
        labelspacing=0.6,
    )

    plt.savefig(out_dir / "fpr_deepec_by_category.png", dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\nResultados salvos em: {args.out_dir}")
    print(f"- fpr_deepec_by_category.png")
    print(f"- summary_deepec.csv")
