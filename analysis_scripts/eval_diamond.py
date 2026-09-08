#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Avalia o baseline DIAMOND/BLASTp contra decoys corrompidos.
Plota FPR por categoria com a mesma lógica do eval_deepec.py:
  - Prediz EC original (FP correto)
  - Prediz EC errado   (FP errado)
  - Sem predição / não-enzima (TN)

Uso:
  python eval_diamond.py \
    Results_DIAMOND/cataliticals.tsv \
    Results_DIAMOND/5A.tsv \
    Results_DIAMOND/10A.tsv \
    Results_DIAMOND/15A.tsv \
    --labels Catalytic 5Å 10Å 15Å \
    --out_dir resultados_diamond
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path


# ── Avaliação ─────────────────────────────────────────────────────────────────

def evaluate_diamond_category(csv_path: str, label: str) -> dict:
    try:
        df = pd.read_csv(csv_path, sep="\t", dtype={"distance_A": "Int64"})
    except Exception as e:
        print(f"Erro ao carregar {csv_path}: {e}")
        return None

    total = len(df)
    if total == 0:
        print(f"Aviso: {csv_path} está vazio.")
        return None

    no_pred    = (~df["has_hit"].astype(bool) | df["ec_pred"].isna()).sum()
    fp_correct = (df["has_hit"].astype(bool) & df["ec_pred"].notna() &  df["match_full"].astype(bool)).sum()
    fp_wrong   = (df["has_hit"].astype(bool) & df["ec_pred"].notna() & ~df["match_full"].astype(bool)).sum()

    print(f"  {label}: total={total:,} | "
          f"original EC={fp_correct} ({fp_correct/total*100:.1f}%) | "
          f"wrong EC={fp_wrong} ({fp_wrong/total*100:.1f}%) | "
          f"no pred={no_pred} ({no_pred/total*100:.1f}%)")

    return {
        "Category":                  label,
        "Predicts original EC (FP)": fp_correct / total * 100,
        "Predicts wrong EC (FP)":    fp_wrong   / total * 100,
        "No prediction (TN)":        no_pred    / total * 100,
    }


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="FPR DIAMOND/BLASTp por Categoria")
    parser.add_argument("tsvs",     nargs="+", help="Arquivos TSV das categorias")
    parser.add_argument("--labels", nargs="+", help="Labels das categorias")
    parser.add_argument("--out_dir", default="resultados_diamond", help="Diretório de saída")

    args = parser.parse_args()

    if args.labels and len(args.labels) != len(args.tsvs):
        print("Erro: Número de labels deve ser igual ao de TSVs.")
        exit(1)

    labels  = args.labels if args.labels else [f"Cat {i+1}" for i in range(len(args.tsvs))]
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Processando {len(args.tsvs)} categorias...")

    results = []
    for path, label in zip(args.tsvs, labels):
        r = evaluate_diamond_category(path, label)
        if r:
            results.append(r)

    summary_df = pd.DataFrame(results)
    summary_df.to_csv(out_dir / "summary_diamond.csv", index=False)

    # ── Plot ──────────────────────────────────────────────────────────────────

    sns.set_style("whitegrid")

    categories = summary_df["Category"].tolist()
    col_orig   = "Predicts original EC (FP)"
    col_wrong  = "Predicts wrong EC (FP)"
    col_nopred = "No prediction (TN)"

    vals_orig   = summary_df[col_orig].values
    vals_wrong  = summary_df[col_wrong].values
    vals_nopred = summary_df[col_nopred].values

    x     = np.arange(len(categories))
    width = 0.25

    colors = {
        col_orig:   "#C0392B",  # vermelho escuro — FP correto (pior caso)
        col_wrong:  "#E67E22",  # laranja         — FP errado
        col_nopred: "#2471A3",  # azul escuro     — TN (correto)
    }

    # Figura um pouco mais larga para abrir espaço à direita para a legenda
    fig, ax = plt.subplots(figsize=(13, 7))

    bars1 = ax.bar(x - width, vals_orig,   width, label=col_orig,   color=colors[col_orig],   edgecolor='white', linewidth=0.8)
    bars2 = ax.bar(x,          vals_wrong,  width, label=col_wrong,  color=colors[col_wrong],  edgecolor='white', linewidth=0.8)
    bars3 = ax.bar(x + width,  vals_nopred, width, label=col_nopred, color=colors[col_nopred], edgecolor='white', linewidth=0.8)

    # Anotações
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

    # Eixos
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=16)
    ax.set_xlabel("Decoy category", fontsize=16, labelpad=10)
    ax.set_ylabel("Sequences (%)", fontsize=16, labelpad=10)
    ax.tick_params(axis='y', labelsize=16)
    ax.set_ylim(0, 105)

    # Grid sutil
    ax.yaxis.grid(True, linestyle='--', alpha=0.5, color='gray')
    ax.set_axisbelow(True)
    sns.despine(left=False, bottom=False)

    # Legenda — movida para fora da área do gráfico, à direita, para não
    # sobrepor as barras
    ax.legend(
        title="Prediction outcome",
        fontsize=15,
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

    # bbox_inches='tight' garante que a legenda externa não seja cortada
    plt.savefig(out_dir / "fpr_diamond_by_category.png", dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\nResultados salvos em: {args.out_dir}")
    print(f"- fpr_diamond_by_category.png")
    print(f"- summary_diamond.csv")
