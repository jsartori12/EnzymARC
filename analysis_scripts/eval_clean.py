#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Avalia a taxa de Falsos Positivos (FPR) por categoria e por Classe Principal de EC (Nível 1).

Análises incluídas:
  - FPR geral por categoria de decoy.
  - Distribuição das Classes Principais de EC (Nível 1) nos decoys.
  - Taxa de erro (FPR) específica para cada Classe Principal (Nível 1).
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# ── Funções de Parsing ──────────────────────────────────────────────────────

def parse_ec(ec_str: str) -> list[str]:
    """'ID|cataliticals|EC1; EC2' → ['EC1', 'EC2']"""
    if pd.isna(ec_str):
        return []
    parts = str(ec_str).split("|")
    if len(parts) < 3:
        return []
    return [e.strip() for e in parts[2].split(";") if e.strip()]

def parse_pred(pred_str) -> tuple:
    """'EC:1.2.3.4/0.9983' → ('1.2.3.4', 0.9983)"""
    if pd.isna(pred_str) or str(pred_str).strip() in ("", "None"):
        return None, None
    s = str(pred_str).strip()
    if s.startswith("EC:"):
        s = s[3:]
    if "/" in s:
        ec, score = s.rsplit("/", 1)
        try:
            return ec.strip(), float(score)
        except ValueError:
            return ec.strip(), None
    return s.strip(), None

def ec_prefix(ec: str, n: int):
    """'1.2.3.4', 1 → '1'"""
    if ec is None:
        return None
    parts = [p for p in ec.replace("-", "").split(".") if p]
    return ".".join(parts[:n]) if len(parts) >= n else None

def match_level(original_ecs: list[str], pred_ec, level: int) -> bool:
    """True se pred_ec bate com qualquer EC original no nível dado."""
    if pred_ec is None:
        return False
    pred_t = ec_prefix(pred_ec, level)
    if pred_t is None:
        return False
    return any(ec_prefix(r, level) == pred_t for r in original_ecs
               if ec_prefix(r, level) is not None)

def load_csv(path: str) -> pd.DataFrame:
    try:
        with open(path) as f:
            max_cols = max(len(line.split(",")) for line in f)
        df = pd.read_csv(path, sep=",", header=None,
                         names=range(max_cols), engine="python")
        return df
    except Exception as e:
        print(f"Erro ao carregar {path}: {e}")
        return pd.DataFrame()

# ── Lógica de Avaliação ──────────────────────────────────────────────────────
def evaluate_with_scores(csv_path, label, level=4):
    df = load_csv(csv_path)
    records = []
    for _, row in df.iterrows():
        original_ecs = parse_ec(row[0])
        pred_ec, score = parse_pred(row[1] if 1 in row.index else None)
        is_fp = match_level(original_ecs, pred_ec, level)
        records.append({
            "category": label,
            "score": score,
            "is_fp": is_fp,
            "result": "FP (wrong EC)" if is_fp else "TN (correct rejection)"
        })
    return pd.DataFrame(records)

def evaluate_category(csv_path: str, label: str, level: int = 4):
    """
    Avalia uma categoria e retorna estatísticas gerais e por EC Nível 1.
    """
    df = load_csv(csv_path)
    if df.empty:
        return None, None
    
    records = []
    for _, row in df.iterrows():
        original_ecs = parse_ec(row[0])
        pred_str = row[1] if 1 in row.index else None
        pred_ec, _ = parse_pred(pred_str)
        
        # O erro (FP) ainda é avaliado no nível desejado (padrão 4)
        is_fp = match_level(original_ecs, pred_ec, level)
        
        # Para análise de distribuição, pegamos apenas o Nível 1 do primeiro EC original
        primary_ec = original_ecs[0] if original_ecs else None
        ec_level1 = ec_prefix(primary_ec, 1) if primary_ec else "Unknown"
        
        records.append({
            "original_ec_lvl1": ec_level1,
            "is_fp": int(is_fp),
            "category": label
        })
        
    res_df = pd.DataFrame(records)
    
    # Estatísticas por EC Nível 1
    ec_stats = res_df.groupby("original_ec_lvl1").agg(
        count=("is_fp", "count"),
        fps=("is_fp", "sum")
    ).reset_index()
    ec_stats["fpr"] = (ec_stats["fps"] / ec_stats["count"]) * 100
    ec_stats["category"] = label
    
    return res_df, ec_stats

# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="FPR por Categoria e Classe Principal (Nível 1)")
    parser.add_argument("csvs", nargs="+", help="Arquivos CSV das categorias")
    parser.add_argument("--labels", nargs="+", help="Labels das categorias")
    parser.add_argument("--level", type=int, default=4, help="Nível EC para considerar FP (1-4, default 4)")
    parser.add_argument("--out_dir", default="resultados_nivel1", help="Diretório de saída")
    
    args = parser.parse_args()
    
    if args.labels and len(args.labels) != len(args.csvs):
        print("Erro: Número de labels deve ser igual ao de CSVs.")
        exit(1)
        
    labels = args.labels if args.labels else [f"Cat {i+1}" for i in range(len(args.csvs))]
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    all_data = []
    all_ec_stats = []
    category_summary = []
    
    print(f"Processando {len(args.csvs)} categorias (Comparação FP no nível {args.level})...")
    
    for path, label in zip(args.csvs, labels):
        raw_data, ec_stats = evaluate_category(path, label, args.level)
        if raw_data is not None:
            all_data.append(raw_data)
            all_ec_stats.append(ec_stats)
            fpr_cat = (raw_data["is_fp"].sum() / len(raw_data)) * 100
            category_summary.append({"Category": label, "FPR (%)": fpr_cat})
            print(f"  {label}: FPR = {fpr_cat:.2f}%")

    full_data = pd.concat(all_data)
    full_ec_stats = pd.concat(all_ec_stats)
    summary_df = pd.DataFrame(category_summary)
    
    # Mapeamento de nomes das classes de enzimas (Nível 1)
    ec_names = {
        "1": "1. Oxidoreductases",
        "2": "2. Transferases",
        "3": "3. Hydrolases",
        "4": "4. Lyases",
        "5": "5. Isomerases",
        "6": "6. Ligases",
        "7": "7. Translocases",
        "Unknown": "Unknown"
    }
    full_data["EC_Name"] = full_data["original_ec_lvl1"].map(lambda x: ec_names.get(x, x))
    full_ec_stats["EC_Name"] = full_ec_stats["original_ec_lvl1"].map(lambda x: ec_names.get(x, x))

    sns.set_style("whitegrid")
    
    # 1. Plot Geral por Categoria
    plt.figure(figsize=(10, 6))
    ax = sns.barplot(x="Category", y="FPR (%)", data=summary_df, hue="Category", palette="viridis", legend=True)
    #plt.legend(title="Category", bbox_to_anchor=(1.05, 1), loc='upper left')
    #plt.title(f"False Positive Rates (FPR) per Category\n(EC level {args.level})", fontsize=14)
    for p in ax.patches:
        if p.get_height() > 0:  # ignora patches fantasmas
            ax.annotate(f'{p.get_height():.1f}%', (p.get_x() + p.get_width() / 2., p.get_height()), 
                    ha='center', va='center', xytext=(0, 9), textcoords='offset points', fontweight='bold', fontsize=14)
    ax.tick_params(axis = 'y', labelsize = 16)
    ax.tick_params(axis = 'x', labelsize = 16)
    ax.set_xlabel('Category', fontsize=16)
    ax.set_ylabel('FPR(%)', fontsize=16)
    plt.savefig(out_dir / "fpr_by_category.png", dpi=300)
    plt.close()

    # 2. Distribuição das Classes Principais (Nível 1)
    plt.figure(figsize=(12, 6))
    dist = full_data["EC_Name"].value_counts().sort_index()
    sns.barplot(x=dist.values, y=dist.index, hue=dist.index, palette="viridis", legend=False)
    #plt.title("EC classes distribution in decoys", fontsize=14)
    plt.xlabel("Number of sequences")
    plt.ylabel("EC Class")
    plt.savefig(out_dir / "ec_level1_distribution.png", dpi=300)
    plt.close()

    # 3. Taxa de Erro (FPR) por Classe Principal (Nível 1)
    plt.figure(figsize=(10, 5))
    
    # FIX: Assigning the barplot output to 'ax' so tick_params can find it
    ax = sns.barplot(
        x="fpr", 
        y="EC_Name", 
        hue="category", 
        data=full_ec_stats.sort_values("EC_Name"), 
        palette="magma"
    )
    
    # Customizing labels
    plt.xlabel("False positives (%)", fontsize=14)
    plt.ylabel("Class", fontsize=14)
    
    # Customizing legend
    plt.legend(title="Category", bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=11, title_fontsize=12)
    
    # Modifying tick parameters through the correctly assigned ax object
    ax.tick_params(axis='y', labelsize=12)
    ax.tick_params(axis='x', labelsize=12)
    
    # CRITICAL: tight_layout must be called AFTER all labels, legends, and ticks are defined
    plt.tight_layout()
    
    # Saving with tight bounding box to prevent clipping the legend
    plt.savefig(out_dir / "fpr_by_ec_level1.png", dpi=300, bbox_inches='tight')
    plt.close()
        
        # Violin plot: distribuição de scores por resultado e categoria
    plt.figure(figsize=(12, 6))
    sns.violinplot(x="category", y="score", hue="result", data=all_scores_df,
                   split=True, palette={"FP (wrong EC)": "salmon", "TN (correct rejection)": "steelblue"},
                   inner="quartile")
    plt.xlabel("Category", fontsize=14)
    plt.ylabel("Prediction confidence score", fontsize=14)
    plt.legend(title="Prediction result", fontsize=12)
    plt.tight_layout()
    # Salvar CSVs
    summary_df.to_csv(out_dir / "summary_categories.csv", index=False)
    full_ec_stats.to_csv(out_dir / "detailed_fpr_by_ec_level1.csv", index=False)
    
    print(f"\nResultados salvos no diretório: {args.out_dir}")
    print(f"- fpr_by_category.png: Comparação entre categorias")
    print(f"- ec_level1_distribution.png: Distribuição das classes principais")
    print(f"- fpr_by_ec_level1.png: Erros por classe principal")
    print(f"- detailed_fpr_by_ec_level1.csv: Tabela completa de dados")

