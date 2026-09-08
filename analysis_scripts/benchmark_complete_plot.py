#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 08:06:11 2026
@author: joao

Versão alternativa do gráfico comparativo de FPR: em vez de barras agrupadas,
mostra a tendência de cada preditor (CLEAN, DeepEC, DIAMOND) como uma curva
única, interpolando os pontos de Catalytic → 5Å → 10Å → 15Å.

Isso deixa mais claro o comportamento de degradação (ou robustez) de cada
método à medida que o raio de mutação aumenta, num único painel.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator

# ── Dados ────────────────────────────────────────────────────────────────────

dict_results = {
    "CLEAN":   [92.9, 93.1, 76.38, 30.81],
    "DeepEC":  [75.3, 67.8, 19.0,  4.4],
    "DIAMOND": [92.1, 92.4, 92.1,  85.7],
}
decoy_categories = ["Catalytic", "5Å", "10Å", "15Å"]
# Posição numérica de cada categoria no eixo X (distância aproximada, em Å,
# do sítio catalítico; "Catalytic" tratado como 0Å para fins de espaçamento)
x_positions = np.array([0, 5, 10, 15])

df = pd.DataFrame(dict_results, index=decoy_categories)

# ── Estilo tipográfico ───────────────────────────────────────────────────────

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Liberation Serif", "Times New Roman", "Times", "DejaVu Serif"],
    "axes.linewidth": 0.0,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "text.color": "#1a1a1a",
    "axes.edgecolor": "#1a1a1a",
})

COLORS = {
    "CLEAN":   "#7F8C8D",
    "DeepEC":  "#D55E00",
    "DIAMOND": "#0072B2",
}
MARKERS = {
    "CLEAN":   "o",
    "DeepEC":  "s",
    "DIAMOND": "^",
}

# ── Figura ───────────────────────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(12, 5.5))
fig.patch.set_facecolor("white")
ax.set_facecolor("#FAFAF8")

ax.yaxis.grid(True, linestyle="-", linewidth=0.6, color="#DDDDD8", zorder=1)
ax.set_axisbelow(True)

# Malha fina para desenhar a curva suave (interpolação monotônica, sem
# overshoot artificial entre os pontos reais — PCHIP preserva a forma)
x_smooth = np.linspace(x_positions.min(), x_positions.max(), 300)

for method in dict_results:
    y = df[method].values
    interp = PchipInterpolator(x_positions, y)
    y_smooth = interp(x_smooth)

    ax.plot(
        x_smooth, y_smooth,
        color=COLORS[method], linewidth=2.2, zorder=3, alpha=0.9,
    )
    ax.plot(
        x_positions, y,
        marker=MARKERS[method], markersize=9,
        linestyle="none",
        color=COLORS[method],
        markeredgecolor="white", markeredgewidth=1.2,
        zorder=4,
    )

    # Rótulo do método ancorado no último ponto (15Å), mais legível que legenda
    ax.annotate(
        method,
        xy=(x_positions[-1], y[-1]),
        xytext=(10, 0),
        textcoords="offset points",
        ha="left", va="center",
        fontsize=11.5, fontweight="bold",
        color=COLORS[method],
    )

    # Valores anotados em cada ponto real.
    # CLEAN e DIAMOND ficam muito próximos em Catalytic/5Å, então CLEAN
    # é anotado abaixo do ponto e DIAMOND acima, para não sobrepor.
    va_offset = (0, -12) if method == "CLEAN" else (0, 9)
    va_align  = "top" if method == "CLEAN" else "bottom"
    for xi, yi in zip(x_positions, y):
        ax.annotate(
            f"{yi:.1f}",
            xy=(xi, yi),
            xytext=va_offset,
            textcoords="offset points",
            ha="center", va=va_align,
            fontsize=8.5, fontweight="bold",
            color=COLORS[method],
        )

# ── Eixos ────────────────────────────────────────────────────────────────────

ax.set_xticks(x_positions)
ax.set_xticklabels(decoy_categories, fontsize=15, fontweight="bold")
ax.set_xlim(x_positions.min() - 1, x_positions.max() + 3.2)  # espaço p/ rótulo
ax.set_ylim(0, 108)
ax.set_ylabel("FPR (%)", fontsize=15, labelpad=8)
ax.set_xlabel("Decoy category (mutation radius)", fontsize=15, labelpad=10)
ax.tick_params(axis="y", labelsize=10.5)
ax.tick_params(axis="x", length=0)
ax.tick_params(axis="y", length=3, color="#999999")

for spine in ["top", "right", "left"]:
    ax.spines[spine].set_visible(False)
ax.spines["bottom"].set_color("#999999")
ax.spines["bottom"].set_linewidth(0.8)

plt.tight_layout()
plt.savefig("geral_benchmark.png", dpi=600, bbox_inches="tight", facecolor="white")
plt.savefig("geral_benchmark.pdf", bbox_inches="tight", facecolor="white")
plt.close()

print("Salvo: geral_benchmark.png / geral_benchmark.pdf")
