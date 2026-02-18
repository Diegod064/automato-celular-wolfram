"""
Autômatos Celulares Unidimensionais de Wolfram
==============================================
Implementação das 256 regras elementares de Wolfram para autômatos celulares
1D, com foco nas regras 30, 90, 110 e 250. Gera visualizações e analisa
propriedades estatísticas dos padrões gerados.

Referência: Wolfram, S. (2002). A New Kind of Science.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap
from collections import Counter
import os

# ─────────────────────────────────────────────────────────────────
# NÚCLEO DO AUTÔMATO
# ─────────────────────────────────────────────────────────────────

def aplicar_regra(celulas: np.ndarray, regra: int) -> np.ndarray:
    """
    Aplica uma regra de Wolfram (0-255) a uma linha de células.

    Args:
        celulas: Array 1D de 0s e 1s representando o estado atual.
        regra:   Número da regra de Wolfram (0–255).

    Returns:
        Novo array com o estado atualizado.
    """
    n = len(celulas)
    nova = np.zeros(n, dtype=np.int8)
    for i in range(n):
        esq   = celulas[(i - 1) % n]  # vizinho esquerdo (bordas periódicas)
        centro = celulas[i]
        dir_  = celulas[(i + 1) % n]  # vizinho direito
        indice = 4 * esq + 2 * centro + dir_
        nova[i] = np.int8((regra >> int(indice)) & 1)
    return nova


def simular(regra: int, n_celulas: int = 256, n_passos: int = 128,
            semente: str = "centro") -> np.ndarray:
    """
    Simula um autômato celular 1D por n_passos iterações.

    Args:
        regra:     Número da regra de Wolfram (0–255).
        n_celulas: Número de células na grade.
        n_passos:  Número de iterações (gerações).
        semente:   "centro" (célula do meio = 1) ou "aleatoria".

    Returns:
        Matriz (n_passos × n_celulas) com o histórico de estados.
    """
    historico = np.zeros((n_passos, n_celulas), dtype=np.int8)

    if semente == "centro":
        historico[0, n_celulas // 2] = 1
    else:
        rng = np.random.default_rng(42)
        historico[0] = rng.integers(0, 2, size=n_celulas, dtype=np.int8)

    for t in range(1, n_passos):
        historico[t] = aplicar_regra(historico[t - 1], regra)

    return historico


# ─────────────────────────────────────────────────────────────────
# ANÁLISE ESTATÍSTICA
# ─────────────────────────────────────────────────────────────────

def calcular_entropia(historico: np.ndarray) -> float:
    """
    Calcula a entropia de Shannon (bits) da distribuição de padrões
    de 3 células (8 possíveis) ao longo de toda a simulação.
    """
    contagens: Counter = Counter()
    for linha in historico:
        for i in range(len(linha)):
            padrao = (int(linha[(i - 1) % len(linha)]),
                      int(linha[i]),
                      int(linha[(i + 1) % len(linha)]))
            contagens[padrao] += 1
    total = sum(contagens.values())
    entropia = 0.0
    for cnt in contagens.values():
        p = cnt / total
        if p > 0:
            entropia -= p * np.log2(p)
    return entropia


def calcular_densidade(historico: np.ndarray) -> np.ndarray:
    """Retorna a fração de células vivas (=1) em cada geração."""
    return historico.mean(axis=1)


def autocorrelacao_espacial(linha: np.ndarray) -> np.ndarray:
    """Calcula a autocorrelação normalizada de uma linha de células."""
    linha_c = linha - linha.mean()
    n = len(linha_c)
    resultado = np.correlate(linha_c, linha_c, mode='full')
    resultado = resultado[n - 1:]  # metade positiva
    if resultado[0] != 0:
        resultado = resultado / resultado[0]
    return resultado


# ─────────────────────────────────────────────────────────────────
# VISUALIZAÇÕES
# ─────────────────────────────────────────────────────────────────

CMAP_PRETO_BRANCO = ListedColormap(["white", "black"])


def plotar_regra(regra: int, n_celulas: int = 256, n_passos: int = 200,
                 semente: str = "centro", ax=None, titulo: str = None):
    """Plota o espaço-tempo de uma regra."""
    historico = simular(regra, n_celulas, n_passos, semente)
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 5))
    ax.imshow(historico, cmap=CMAP_PRETO_BRANCO, interpolation="nearest",
              aspect="auto")
    ax.set_title(titulo or f"Regra {regra}", fontsize=13, fontweight="bold")
    ax.set_xlabel("Célula", fontsize=10)
    ax.set_ylabel("Geração (tempo →)", fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    return historico


def figura_comparativa(regras=(30, 90, 110, 250), n_celulas=300,
                       n_passos=200, salvar_em="comparativo_regras.png"):
    """
    Gera figura 2×2 comparando quatro regras de Wolfram.
    Salva em PNG e retorna o caminho do arquivo.
    """
    descricoes = {30: "Regra 30 – Caos / Pseudo-aleatoriedade",
                  90: "Regra 90 – Autossimilaridade (Sierpiński)",
                  110: "Regra 110 – Estrutura Localizada (Turing-completa)",
                  250: "Regra 250 – Repetição Periódica"}

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Autômatos Celulares Elementares de Wolfram\n"
                 "Comparativo de Quatro Regras Canônicas",
                 fontsize=15, fontweight="bold", y=1.01)

    historicos = {}
    for ax, regra in zip(axes.flat, regras):
        desc = descricoes.get(regra, f"Regra {regra}")
        hist = plotar_regra(regra, n_celulas, n_passos, "centro", ax, desc)
        historicos[regra] = hist

    plt.tight_layout()
    plt.savefig(salvar_em, dpi=150, bbox_inches="tight")
    print(f"[✓] Figura comparativa salva em: {salvar_em}")
    plt.close()
    return historicos


def figura_analise_regra30(salvar_em="analise_regra30.png"):
    """
    Análise aprofundada da Regra 30: espaço-tempo, densidade,
    autocorrelação e distribuição de padrões de 3-bits.
    """
    n_celulas, n_passos = 400, 300
    hist = simular(30, n_celulas, n_passos, semente="centro")

    fig = plt.figure(figsize=(16, 12))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.35)

    # ── 1. Espaço-tempo ──
    ax1 = fig.add_subplot(gs[0, :2])
    ax1.imshow(hist, cmap=CMAP_PRETO_BRANCO, interpolation="nearest", aspect="auto")
    ax1.set_title("Diagrama Espaço-Tempo – Regra 30", fontsize=13, fontweight="bold")
    ax1.set_xlabel("Posição da célula")
    ax1.set_ylabel("Geração")

    # ── 2. Coluna central (pseudo-RNG) ──
    ax2 = fig.add_subplot(gs[0, 2])
    coluna = hist[:, n_celulas // 2]
    ax2.plot(coluna, np.arange(len(coluna)), 'k|', markersize=4)
    ax2.set_xlim(-0.2, 1.2)
    ax2.set_title("Coluna Central\n(sequência 0/1)", fontsize=11, fontweight="bold")
    ax2.set_xlabel("Estado (0 ou 1)")
    ax2.set_ylabel("Geração")
    ax2.invert_yaxis()
    frac_uns = coluna.mean()
    ax2.text(0.5, 0.02, f"Fração de 1s: {frac_uns:.3f}",
             transform=ax2.transAxes, ha="center", fontsize=9,
             bbox=dict(boxstyle="round", facecolor="lightyellow"))

    # ── 3. Densidade ao longo do tempo ──
    ax3 = fig.add_subplot(gs[1, 0])
    densidade = calcular_densidade(hist)
    ax3.plot(densidade, color="steelblue", linewidth=1)
    ax3.axhline(0.5, color="red", linestyle="--", linewidth=0.8, label="0.5")
    ax3.set_title("Densidade de Células Ativas\npor Geração", fontsize=11, fontweight="bold")
    ax3.set_xlabel("Geração")
    ax3.set_ylabel("Fração de 1s")
    ax3.set_ylim(0, 1)
    ax3.legend(fontsize=8)

    # ── 4. Autocorrelação espacial (última geração) ──
    ax4 = fig.add_subplot(gs[1, 1])
    acorr = autocorrelacao_espacial(hist[-1].astype(float))
    lags = np.arange(min(80, len(acorr)))
    ax4.plot(lags, acorr[:len(lags)], color="darkgreen", linewidth=1)
    ax4.axhline(0, color="gray", linestyle="--", linewidth=0.7)
    ax4.set_title("Autocorrelação Espacial\n(Geração Final)", fontsize=11, fontweight="bold")
    ax4.set_xlabel("Lag (células)")
    ax4.set_ylabel("Correlação normalizada")

    # ── 5. Distribuição de padrões de 3-bits ──
    ax5 = fig.add_subplot(gs[1, 2])
    contagens: Counter = Counter()
    for linha in hist:
        for i in range(len(linha)):
            p = (int(linha[(i-1) % len(linha)]),
                 int(linha[i]),
                 int(linha[(i+1) % len(linha)]))
            contagens[p] += 1
    labels = [f"{a}{b}{c}" for a, b, c in sorted(contagens)]
    valores = [contagens[k] for k in sorted(contagens)]
    cores = ["#e15759" if int("".join(map(str, k)), 2) in
             [v for v in range(8) if (30 >> v) & 1] else "#4e79a7"
             for k in sorted(contagens)]
    ax5.bar(labels, valores, color=cores)
    ax5.set_title("Frequência de Padrões de 3-bits\n(■ ativo = gera 1)", fontsize=11, fontweight="bold")
    ax5.set_xlabel("Padrão (esq-centro-dir)")
    ax5.set_ylabel("Ocorrências")

    from matplotlib.patches import Patch
    ax5.legend(handles=[Patch(color="#e15759", label="→ saída 1"),
                         Patch(color="#4e79a7", label="→ saída 0")], fontsize=8)

    entropia = calcular_entropia(hist)
    fig.text(0.5, -0.01,
             f"Entropia de Shannon dos padrões de 3-bits: {entropia:.4f} bits  "
             f"(máx teórico = 3.000 bits)",
             ha="center", fontsize=10, style="italic",
             bbox=dict(boxstyle="round", facecolor="#f0f0f0"))

    plt.savefig(salvar_em, dpi=150, bbox_inches="tight")
    print(f"[✓] Análise da Regra 30 salva em: {salvar_em}")
    plt.close()
    return hist


def figura_comparacao_sementes(salvar_em="sensibilidade_condicoes_iniciais.png"):
    """
    Demonstra a sensibilidade às condições iniciais da Regra 30
    com três sementes ligeiramente diferentes.
    """
    n_celulas, n_passos = 200, 150
    sementes_desc = [
        ("Célula central = 1",            "centro"),
        ("Células centrais 127–128 = 1",  "duas"),
        ("50 células aleatórias = 1",     "aleatoria"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(15, 6))
    fig.suptitle("Regra 30 – Sensibilidade às Condições Iniciais",
                 fontsize=13, fontweight="bold")

    for ax, (desc, tipo) in zip(axes, sementes_desc):
        if tipo == "duas":
            hist = simular(30, n_celulas, n_passos, "centro")
            # sobrescreve linha 0
            hist[0] = 0
            hist[0, n_celulas // 2]     = 1
            hist[0, n_celulas // 2 + 1] = 1
            for t in range(1, n_passos):
                hist[t] = aplicar_regra(hist[t - 1], 30)
        else:
            hist = simular(30, n_celulas, n_passos, tipo)

        ax.imshow(hist, cmap=CMAP_PRETO_BRANCO, interpolation="nearest", aspect="auto")
        ax.set_title(desc, fontsize=10, fontweight="bold")
        ax.set_xlabel("Célula")
        ax.set_ylabel("Geração" if ax == axes[0] else "")
        ax.set_xticks([])
        ax.set_yticks([])

    plt.tight_layout()
    plt.savefig(salvar_em, dpi=150, bbox_inches="tight")
    print(f"[✓] Figura de sensibilidade salva em: {salvar_em}")
    plt.close()


def figura_espaco_regras(salvar_em="espaco_256_regras.png"):
    """Visualiza um subconjunto de 32 regras em miniatura."""
    n_cel, n_pas = 101, 50
    regras_amostra = list(range(0, 256, 8))  # 32 regras

    fig, axes = plt.subplots(4, 8, figsize=(18, 9))
    fig.suptitle("Espaço das 256 Regras de Wolfram (amostra, passo = 8)",
                 fontsize=13, fontweight="bold")

    for ax, r in zip(axes.flat, regras_amostra):
        hist = simular(r, n_cel, n_pas, "centro")
        ax.imshow(hist, cmap=CMAP_PRETO_BRANCO, interpolation="nearest", aspect="auto")
        ax.set_title(f"{r}", fontsize=7)
        ax.axis("off")

    plt.tight_layout()
    plt.savefig(salvar_em, dpi=120, bbox_inches="tight")
    print(f"[✓] Espaço de regras salvo em: {salvar_em}")
    plt.close()


# ─────────────────────────────────────────────────────────────────
# PONTO DE ENTRADA
# ─────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("  AUTÔMATOS CELULARES ELEMENTARES DE WOLFRAM")
    print("=" * 60)

    os.makedirs("resultados", exist_ok=True)

    print("\n[1/4] Gerando comparativo das quatro regras canônicas...")
    historicos = figura_comparativa(
        salvar_em="resultados/comparativo_regras.png")

    print("\n[2/4] Gerando análise aprofundada da Regra 30...")
    figura_analise_regra30(salvar_em="resultados/analise_regra30.png")

    print("\n[3/4] Gerando análise de sensibilidade às condições iniciais...")
    figura_comparacao_sementes(
        salvar_em="resultados/sensibilidade_condicoes_iniciais.png")

    print("\n[4/4] Gerando visão do espaço de 256 regras...")
    figura_espaco_regras(salvar_em="resultados/espaco_256_regras.png")

    # ── Métricas da Regra 30 ──
    print("\n── Métricas da Regra 30 ─────────────────────────────────")
    hist30 = historicos[30]
    entropia = calcular_entropia(hist30)
    densidade_media = calcular_densidade(hist30).mean()
    coluna_central = hist30[:, hist30.shape[1] // 2]
    balance = coluna_central.mean()
    # Comprimentos de runs (para randomicidade)
    from itertools import groupby
    runs = [sum(1 for _ in g) for _, g in groupby(coluna_central)]
    print(f"  Entropia de Shannon (padrões 3-bits): {entropia:.4f} bits")
    print(f"  Máximo teórico: 3.0000 bits")
    print(f"  Densidade média de células ativas:    {densidade_media:.4f}")
    print(f"  Fração de 1s na coluna central:       {balance:.4f}  (ideal: 0.5)")
    print(f"  Comprimento médio de runs (coluna):   {np.mean(runs):.3f}")
    print(f"  Total de gerações simuladas:          {hist30.shape[0]}")
    print(f"  Total de células simuladas:           {hist30.size}")
    print("-" * 60)
    print("Todos os resultados foram salvos em ./resultados/")


if __name__ == "__main__":
    main()
