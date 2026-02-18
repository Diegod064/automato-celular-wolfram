# 🧬 Autômatos Celulares Elementares de Wolfram

Implementação em Python dos **autômatos celulares unidimensionais elementares** de Stephen Wolfram, com foco especial na **Regra 30** — um sistema determinístico que exibe comportamento caótico e é utilizado em geração de números pseudo-aleatórios, criptografia e modelagem de padrões biológicos.

---

## 📌 Sobre o Projeto

Autômatos celulares são sistemas computacionais discretos compostos por células dispostas em uma grade. A cada passo de tempo, o estado de cada célula é atualizado segundo uma **regra local** que depende apenas do estado dos seus vizinhos imediatos.

Wolfram catalogou as **256 regras elementares** (denominadas pelo número binário que descreve seu comportamento) e as classificou em quatro classes:

| Classe | Comportamento | Exemplo |
|--------|--------------|---------|
| I      | Evolui para estado uniforme | Regra 0 |
| II     | Padrões periódicos/estáveis | Regra 250 |
| III    | Comportamento caótico/pseudo-aleatório | **Regra 30** ⭐ |
| IV     | Estrutura localizada complexa | Regra 110 |

Este projeto gera visualizações e análises quantitativas para as regras **30, 90, 110 e 250**.

---

## 📂 Estrutura do Projeto

```
automato_celular/
├── automato_celular.py     # Código-fonte principal
├── requirements.txt        # Dependências Python
├── README.md               # Este arquivo
└── resultados/             # Gerado automaticamente ao executar
    ├── comparativo_regras.png
    ├── analise_regra30.png
    ├── sensibilidade_condicoes_iniciais.png
    └── espaco_256_regras.png
```

---

## ⚙️ Como Rodar

### 1. Clone o repositório

```bash
git clone https://github.com/Diegod064/automato-celular-wolfram.git
cd automato-celular-wolfram
```

### 2. Crie e ative um ambiente virtual (opcional, mas recomendado)

```bash
python -m venv venv
source venv/bin/activate        # Linux/macOS
venv\Scripts\activate           # Windows
```

### 3. Instale as dependências

```bash
pip install -r requirements.txt
```

### 4. Execute

```bash
python automato_celular.py
```

Os resultados serão salvos na pasta `resultados/`.

---

## 📊 O que é Gerado

| Arquivo | Descrição |
|---------|-----------|
| `comparativo_regras.png` | Diagrama espaço-tempo das 4 regras canônicas lado a lado |
| `analise_regra30.png` | Análise detalhada: entropia, densidade, autocorrelação e distribuição de padrões |
| `sensibilidade_condicoes_iniciais.png` | Comparação de 3 sementes distintas para a Regra 30 |
| `espaco_256_regras.png` | Mosaico com 32 regras amostradas do espaço total de 256 |

---

## 🔬 Métricas Computadas (Regra 30)

- **Entropia de Shannon** dos padrões de 3-bits (obtido: ~2.55 bits; máximo teórico: 3.0 bits)
- **Densidade média** de células ativas por geração (~0.317)
- **Balanço** da coluna central (~0.525 de 1s — próximo do ideal 0.5)
- **Autocorrelação espacial** — mede a estrutura de longo alcance
- **Comprimento médio de runs** na coluna central (~2.15)

---

## 📐 Detalhes da Implementação

A função central `aplicar_regra(celulas, regra)` implementa a codificação binária de Wolfram:

```python
indice = 4 * vizinho_esquerdo + 2 * celula_central + vizinho_direito
novo_estado = (regra >> indice) & 1
```

Para a **Regra 30** (`00011110` em binário), as combinações que resultam em célula ativa são:
`100`, `011`, `010`, `001` — produzindo um padrão intrinsecamente assimétrico e caótico.

A grade usa **bordas periódicas** (topologia toroidal), evitando efeitos de borda artificiais.

---

## 📚 Referências

- Wolfram, S. (2002). *A New Kind of Science*. Wolfram Media.
- Princeton University COS 126 — General Computer Science (Fall 2002). *Lecture P4: Cellular Automata*.
- Wolfram, S. (1984). Universality and complexity in cellular automata. *Physica D*, 10(1–2), 1–35.
- Cook, M. (2004). Universality in elementary cellular automata. *Complex Systems*, 15(1), 1–40.

---

## 📄 Licença

MIT License — sinta-se livre para usar, modificar e distribuir.
