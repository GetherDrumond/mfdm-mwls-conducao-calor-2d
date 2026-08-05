# mfdm-mwls-conducao-calor-2d

Método **MFDM/MWLS** (*Meshless Finite Difference Method* / *Moving Weighted Least Squares*) aplicado à condução de calor 2D em domínios irregulares.

Repositório de código, geometrias e malhas do Relatório Final de Iniciação Científica — **UFES**, Programa Institucional de Iniciação Científica, Edital PIIC 2025-2/2026-1.

- **Autor:** Gether Filipe Pita Drumond
- **Orientador:** Prof. Juan Sérgio Romero Saenz
- **Método de referência:** Milewski, *Numerical Algorithms* 63 (2013), adaptado para condução de calor

---

## Estrutura

O repositório é organizado por geometria estudada, na mesma ordem em que os casos aparecem no relatório:

```
Perfil_I/              perfil retangular (trilho ferroviário)
Perfil_U/              perfil em canal
Aleta_Trapezoidal/     aleta trapezoidal com dois furos circulares
├── geometria/         aleta.geo (Gmsh, kernel OpenCASCADE)
├── malhas/            aleta.msh, aleta_M*.msh, fino_r*.msh
├── octave/            sete funções do pipeline MFDM
└── python/            implementação NumPy/SciPy + solver FEM P1
```

## Correspondência com o relatório

| Caso | Pasta | Resultados |
|---|---|---|
| Perfil I — trilho | `Perfil_I/` | campos de temperatura e fluxo; verificação contra solução analítica |
| Perfil U — canal | `Perfil_U/` | idem, com correção da visualização em domínio não convexo |
| Aleta trapezoidal | `Aleta_Trapezoidal/` | verificação, campo físico, convergência, fator de concentração e comparação MFDM × FEM |

## As duas implementações

O trabalho emprega duas implementações do mesmo método, com papéis complementares:

- **OCTAVE** — código original analisado e corrigido ao longo da IC. Resolve **uma discretização por execução**; produz a verificação e os campos de cada geometria.
- **Python 3 (NumPy/SciPy)** — mesma formulação, mesmos parâmetros (`m=16`, `p=2`, `g=0,1`, `λ=1e-12`) e as mesmas correções. Desenvolvida para os estudos que exigem **varredura sobre múltiplas malhas**: convergência, análise paramétrica do raio e comparação com o FEM P1.

Ambas reproduzem o campo linear em precisão de máquina sobre a mesma malha (erro L2 da ordem de 1e-10).

---

## Como executar — Aleta Trapezoidal

### 1. Gerar a malha (fora do Octave)

```bash
gmsh Aleta_Trapezoidal/geometria/aleta.geo -2
```

Variações úteis:

```bash
gmsh aleta.geo -2 -clscale 0.75          # malha mais fina
gmsh -setnumber r 0.30 aleta.geo -2      # outro raio de furo
```

> Use `-clscale` em vez de alterar o `lc`: ele escala também o campo de refino junto aos furos, preservando a proporção. Mexer só no `lc` deixaria o refino do furo fixo e contaminaria o estudo de convergência.

### 2. Executar no Octave

Copie o `aleta.msh` para a pasta `octave/`, abra o Octave nela e rode:

```octave
TEST_HEAT_aleta
```

O script executa dois casos em sequência:

- **Caso A — verificação.** Impõe Dirichlet linear em todo o contorno, inclusive nas bordas dos furos. Como funções lineares são harmônicas, a solução exata é conhecida e os furos não perturbam o campo: qualquer desvio é erro de implementação.
- **Caso B — caso físico.** Contorno externo com Dirichlet e furos isolados (Neumann homogêneo).

### 3. Conferir o teste de aceitação

| Grandeza | Esperado | Tolerância |
|---|---|---|
| Tmin / Tmax | 0,00 / 500,00 °C | exato |
| qx médio | −83,333 W/m² | erro < 0,01% |
| erro L2 em T | ~1e-10 | < 1e-6 |

O script imprime o veredito ao final:

```
>>> TESTE DE ACEITACAO: APROVADO
```

Se sair `REPROVADO`, o problema está na implementação e não adianta interpretar o Caso B. Os dois pontos a conferir primeiro são o `diag(esc)*Mn` no final do `mwls.m` (reescalonamento das derivadas) e a indexação dos nós no `msh2data.m`.

### 4. Estudos com múltiplas malhas

Convergência, fator de concentração e comparação com o FEM exigem laços sobre várias discretizações e são executados pela implementação em Python:

```bash
python Aleta_Trapezoidal/python/rodar_tudo.py
```

---

## O pipeline em Octave

Os sete arquivos `.m` não são scripts independentes: formam uma cadeia em que a saída de um é a entrada do seguinte.

```
TEST_HEAT_aleta.m
 ├── msh2data.m                le a malha, atribui os codigos de contorno
 ├── normais_contorno.m        normais externas nodais
 ├── localMFDM_heat_aleta.m    monta e resolve A*T = B
 │    ├── star_visivel.m       estrela de vizinhos de cada no
 │    └── mwls.m               operadores locais
 └── fluxos_aleta.m            pos-processamento (Lei de Fourier)
      ├── star_visivel.m       (mesma estrela)
      └── mwls.m               (mesmo operador)
```

`star_visivel` e `mwls` são chamadas **duas vezes por nó** — uma na montagem do sistema, outra no pós-processamento. É por isso que dominam o tempo de execução.

### A matriz M

Quase tudo gira em torno da matriz devolvida por `mwls.m`. Ela tem 6 linhas, cada uma aproximando um operador diferencial a partir dos `m` valores nodais da estrela:

| Linha | Operador | Uso |
|---|---|---|
| `M(1,:)` | T (valor) | não usada |
| `M(2,:)` | ∂T/∂x | fluxo qx e Neumann |
| `M(3,:)` | ∂T/∂y | fluxo qy e Neumann |
| `M(4,:)` | ∂²T/∂x² | Laplaciano |
| `M(5,:)` | ∂²T/∂x∂y | não usada |
| `M(6,:)` | ∂²T/∂y² | Laplaciano |

Multiplicar uma linha de `M` pelo vetor de temperaturas dos vizinhos (`M(2,:)*T(idx)`) devolve a derivada naquele nó. É esse mecanismo que substitui a malha estruturada do método de diferenças finitas clássico.

### Códigos de contorno

O vetor `BC` recebe o número da *Physical Curve* definida no `.geo`:

| Código | Fronteira | Condição no Caso B |
|---|---|---|
| 0 | nó interno | resolve o Laplaciano |
| 1 | topo | Dirichlet |
| 2 | base | Dirichlet |
| 3 | lados inclinados | Dirichlet |
| 4 | borda dos furos | Neumann (isolado) |

Nos vértices em que duas fronteiras se encontram vale o **menor código**, de modo que Dirichlet prevalece sobre Neumann. Isso evita um sistema mal-posto: se um canto ficasse só com Neumann, a temperatura ali não teria referência.

---

## As três correções do domínio com furos

Nenhuma delas se manifestava nos Perfis I e U. Todas apareceram ao passar para o domínio multiplamente conexo.

| # | Arquivo | Sintoma sem a correção | O que foi feito |
|---|---|---|---|
| 5 | `mwls.m` | O erro **cresce** sob refino (2,0e-4 → 8,4e-3). Falha silenciosa. | Adimensionalizar por `dm = max(d)` e reescalar as derivadas por `dm⁻¹` e `dm⁻²` |
| 6 | `normais_contorno.m` | K = 50 em vez de 2; normal do furo invertida | Orientar pelo triângulo adjacente, no sentido oposto ao vértice oposto |
| 7 | `star_visivel.m` | K = 5,58 para r = 0,15 m; a estrela atravessava o furo | Rejeitar o vizinho cujo segmento até o nó central cruze o interior de um furo |

A **Correção 5** é a mais importante das três porque não produz erro visível: o código roda, gera gráficos plausíveis, e só se revela ao comparar malhas de refinos diferentes.

---

## Parâmetros numéricos

| Parâmetro | Valor | Controla | Efeito ao alterar |
|---|---|---|---|
| `m` | 16 | vizinhos na estrela | menor → operador mais local, mais sensível à distribuição dos nós; maior → mais suave e mais caro |
| `p` | 2 | ordem da base polinomial | `p=2` dá as derivadas segundas do Laplaciano; `p=1` só daria as primeiras |
| `g` | 0,1 | suavização dos pesos | `g=0` dá peso singular (1/d⁶) e pode gerar matriz mal-condicionada |
| `lambda` | 1e-12 | regularização de Tikhonov | estabiliza a inversão local |
| `kcond` | 1,0 | condutividade térmica | escala linearmente os fluxos |
| `qvol` | 0,0 | geração volumétrica | `0` → Laplace; `≠0` → Poisson |

---

## Requisitos

- **GNU Octave** 9.x (ou MATLAB)
- **Python** 3.x com NumPy e SciPy
- **Gmsh** 4.x (apenas para regerar as malhas; os `.msh` já estão no repositório)

## Licença

MIT — ver [LICENSE](LICENSE).
