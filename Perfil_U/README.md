# Perfil U – MFDM/MWLS (Octave/MATLAB) — v2 corrigida

## Correções em relação à v1

| Problema (v1)                                | Correção (v2)                                              |
|---------------------------------------------|------------------------------------------------------------|
| Geometria plotada como retângulo cheio       | Triangulação filtrada: triângulos no vazio são removidos   |
| `gerar_pontos_U` gerava pontos sobrepostos   | Grade regular com máscara do domínio U                     |
| Resultados diferentes do Python             | MWLS com regularização de Tikhonov e n_viz=15              |
| `tricontour` não disponível no Octave        | Removido; contorno do U desenhado com `plot`               |
| Apenas 1 gráfico aparecia no subplot         | subplot 2×2 com `colormap(ax, ...)` por eixo               |

## Arquivos

| Arquivo                  | Descrição                                        |
|--------------------------|--------------------------------------------------|
| `perfilU_main.m`         | Script principal — executa tudo                  |
| `gerar_pontos_U.m`       | Gera a nuvem de pontos do Perfil U               |
| `montar_sistema_mfdm.m`  | Monta o sistema linear MFDM/MWLS                 |

## Como executar

```octave
cd /caminho/para/perfilU_octave_v2
perfilU_main
```

## Parâmetros (em `perfilU_main.m`)

| Parâmetro | Valor | Descrição                        |
|-----------|-------|----------------------------------|
| `W`       | 4.0   | Largura total [m]                |
| `H`       | 3.0   | Altura total [m]                 |
| `th`      | 0.6   | Espessura paredes/base [m]       |
| `n_base`  | 14    | Densidade da malha               |
| `k_cond`  | 1.0   | Condutividade térmica [W/m·K]    |

## Condição de contorno

`T(x) = 500 × (x − xmin)/(xmax − xmin)`  → 0 °C na esquerda, 500 °C na direita

## Gráficos gerados

- **Figura 1** — Campos térmicos: T, qx, qy, |q| (2×2)
- **Figura 2** — Malha Delaunay com nós de contorno destacados
- **Figura 3** — Perfis T(x) em 3 seções transversais + referência analítica

---
**UFES – Iniciação Científica – MFDM para transferência de calor**
