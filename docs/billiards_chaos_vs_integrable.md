# Billares cuánticos — caos vs integrabilidad (Bohigas–Giannoni–Schmit)

_Demostración del open_math_lab. NO es teoría nueva._
_Generado por `scripts/run_billiards_pair.py` (estadio h=0.015, 700 modos, semilla 2026)._

========================================================================
  Billares cuánticos — par controlado: integrable vs caótico (BGS)
========================================================================

## Reglas pre-registradas (antes de calcular)
- **Baseline del estadio = GOE** (hay simetría de inversión temporal), NO GUE: repulsión lineal P(s)∼s, var≈0.27. Comparar contra GUE mostraría una desviación ESPERADA, no una anomalía.
- **Desimetrizar:** se trabaja en un solo sector (cuarto de estadio, Dirichlet en todo el borde, ejes de simetría incluidos). Mezclar sectores destruye la repulsión y finge Poisson — es el confounder principal del experimento.
- **Unfolding por Weyl con término de perímetro:** N(E)≈(A/4π)E−(P/4π)√E, no sólo el área.
- **Tamaño finito:** los modos bajos se desvían; el match mejora alto en el espectro (se descartan los primeros 100).

## Dominios (mismo solver de Helmholtz–Dirichlet)
- rectángulo 1.0×1.6180 (razón áurea, cuadrado irracional → evita degeneraciones aritméticas): 2000 modos analíticos; A=1.618, P=5.236.
- cuarto de estadio a=1.0, r=1.0 (FD h=0.015): 700 modos; A=1.785, P=5.571; E_max=4742.
- media de espaciado pre-normalización (Weyl): rectángulo 1.000, estadio 0.909. La leve desviación del estadio refleja la frontera escalonada del FD; se normaliza a media 1 (paso estándar de unfolding).

## Estadística de espaciados (misma maquinaria del lab de Riemann)
| dominio | n | var | P(s<0.5) | P(s<0.1) | KS vs Poisson | KS vs GOE |
|---------|---|-----|----------|----------|---------------|-----------|
| rectángulo (integrable) | 1899 | 0.9371 | 0.3939 | 0.1195 | 0.0328 | 0.2215 |
| estadio (caótico) | 599 | 0.3022 | 0.2104 | 0.0117 | 0.2019 | 0.0349 |
| — Poisson (ref.) | 7170 | 1.0292 | 0.4021 | 0.0983 | — | — |
| — GOE (ref.) | 7170 | 0.2811 | 0.1854 | 0.0075 | — | — |

Lectura: el rectángulo está más cerca de Poisson que de GOE; el estadio, al revés. P(s<0.1) ≈ 0 en el estadio (hoyo de repulsión); grande en el rectángulo.

## Figuras
- docs/figures/riemann_billiards/par_espaciados.png
- docs/figures/riemann_billiards/cdf_par.png
- docs/figures/riemann_billiards/weyl_estadio.png

## Veredicto (lenguaje prudente)
- rectángulo ≈ Poisson (no GOE): SÍ (KS vs Poisson 0.0328 < KS vs GOE 0.2215; var≈1)
- estadio ≈ GOE (no Poisson): SÍ (KS vs GOE 0.0349 < KS vs Poisson 0.2019; var≈0.27)
- la repulsión aparece SÓLO al volver caótica la mesa: SÍ (P(s<0.5): 0.3939 → 0.2104; var: 0.9371 → 0.3022)

**SE REPRODUCE BGS.** Con el MISMO solver y la MISMA maquinaria estadística, cambiar la geometría de integrable a caótica hace aparecer la repulsión de niveles tipo GOE; el caso integrable se queda en Poisson. La conexión caos → repulsión que el lab de Riemann tomó prestada queda GENERADA aquí desde cero. Prudencia: es una demostración de Bohigas–Giannoni–Schmit, no teoría nueva, y el match es de tamaño finito (mejora alto en el espectro).

========================================================================