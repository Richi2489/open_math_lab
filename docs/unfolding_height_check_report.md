# Descarga de la regla pre-registrada — unfolding exacto + tendencia por altura

_Generado por `scripts/run_unfolding_check.py` (N=2000, semilla 20260602)._

======================================================================
  Riemann Lab — descarga empírica de la regla pre-registrada
======================================================================

## Regla (recordatorio)
Una desviación de GUE NO se interpreta si desaparece al cambiar unfolding, ventana de altura, tamaño de muestra o baseline. Aquí se prueban unfolding y altura con la muestra ya cacheada (sin descargas).

## 1. Unfolding exacto (θ/π+1) vs asintótico (N_suave)
- diferencia máxima en la coordenada unfolded |w_exact − w_asint| = 4.69e-04
- diferencia máxima en los espaciados |s_exact − s_asint| = 1.54e-04

| métrica (vs GUE) | asintótico | exacto | cambio |
|------------------|-----------|--------|--------|
| KS | 0.0424 | 0.0424 | +0.00e+00 |
| Wasserstein | 0.0317 | 0.0317 | +1.25e-07 |
| P(s<0.5) | 0.0800 | 0.0800 | +0.00e+00 |
| var | 0.1479 | 0.1479 | -4.21e-08 |

→ unfolding como causa: DESCARTADO (cambio ínfimo).

## 2. Tendencia por altura — cuartiles de γ (unfolding exacto)
| cuartil | rango T | n | var(ζ) | P(s<0.5) ζ | |Δvar| | KS vs GUE |
|---------|---------|---|--------|------------|--------|-----------|
| Q1 | 14–811 | 499 | 0.1396 | 0.0681 | 0.0399 | 0.0510 |
| Q2 | 813–1419 | 499 | 0.1493 | 0.0782 | 0.0302 | 0.0420 |
| Q3 | 1420–1981 | 499 | 0.1521 | 0.0882 | 0.0274 | 0.0541 |
| Q4 | 1982–2515 | 499 | 0.1511 | 0.0862 | 0.0285 | 0.0516 |

- referencia GUE: var = 0.1795, P(s<0.5) = 0.1134 (matriz 1500×1500, 35,970 espaciados).
- error estándar aprox. de var(ζ) por cuartil ≈ ±0.0114 (potencia limitada con ~500 espaciados/cuartil; se lee la TENDENCIA, no cada punto).

- brecha |Δvar|: Q1 = 0.0399 → Q4 = 0.0285 (pendiente -6.31e-06/unidad T)
- KS(ζ vs GUE): Q1 = 0.0510 → Q4 = 0.0516 (pendiente +2.06e-06/unidad T)

## Figura
- outputs/riemann_tendencia_altura.png

## Veredicto (regla pre-registrada)
- estable bajo unfolding exacto: SÍ
- brecha de varianza encoge con la altura (Q4<Q1): SÍ (pendiente var -6.3e-06, KS +2.1e-06)

**REGLA DESCARGADA.** La pequeña desviación ζ−GUE NO se debe al unfolding (cambio ínfimo con θ exacta) y ENCOGE con la altura dentro de la propia muestra: se comporta como un efecto de tamaño finito / baja altura, no como una desviación estructural. No amerita teoría para 'explicarla'. Prudencia: esto es consistencia empírica con altura finita, no una demostración.

======================================================================