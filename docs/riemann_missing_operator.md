# El operador que falta — por qué xp "casi" funciona (Touching the Wall)

_Demostración expositiva del Riemann Statistical Lab. NO es un intento de RH._
_Generado por `scripts/run_xp_demo.py` (grid 600, L=30.0, semilla 2026)._

========================================================================
  Riemann Lab — "Touching the Wall": el operador xp de Berry–Keating
  (cierre conceptual; NO es un intento de RH)
========================================================================

## A.1 — xp captura la SILUETA (densidad suave)
| E | xp semiclásico N(E) | Riemann–von Mangoldt suave | leading (E/2π)(logE−1) |
|---|---------------------|----------------------------|------------------------|
| 50 | 9.4228 | 9.4228 | 23.1731 |
| 200 | 79.1932 | 79.1932 | 136.8197 |
| 1000 | 648.6162 | 648.6162 | 940.2485 |
| 2500 | 1984.8086 | 1984.8086 | 2715.2016 |

- diferencia xp-semiclásico vs RvM-suave: 0.00e+00 (son la MISMA fórmula).
- Esta es exactamente la densidad que usa el unfolding del repo (`spacing.N_suave`). De ahí que xp sea el sospechoso natural.

## A.2 — el generador de dilatación discretizado: un PEINE uniforme
- caja u∈[0,30.0], grid 600 puntos; H = −i d/du (momento en u).
- primeros niveles positivos: [0.0, 0.2094, 0.4189, 0.6283, 0.8378, 1.0472]
- paso del peine ≈ 0.2094  (teórico 2π/L = 0.2094)
- desviación estándar de los espaciados unfolded: 0.00000 (≈0 ⇒ peine rígido, NO los t_n).

## A.3 — niveles y estadística: xp vs ζ vs GUE
| espectro | var espaciados | P(s<0.5) | repulsión (P(s<0.25)) |
|----------|----------------|----------|------------------------|
| xp (peine) | 0.0000 | 0.0000 | 0.0000 |
| GUE | 0.1808 | 0.1148 | 0.0166 |
| ζ (ceros) | 0.1477 | 0.0800 | 0.0080 |

- El peine xp tiene varianza de espaciados ≈ 0 (rígido), mientras GUE y los ceros comparten la distribución de Wigner con repulsión. xp es dinámicamente demasiado pobre para generar estadística GUE.

## Figuras
- docs/figures/riemann_xp/A1_silueta_conteo.png
- docs/figures/riemann_xp/A2_peine.png
- docs/figures/riemann_xp/A3_niveles_vs_zeros.png
- docs/figures/riemann_xp/A3_espaciados.png

## Veredicto A
xp reproduce la densidad media (silueta) EXACTAMENTE, pero su espectro discreto es un peine uniforme: no da los niveles individuales t_n, ni la repulsión tipo GUE, ni las correcciones aritméticas. **Toca la silueta, no el rostro.**

## Tesis del reporte
El operador buscado tendría que unificar TRES capas:
  1. la densidad suave  → la da xp (silueta);
  2. la estadística GUE → la dan las matrices aleatorias (Montgomery–Odlyzko), pero como modelo estadístico, no como un operador con esos autovalores;
  3. la huella aritmética de los primos → vive en la fórmula explícita, no en xp.
Los modelos existentes capturan PIEZAS, no el objeto. El programa de Hilbert–Pólya no está bloqueado por falta de cómputo, sino porque el objeto/geometría que tendría esas tres capas a la vez no está definido. El código ilumina el muro; no lo atraviesa. (Lenguaje prudente: esto es expositivo, no una contribución a RH.)

========================================================================