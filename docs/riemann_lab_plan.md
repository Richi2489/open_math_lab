# Riemann Statistical Lab — plan e iteración 1

> **NO es proof-hunting de la Hipótesis de Riemann.** Es un laboratorio estadístico
> reproducible de la conexión **Montgomery–Odlyzko**: ¿los gaps normalizados (*unfolded*)
> de los ceros no triviales de ζ exhiben **repulsión tipo GUE** (Gaussian Unitary
> Ensemble)? Misma filosofía anti-autoengaño que el lab de Collatz. Lenguaje prudente.

## El paralelo con Collatz

En Collatz, el confounder fatal era la **magnitud** (mezclada con la longitud). Aquí el
confounder fatal es el **unfolding**: la densidad de ceros crece con la altura `T`, así que
normalizar dividiendo por la media **global** fabrica "desviaciones de GUE" que son puro
artefacto de densidad. La lección es la misma en terreno nuevo: si normalizas mal, te
engañas solo.

## Regla anti-autoengaño (PRE-REGISTRADA, fijada antes de ver resultados)

> **No se declarará una desviación de GUE si esa desviación desaparece al cambiar:**
> **(a) la ventana de altura `T`, (b) el método de unfolding, (c) el tamaño de muestra,**
> **o (d) el baseline.** Toda diferencia contra GUE se compara contra **GUE de tamaño
> finito simulado** (mismo `N`), **no** contra la fórmula asintótica ideal.

Corolario operativo: el "éxito" de la iteración 1 NO es encontrar una desviación
interesante; es **reproducir limpiamente la repulsión tipo GUE** y entender qué parte de
cualquier discrepancia es artefacto de tamaño finito, unfolding o ventana.

## Iteración 1 — objetivo (sanity check)

> Reproducir limpiamente la señal Montgomery–Odlyzko: que los gaps *unfolded* de los ceros
> de ζ exhiban repulsión tipo GUE, contrastado contra **GUE finito simulado** y contra
> **Poisson** — con unfolding correcto y controles de tamaño finito desde el día 1.

Métricas: distribución de espaciados de vecino más cercano (vs sorpresa de Wigner),
correlación de pares R₂ (la firma Montgomery–Odlyzko) vs `1 − (sin(πx)/(πx))²`, distancias
KS y Wasserstein (ζ vs GUE y ζ vs Poisson), momentos, y estabilidad por ventana de altura.

## Reproducibilidad

- Ceros vía `mpmath.zetazero(n)` (sin descargas externas), cacheados en `data/`.
- Hook dejado para ingerir tablas de Odlyzko de alta altura `T` en una iteración futura
  (ahí el match GUE es más limpio); **no se descarga nada todavía**.
- Semilla fija para GUE/Poisson. Parámetros (N, n de matriz, nº de matrices, ventanas) se
  registran en el reporte `docs/spacing_replication_report.md`.
