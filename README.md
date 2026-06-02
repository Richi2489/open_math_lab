# open_math_lab

**Un laboratorio de *método*: cómo investigar una heurística matemática sin
autoengañarse.**

Este repositorio no persigue demostrar teoremas ni cazar contraejemplos. Toma señales
empíricas en objetos matemáticos clásicos y las somete a una disciplina de falsación:
hipótesis explícita, **regla de decisión pre-registrada** (escrita *antes* de ver los
resultados), corrección de sesgos de muestra finita, pruebas de permutación, y un cuidado
obsesivo por los **confounders** — esos efectos de escala que fabrican "señales" que en
realidad son artefactos. El producto típico es un **resultado negativo bien hecho y
reproducible**, que vale tanto como uno positivo cuando el método es honesto.

Cada subcarpeta-paquete es un "lab" independiente con su propia pregunta, sus scripts
reproducibles (semilla fija) y su reporte en `docs/`.

## Labs

| Lab | Pregunta | Estado |
|-----|----------|--------|
| [`collatz/`](collatz/) | ¿Hay dependencia de corto alcance (lag-1) en la dinámica impar de Collatz, más allá de la heurística i.i.d.? | **Cerrado** — artefacto de tamaño finito (ver abajo) |
| [`riemann/`](riemann/) | ¿Reproducimos limpiamente la repulsión tipo GUE (Montgomery–Odlyzko) en los gaps *unfolded* de los ceros de ζ? | Iteración 1 |

---

# Lab: Collatz

Laboratorio empírico sobre la **conjetura de Collatz**.

> **El objetivo NO es buscar contraejemplos.** La conjetura está verificada por
> supercómputo hasta 2^71. El objetivo es **medir empíricamente la heurística
> probabilística del 3/4** y, sobre todo, **cuantificar la dependencia entre pasos**:
> la heurística clásica asume que los pasos son independientes, y no lo son.

## Estado del lab Collatz

La hipótesis investigada (una autocorrelación de **lag-1 positiva** como señal de
correlación de corto alcance) se **cerró** tras una falsación con regla pre-registrada y
un des-confound magnitud-vs-longitud: el efecto decae hacia cero y es un artefacto de
tamaño finito, no señal nueva. Es un **resultado negativo bien hecho**, reproducible.

- Post-mortem (arco completo y honesto): [`docs/collatz_lag1_closure.md`](docs/collatz_lag1_closure.md)
- Decisión formal (ADR): [`docs/decisions/0001-cierre-veta-lag1.md`](docs/decisions/0001-cierre-veta-lag1.md)
- Reportes reproducibles: [`docs/iteracion2_reporte.md`](docs/iteracion2_reporte.md), [`docs/confirmatorio_reporte.md`](docs/confirmatorio_reporte.md)

## La pregunta de la iteración 1 (Collatz)

El mapa acelerado sobre los impares es

```
T(n) = (3n + 1) / 2^v,   con v = v_2(3n + 1)   (siempre v >= 1)
```

La heurística modela cada `v` como una variable geométrica independiente con
`P(v = k) = 2^-k`. De ahí salen dos predicciones limpias:

1. **E[v] ≈ 2.0** — la mitad de la heurística.
2. **Drift geométrico ≈ 0.75 = 3/4** por paso — `exp(E[log(T/n)]) = exp(log 3 − E[v]·log 2) = 3/4`. La otra mitad: en media geométrica cada paso **encoge** el número, lo que explica por qué las trayectorias tienden a 1.

La parte jugosa es la tercera:

3. **Autocorrelación de la secuencia de `v`s** de una trayectoria. Si los `v`s fueran
   independientes (como asume la heurística), la autocorrelación sería ≈ 0 para todo
   lag ≥ 1. Medir **cuánto se aleja de 0 y qué tan rápido decae** es exactamente donde
   una mirada actuarial/estocástica tiene algo que aportar.

**Criterio de "terminado" para la iteración 1:** tres números (E[v], drift, ACF) y una
gráfica de autocorrelación.

## Estructura

```
open_math_lab/
├── collatz/
│   ├── __init__.py
│   ├── engine.py     # mapa clásico y acelerado, trayectorias, tiempos de parada
│   ├── stats.py      # E[v], drift geométrico, autocorrelación, distribución de tiempos
│   └── viz.py        # gráficas (trayectorias log, histograma, autocorrelación)
├── scripts/
│   └── run_baseline.py
├── tests/
│   └── test_engine.py
├── requirements.txt
├── .gitignore
└── README.md
```

## Cómo correrlo

```powershell
# 1. Crear y activar el entorno virtual
python -m venv .venv
.\.venv\Scripts\Activate.ps1        # en Windows PowerShell

# 2. Instalar dependencias
pip install -r requirements.txt

# 3. Pruebas
pytest -q

# 4. Corrida base (imprime E[v], drift, autocorrelación; guarda gráficas en outputs/)
python scripts/run_baseline.py

# Versión rápida con menos muestras:
python scripts/run_baseline.py --muestras 200000 --trayectorias 1000
```

> En Linux/macOS la activación es `source .venv/bin/activate`.

Las gráficas se guardan en `outputs/` (ignorada por git).

## Notas

- `numba` es **opcional**; si está instalado, acelera el conteo de tiempos de parada
  sobre valores que caben en `int64`.
- Las funciones escalares usan enteros de Python de precisión arbitraria (sin riesgo de
  desbordamiento); las versiones vectorizadas con numpy usan `int64` y están pensadas
  para muestreo estadístico con valores acotados.

---

# Lab: Riemann Statistical Lab

Laboratorio estadístico reproducible de la conexión **Montgomery–Odlyzko**.

> **NO es proof-hunting de la Hipótesis de Riemann.** La pregunta de la iteración 1:
> ¿los gaps normalizados (*unfolded*) de los ceros no triviales de ζ exhiben **repulsión
> tipo GUE**, contrastados contra GUE finito simulado y contra Poisson? El confounder
> central aquí es el **unfolding** (el análogo de la magnitud en Collatz): normalizar por
> la media global, en vez de por la densidad local `N(T)`, fabrica artefactos de escala.

- Plan y **regla anti-autoengaño pre-registrada**: [`docs/riemann_lab_plan.md`](docs/riemann_lab_plan.md)
- Reporte iteración 1 (mpmath, ~10³): [`docs/spacing_replication_report.md`](docs/spacing_replication_report.md)
- Descarga de la regla (unfolding exacto + altura): [`docs/unfolding_height_check_report.md`](docs/unfolding_height_check_report.md)
- Reporte iteración 2 (Odlyzko, ~10⁶): [`docs/odlyzko_highT_report.md`](docs/odlyzko_highT_report.md)
- Cierre conceptual — el operador `xp` de Berry–Keating ("Touching the Wall"): [`docs/riemann_missing_operator.md`](docs/riemann_missing_operator.md)

**Resultado iteración 1:** repulsión tipo GUE **reproducida limpiamente** en tres métricas
(espaciados de vecino más cercano vs sorpresa de Wigner; KS/Wasserstein mucho más cerca de
GUE que de Poisson; correlación de pares R₂ con el hoyo de repulsión en `r→0`).

**Resultado iteración 2 — réplica confirmada:** subiendo la altura de ~10³ a ~10⁶ con los
2,001,052 ceros de Odlyzko (`zeros6`), la pequeña desviación ζ−GUE de baja altura **encoge**
monótonamente (|Δvar| 0.034 → 0.016; Wasserstein 0.035 → 0.016), Poisson queda descartado por
amplio margen, y la desviación decrece por ventanas de altura. El p-valor minúsculo del KS a
2M es la **trampa de muestra grande** (pre-registrada), no una anomalía: es la corrección
conocida de altura finita al GUE. La conexión Montgomery–Odlyzko queda **reproducida
limpiamente** (réplica empírica de un resultado conocido, no evidencia sobre RH).

```powershell
# Iteración 1 (N=2000 ceros vía mpmath; cachea en data/)
python scripts/run_spacing_replication.py --n-ceros 2000

# Iteración 2 (alta altura; requiere el archivo local de Odlyzko zeros6 en data/odlyzko/)
python scripts/run_odlyzko_highT.py --n 2000000 --gue-samples 40
```

Los ceros de mpmath se cachean en `data/`. El dataset de Odlyzko (`zeros6`, ~36 MB) **no se
versiona** (está en `.gitignore`); se obtiene de las
[tablas de Odlyzko](https://www-users.cse.umn.edu/~odlyzko/zeta_tables/) y se coloca en
`data/odlyzko/`. El script lee el archivo local; no descarga en runtime.
