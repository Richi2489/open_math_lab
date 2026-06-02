# open_math_lab — Collatz

Laboratorio empírico de investigación sobre la **conjetura de Collatz**.

> **El objetivo NO es buscar contraejemplos.** La conjetura está verificada por
> supercómputo hasta 2^71. El objetivo es **medir empíricamente la heurística
> probabilística del 3/4** y, sobre todo, **cuantificar la dependencia entre pasos**:
> la heurística clásica asume que los pasos son independientes, y no lo son.

## La pregunta de la iteración 1

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
