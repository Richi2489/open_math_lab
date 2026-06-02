"""
open_math_lab — collatz
=======================

Laboratorio empírico sobre la conjetura de Collatz.

El objetivo NO es buscar contraejemplos (verificado por supercómputo hasta 2^71),
sino medir empíricamente la heurística probabilística del 3/4 y, sobre todo,
CUANTIFICAR LA DEPENDENCIA entre pasos (la heurística asume independencia y no la hay).

Submódulos:
    engine : el mapa de Collatz (clásico y acelerado), trayectorias y tiempos de parada.
    stats  : E[v], drift geométrico, autocorrelación de v's, distribución de tiempos.
    viz    : gráficas (trayectorias log, histograma de tiempos, autocorrelación).
"""

from .engine import (
    trayectoria,
    tiempo_total_parada,
    valuacion_2,
    paso_acelerado,
    paso_acelerado_vector,
    secuencia_v,
)

__all__ = [
    "trayectoria",
    "tiempo_total_parada",
    "valuacion_2",
    "paso_acelerado",
    "paso_acelerado_vector",
    "secuencia_v",
]

__version__ = "0.1.0"
