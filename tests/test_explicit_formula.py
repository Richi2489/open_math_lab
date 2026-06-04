"""
test_explicit_formula.py — pruebas de la reconstrucción de primos desde los ceros
=================================================================================
"""

import math

import numpy as np
import pytest

from riemann import explicit_formula as ef
from riemann import zeros


# ---------------------------------------------------------------------------
# von Mangoldt y ψ exacta
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("n,esperado", [
    (1, 0.0), (4, math.log(2)), (6, 0.0), (7, math.log(7)),
    (9, math.log(3)), (12, 0.0), (27, math.log(3)),
])
def test_von_mangoldt(n, esperado):
    assert abs(ef.von_mangoldt(n) - esperado) < 1e-12


@pytest.mark.parametrize("X", [10, 20, 30, 50, 97])
def test_psi_exacta_coincide_suma_directa(X):
    directo = sum(ef.von_mangoldt(n) for n in range(2, X + 1))
    assert abs(ef.psi_exacta(X) - directo) < 1e-9


def test_psi_exacta_es_escalera_con_saltos_en_primos():
    # justo antes y después de 7 hay un salto de ln(7).
    salto = ef.psi_exacta(7.0) - ef.psi_exacta(6.999)
    assert abs(salto - math.log(7)) < 1e-9


def test_pnt_psi_sobre_x_tiende_a_uno():
    assert abs(ef.psi_exacta(1e5) / 1e5 - 1.0) < 0.05


# ---------------------------------------------------------------------------
# Fórmula explícita: convergencia con N
# ---------------------------------------------------------------------------
def test_psi_aprox_converge_con_N():
    g = zeros.cargar_ceros(500)
    xg = np.linspace(2.0, 50.0, 800)
    psi_ex = ef.psi_exacta(xg, 50.0)
    err_pocos = ef.error_medio(xg, g[:20], psi_ex)
    err_muchos = ef.error_medio(xg, g[:400], psi_ex)
    assert err_muchos < err_pocos
    assert err_muchos < 0.2  # con 400 ceros la reconstrucción ya es buena


def test_psi_aprox_termino_dominante_es_x():
    # sin ceros, ψ_0(x) ≈ x − ln(2π) − ½ln(1−x⁻²): cerca de x para x grande.
    x = 80.0
    assert abs(ef.psi_aprox(x, np.empty(0)) - x) < 3.0
