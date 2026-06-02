"""
test_engine.py — pruebas del motor de Collatz
=============================================

Casos conocidos y propiedades estructurales del mapa.
El caso estrella es n=27: trayectoria con máximo 9232 y tiempo total de parada 111.
"""

import numpy as np
import pytest

from collatz import engine


# ---------------------------------------------------------------------------
# Caso conocido: n = 27
# ---------------------------------------------------------------------------
def test_27_maximo_9232():
    seq = engine.trayectoria(27)
    assert max(seq) == 9232


def test_27_tiempo_parada_111():
    assert engine.tiempo_total_parada(27) == 111


def test_27_trayectoria_empieza_y_acaba_bien():
    seq = engine.trayectoria(27)
    assert seq[0] == 27
    assert seq[-1] == 1
    assert len(seq) == 112  # 111 pasos => 112 puntos


# ---------------------------------------------------------------------------
# Casos triviales
# ---------------------------------------------------------------------------
def test_uno_es_punto_fijo():
    assert engine.trayectoria(1) == [1]
    assert engine.tiempo_total_parada(1) == 0


@pytest.mark.parametrize("n,esperado", [(2, 1), (4, 2), (8, 3), (16, 4)])
def test_potencias_de_dos(n, esperado):
    # 2^k llega a 1 en exactamente k pasos (solo divisiones).
    assert engine.tiempo_total_parada(n) == esperado


def test_entrada_invalida():
    with pytest.raises(ValueError):
        engine.trayectoria(0)
    with pytest.raises(ValueError):
        engine.tiempo_total_parada(-5)


# ---------------------------------------------------------------------------
# Valuación 2-ádica
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("m,v", [(1, 0), (2, 1), (3, 0), (8, 3), (40, 3), (96, 5)])
def test_valuacion_2(m, v):
    assert engine.valuacion_2(m) == v


# ---------------------------------------------------------------------------
# Mapa acelerado
# ---------------------------------------------------------------------------
def test_paso_acelerado_27():
    # 3*27+1 = 82 = 2 * 41  => v=1, siguiente=41
    siguiente, v = engine.paso_acelerado(27)
    assert (siguiente, v) == (41, 1)


def test_paso_acelerado_devuelve_impar():
    for n in (1, 3, 5, 7, 27, 41, 999):
        siguiente, v = engine.paso_acelerado(n)
        assert siguiente % 2 == 1
        assert v >= 1  # 3n+1 con n impar siempre es par


def test_secuencia_v_termina_en_uno():
    vs = engine.secuencia_v(27)
    assert len(vs) > 0
    assert vs.dtype == np.int64
    assert (vs >= 1).all()


# ---------------------------------------------------------------------------
# Coherencia vectorizado vs escalar
# ---------------------------------------------------------------------------
def test_vectorizado_coincide_con_escalar():
    impares = np.array([1, 3, 5, 7, 27, 41, 12345], dtype=np.int64)
    sig_v, v_v = engine.paso_acelerado_vector(impares)
    for i, n in enumerate(impares):
        sig_s, v_s = engine.paso_acelerado(int(n))
        assert sig_s == sig_v[i]
        assert v_s == v_v[i]
