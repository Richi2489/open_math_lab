"""
test_riemann.py — pruebas del Riemann Statistical Lab
=====================================================

- Unfolding da media de espaciado ≈ 1 (round-trip por N_suave).
- GUE simulado reproduce repulsión (pocos gaps chicos, forma tipo Wigner).
- Poisson da gaps exponenciales (sin repulsión).
- Primeros ceros de ζ coinciden con valores conocidos.
"""

import numpy as np
import pytest

from riemann import gue, metrics, spacing, zeros


# ---------------------------------------------------------------------------
# Unfolding
# ---------------------------------------------------------------------------
def test_N_suave_monotona():
    T = np.linspace(20, 5000, 500)
    assert np.all(np.diff(spacing.N_suave(T)) > 0)


def test_unfolding_media_uno_roundtrip():
    # Puntos unfolded con gaps exactamente 1; al des-unfold y volver a unfold,
    # los espaciados deben recuperar media 1.
    w = np.arange(50.0, 250.0)
    gammas = spacing.N_suave_inversa(w)
    s = spacing.gaps_vecino(gammas)
    assert abs(s.mean() - 1.0) < 1e-6


def test_normalizacion_global_es_distinta():
    # El método incorrecto no coincide con el unfolding correcto sobre ceros sintéticos.
    w = np.arange(50.0, 400.0)
    gammas = spacing.N_suave_inversa(w)
    s_ok = spacing.gaps_vecino(gammas)
    s_mal = spacing.gaps_normalizacion_global_INCORRECTA(gammas)
    assert s_ok.var() != pytest.approx(s_mal.var(), rel=1e-3)


# ---------------------------------------------------------------------------
# GUE y Poisson
# ---------------------------------------------------------------------------
def test_gue_media_uno_y_repulsion():
    rng = np.random.default_rng(0)
    s = gue.espaciados_gue(120, rng, n_matrices=20, recorte=0.15)
    assert abs(s.mean() - 1.0) < 0.1
    # repulsión: muy pocos espaciados diminutos
    assert (s < 0.1).mean() < 0.05


def test_poisson_exponencial():
    rng = np.random.default_rng(1)
    s = gue.espaciados_poisson(20000, rng)
    assert abs(s.mean() - 1.0) < 0.05
    # exponencial: P(s<0.3) ≈ 1 - e^-0.3 ≈ 0.259
    assert abs((s < 0.3).mean() - (1 - np.exp(-0.3))) < 0.03


def test_gue_reprime_mas_que_poisson():
    rng = np.random.default_rng(2)
    s_gue = gue.espaciados_gue(120, rng, n_matrices=20)
    s_pois = gue.espaciados_poisson(len(s_gue), rng)
    # GUE tiene muchos menos gaps chicos que Poisson (repulsión).
    assert (s_gue < 0.3).mean() < (s_pois < 0.3).mean()


def test_wigner_y_r2_teoricos():
    # Wigner GUE se anula en s=0 y es positivo en s=1.
    assert gue.wigner_gue(0.0) == 0.0
    assert gue.wigner_gue(1.0) > 0
    # R_2 GUE: hoyo en r->0, ~1 en r=1 (sin(pi)=0).
    assert metrics.r2_teorico_gue(0.01) < 0.1
    assert abs(metrics.r2_teorico_gue(1.0) - 1.0) < 1e-6


# ---------------------------------------------------------------------------
# Ceros de ζ (computa 2; ~0.5 s, luego cacheado)
# ---------------------------------------------------------------------------
def test_primeros_ceros_conocidos():
    g = zeros.cargar_ceros(2)
    assert abs(g[0] - 14.134725) < 1e-3
    assert abs(g[1] - 21.022040) < 1e-3


def test_validar_ceros_detecta_no_monotono():
    with pytest.raises(ValueError):
        zeros.validar_ceros(np.array([14.0, 13.0, 21.0]))


# ---------------------------------------------------------------------------
# Unfolding exacto (theta de Riemann-Siegel)
# ---------------------------------------------------------------------------
def test_theta_unfold_cerca_de_asintotico():
    # A alturas moderadas, θ/π+1 y la forma asintótica coinciden a ~1e-3.
    T = np.array([100.0, 500.0, 1500.0])
    assert np.max(np.abs(spacing.theta_unfold(T) - spacing.N_suave(T))) < 1e-2


def test_unfold_exacto_gaps_media_uno():
    # Round-trip sobre ceros reales: los gaps unfolded exactos tienen media ≈ 1.
    w = zeros.cargar_unfold_exacto(60)
    s = np.diff(w)
    assert abs(s.mean() - 1.0) < 0.05
