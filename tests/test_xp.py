"""
test_xp.py — pruebas del demo del operador xp (Berry–Keating)
=============================================================

Verifican las dos afirmaciones del cierre conceptual: xp da la densidad suave exacta
(silueta) pero su espectro discreto es un peine rígido (no GUE, no los t_n).
"""

import numpy as np
import pytest

from riemann import gue, spacing, xp


def test_conteo_semiclasico_xp_es_N_suave():
    E = np.array([50.0, 200.0, 1000.0, 2500.0])
    assert np.allclose(xp.conteo_semiclasico_xp(E), spacing.N_suave(E))


def test_operador_dilatacion_hermitiano():
    H = xp.operador_dilatacion(100, 20.0)
    assert np.allclose(H, H.conj().T)


def test_operador_dilatacion_requiere_par():
    with pytest.raises(ValueError):
        xp.operador_dilatacion(101, 20.0)


def test_espectro_xp_peine_uniforme():
    L = 20.0
    niveles = xp.niveles_positivos(xp.espectro_xp(400, L), recorte=0.4)
    # paso medio ≈ 2π/L
    assert abs(np.diff(niveles).mean() - 2 * np.pi / L) < 1e-2
    # rígido: los espaciados unfolded casi no varían
    assert xp.espaciados_peine(niveles).std() < 1e-3


def test_xp_no_genera_estadistica_gue():
    niveles = xp.niveles_positivos(xp.espectro_xp(400, 20.0), recorte=0.4)
    var_xp = xp.espaciados_peine(niveles).var()
    rng = np.random.default_rng(0)
    var_gue = gue.espaciados_gue(200, rng, n_matrices=15).var()
    # el peine es rígido (var≈0); GUE tiene varianza claramente mayor
    assert var_xp < 0.01 < var_gue
