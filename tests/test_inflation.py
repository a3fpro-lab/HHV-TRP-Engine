"""
tests/test_inflation.py

Tests for hhv_math.inflation helpers.

Checks:
  1. H_from_r_As and r_from_H_As are numerical inverses.
  2. slow_roll_epsilon_from_r(r) = r / 16.
  3. Larger r -> larger H_I and larger V^(1/4).
"""

import numpy as np

from hhv_math import inflation, cmb_params


def test_H_and_r_roundtrip():
    """
    r -> H_I -> r_back should recover the original r
    within tight numerical tolerance.
    """
    As = cmb_params.A_s_planck_2018
    r_vals = np.array([1e-3, 1e-2, 3e-2], dtype=float)

    H_vals = inflation.H_from_r_As(r_vals, As, Mpl=1.0)
    r_back = inflation.r_from_H_As(H_vals, As, Mpl=1.0)

    assert np.allclose(
        r_back, r_vals, rtol=1e-12, atol=0.0
    ), f"roundtrip r->H->r failed: r={r_vals}, r_back={r_back}"


def test_slow_roll_epsilon_relation():
    """
    slow_roll_epsilon_from_r(r) must equal r / 16.
    """
    r_vals = np.array([1e-4, 1e-3, 1e-2, 3e-2], dtype=float)
    eps_expected = r_vals / 16.0
    eps_actual = inflation.slow_roll_epsilon_from_r(r_vals)

    assert np.allclose(
        eps_actual, eps_expected, rtol=1e-15, atol=0.0
    ), f"epsilon mismatch: expected={eps_expected}, got={eps_actual}"


def test_monotonic_energy_scale_with_r():
    """
    Larger r -> larger H_I -> larger V^(1/4).
    """
    As = cmb_params.A_s_planck_2018
    r_small = 1e-3
    r_large = 3e-2

    H_small = inflation.H_from_r_As(r_small, As, Mpl=1.0)
    H_large = inflation.H_from_r_As(r_large, As, Mpl=1.0)

    Vq_small = inflation.V_quarter_from_H(H_small, Mpl=1.0)
    Vq_large = inflation.V_quarter_from_H(H_large, Mpl=1.0)

    assert H_large > H_small, "Expected H_large > H_small for larger r"
    assert Vq_large > Vq_small, "Expected Vq_large > Vq_small for larger r"
