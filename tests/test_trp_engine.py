"""
tests/test_trp_engine.py

Basic tests for the HHV–TRP engine.

These tests verify:

1. Pivot calibration:
   epsilon_max(Ne_star, HI_star) == eps_star (within numerical tolerance).

2. Monotonic behavior in H_I at fixed Ne:
   Lower H_I (larger area/entropy) allows larger |eps|_max.

All tests rely only on:

- hhv_math.inflation
- hhv_math.cmb_params
- trp_engine.TRPEngine
"""

import numpy as np

from hhv_math import inflation, cmb_params
from trp_engine import TRPEngine


def make_engine_from_planck():
    """
    Helper: build a TRPEngine using Planck-like pivot from hhv_math.
    """
    As = cmb_params.A_s_planck_2018
    r_star = cmb_params.r_upper_approx_0p005
    Ne_star = cmb_params.N_e_pivot

    # Compute H_I at pivot via slow-roll relation
    HI_star = inflation.H_from_r_As(r_star, As, Mpl=1.0)

    eps_star = 0.02  # from Planck-like anisotropy bound

    engine = TRPEngine(
        S0=1e5,
        T_min=3.0,
        k_A=1e3,
        sigma=1.0,
        Ne_star=Ne_star,
        HI_star=HI_star,
        eps_star=eps_star,
    )

    return engine, Ne_star, HI_star, eps_star


def test_pivot_epsilon_matches_eps_star():
    """
    At the calibration pivot, epsilon_max must match eps_star
    (within numerical tolerance), by construction:

      T(Ne*, HI*, eps*) = T_min  =>  C_max(Ne*,HI*) = C(eps*)
      => epsilon_max(Ne*,HI*) = eps*.
    """
    engine, Ne_star, HI_star, eps_star = make_engine_from_planck()

    eps_max_star = engine.epsilon_max(Ne_star, HI_star)

    assert np.isclose(
        eps_max_star, eps_star, rtol=1e-10, atol=1e-12
    ), f"pivot eps_max={eps_max_star:.6g} != eps_star={eps_star:.6g}"


def test_lower_HI_allows_larger_epsilon_max():
    """
    For fixed Ne, lowering H_I increases the area A and entropy S_geom,
    increasing R and hence C_max. This should allow larger |eps|_max.

    So for H_low < H_star < H_high, we expect:

      eps_max(H_low) > eps_max(H_high).
    """
    engine, Ne_star, HI_star, _ = make_engine_from_planck()

    H_low = 0.5 * HI_star
    H_high = 2.0 * HI_star

    eps_low = engine.epsilon_max(Ne_star, H_low)
    eps_high = engine.epsilon_max(Ne_star, H_high)

    assert eps_low > eps_high, (
        f"expected eps_max at lower H_I ({H_low:.3e}) "
        f"to exceed eps_max at higher H_I ({H_high:.3e}), "
        f"got eps_low={eps_low:.6g}, eps_high={eps_high:.6g}"
    )


def test_T_violates_for_too_large_eps():
    """
    If we choose an epsilon much larger than eps_star at the pivot,
    TRP should fall below T_min (i.e., configuration is non-viable).
    """
    engine, Ne_star, HI_star, eps_star = make_engine_from_planck()

    # Take eps = 5x the calibrated eps_star
    eps_big = 5.0 * eps_star

    T_val = engine.T(Ne_star, HI_star, eps_big)

    assert T_val < engine.T_min, (
        f"expected T < T_min for large eps={eps_big:.3g}, "
        f"got T={T_val:.6g} >= T_min={engine.T_min:.6g}"
    )
