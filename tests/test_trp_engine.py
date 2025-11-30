"""
tests/test_trp_engine.py

Basic tests for the HHV-TRP engine.

These tests verify:

1. Pivot calibration:
   epsilon_max(Ne_star, HI_star) == eps_star (within numerical tolerance).

2. Monotonic behavior in H_I at fixed Ne:
   Lower H_I (larger area/entropy) allows larger |eps|_max.

3. TRP suppression for large anisotropy:
   For eps >> eps_star at the pivot, T falls below T_min.

4. TRP factorization:
   decompose_T returns T = R * P and matches T(Ne,HI,eps).
"""

import os
import sys

# Ensure project root is on sys.path so 'hhv_math' and 'trp_engine' can be imported
ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

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


def test_TRP_factorization_decompose_T():
    """
    Check that decompose_T returns R, P, T with:

      T_dec = R * P

    and that this matches the direct T(Ne,HI,eps) evaluation.
    """
    engine, Ne_star, HI_star, eps_star = make_engine_from_planck()

    # Sample a small grid around the pivot
    Ne_list = [Ne_star - 5.0, Ne_star, Ne_star + 5.0]
    HI_factors = [0.7, 1.0, 1.3]
    eps_list = [0.0, 0.5 * eps_star, eps_star]

    for Ne in Ne_list:
        for f in HI_factors:
            HI = HI_star * f
            for eps in eps_list:
                Rval, Pval, T_dec = engine.decompose_T(Ne, HI, eps)
                T_dir = engine.T(Ne, HI, eps)

                # T must equal R * P and match direct evaluation
                assert np.isclose(
                    T_dec, Rval * Pval, rtol=1e-12, atol=0.0
                ), "decompose_T does not satisfy T = R * P"

                assert np.isclose(
                    T_dec, T_dir, rtol=1e-12, atol=0.0
                ), "decompose_T T does not match direct T(Ne,HI,eps)"
