"""
examples/pivot_from_planck.py

Demonstration script:

- Uses Planck-like CMB parameters from hhv_math.cmb_params
- Uses slow-roll relations from hhv_math.inflation
- Calibrates the TRPEngine at the pivot (Ne*, H_I*, eps*)
- Prints derived quantities:
    * H_I in Planck units and GeV
    * V^(1/4) energy scale in GeV
    * max allowed |eps| at the pivot and nearby points
"""

import numpy as np

from hhv_math import inflation, cmb_params, units
from trp_engine import TRPEngine


def main():
    # --- 1. Load Planck-like parameters ---

    As = cmb_params.A_s_planck_2018
    r_star = cmb_params.r_upper_approx_0p005
    Ne_star = cmb_params.N_e_pivot

    print("Planck-like pivot parameters:")
    print(f"  A_s      = {As:.3e}")
    print(f"  r_star   = {r_star:.3e}")
    print(f"  Ne_star  = {Ne_star:.1f}")

    # --- 2. Compute H_I at pivot from slow-roll relation ---

    HI_star = inflation.H_from_r_As(r_star, As, Mpl=1.0)
    print(f"\nDerived H_I from r and A_s:")
    print(f"  H_I* (Planck units) = {HI_star:.3e}")
    print(f"  H_I* (GeV)          = {units.H_in_GeV(HI_star):.3e} GeV")

    # --- 3. Estimate inflationary energy scale V^(1/4) ---

    V_quarter = inflation.V_quarter_from_H(HI_star, Mpl=1.0)
    V_quarter_GeV = units.V_quarter_in_GeV(V_quarter)
    print(f"\nApproximate inflation energy scale:")
    print(f"  V^(1/4) (Planck units) = {V_quarter:.3e}")
    print(f"  V^(1/4) (GeV)          = {V_quarter_GeV:.3e} GeV")

    # --- 4. Instantiate TRPEngine calibrated at this pivot ---

    eps_star = 0.02  # from Planck-like statistical isotropy bound

    engine = TRPEngine(
        S0=1e5,
        T_min=3.0,
        k_A=1e3,
        sigma=1.0,
        Ne_star=Ne_star,
        HI_star=HI_star,
        eps_star=eps_star,
    )

    print(f"\nTRP engine calibrated at pivot:")
    print(f"  mu        = {engine.mu:.3e}")
    print(f"  S0        = {engine.S0:.3e}")
    print(f"  T_min     = {engine.T_min:.3e}")
    print(f"  k_A       = {engine.k_A:.3e}")
    print(f"  eps_star  = {eps_star:.3e}")

    # --- 5. Check |eps|_max at the pivot ---

    eps_max_star = engine.epsilon_max(Ne_star, HI_star)
    print(f"\nMax allowed |eps| at pivot:")
    print(f"  |eps|_max(Ne*={Ne_star:.1f}, H_I*={HI_star:.3e}) = {eps_max_star:.5f}")

    # --- 6. Sample a small grid around the pivot ---

    print("\nSample |eps|_max around pivot (Ne varying, H_I fixed):")
    for Ne in [50.0, 55.0, 60.0, 65.0]:
        eps_max = engine.epsilon_max(Ne, HI_star)
        print(f"  Ne={Ne:5.1f} -> |eps|_max = {eps_max:.5f}")

    print("\nSample |eps|_max around pivot (H_I varying, Ne fixed):")
    for HI in [0.5 * HI_star, HI_star, 2.0 * HI_star]:
        eps_max = engine.epsilon_max(Ne_star, HI)
        print(f"  H_I={HI:.3e} -> |eps|_max = {eps_max:.5f}")


if __name__ == "__main__":
    main()
