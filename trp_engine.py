"""
trp_engine.py

Hartle–Hawking–Vireon TRP engine.

Implements the TRP ansatz

  T(Ne, HI, eps) = R(Ne, HI) * exp( - mu * C(eps) )

with

  R(Ne, HI)   = log( S_geom(Ne,HI) / S0 )
  S_geom      = A(Ne,HI) / 4          (Planck units, 4 G hbar = 1)
  A(Ne,HI)    = k_A * exp(2 Ne) / HI**2
  C(eps)      = eps**2 / (2 sigma**2)

and a viability condition T >= T_min.

The parameter mu is calibrated at a chosen pivot (Ne_star, HI_star, eps_star)
by requiring T(Ne_star, HI_star, eps_star) = T_min.
"""

import numpy as np


class TRPEngine:
    """
    TRPEngine(S0, T_min, k_A, sigma, Ne_star, HI_star, eps_star)

    Parameters
    ----------
    S0 : float
        Reference entropy scale for the minimal recursion engine.
    T_min : float
        Minimal viable TRP capacity (dimensionless).
    k_A : float
        Geometric prefactor in the exit-area relation:
          A(Ne,HI) = k_A * exp(2 Ne) / HI**2.
    sigma : float
        Scale parameter for the curvature-complexity C(eps) = eps^2 / (2 sigma^2).
    Ne_star : float
        Pivot number of e-folds (e.g. ~60 for CMB scales).
    HI_star : float
        Pivot Hubble scale H_I at exit, in Planck units.
    eps_star : float
        Pivot anisotropy parameter, saturating the TRP bound at the pivot.
    """

    def __init__(
        self,
        S0: float = 1e5,
        T_min: float = 3.0,
        k_A: float = 1e3,
        sigma: float = 1.0,
        Ne_star: float = 60.0,
        HI_star: float | None = None,
        eps_star: float = 0.02,
    ):
        if HI_star is None:
            raise ValueError("HI_star must be provided for calibration.")

        self.S0 = float(S0)
        self.T_min = float(T_min)
        self.k_A = float(k_A)
        self.sigma = float(sigma)

        self.Ne_star = float(Ne_star)
        self.HI_star = float(HI_star)
        self.eps_star = float(eps_star)

        # Calibrate mu from the pivot condition
        self.mu = self._calibrate_mu()

    # ---- core geometric helpers ----

    def area(self, Ne, HI):
        """
        Exit-surface area A(Ne, HI) in Planck units:

          A = k_A * exp(2 Ne) / HI**2

        Accepts scalars or numpy arrays.
        """
        Ne = np.asarray(Ne, dtype=float)
        HI = np.asarray(HI, dtype=float)
        return self.k_A * np.exp(2.0 * Ne) / (HI ** 2)

    def S_geom(self, Ne, HI):
        """
        Geometric entropy S_geom = A / 4 in Planck units,
        using 4 G hbar = 1.
        """
        return self.area(Ne, HI) / 4.0

    def R(self, Ne, HI):
        """
        Reality bandwidth:

          R = log( S_geom / S0 ) = log( A / (4 S0) ).
        """
        Sg = self.S_geom(Ne, HI)
        return np.log(Sg / self.S0)

    # ---- curvature-complexity and TRP ----

    def C(self, eps):
        """
        Curvature-complexity ansatz:

          C(eps) = eps^2 / (2 sigma^2).
        """
        eps = np.asarray(eps, dtype=float)
        return (eps ** 2) / (2.0 * self.sigma ** 2)

    def T(self, Ne, HI, eps):
        """
        Total recursive processing capacity:

          T(Ne,HI,eps) = R(Ne,HI) * exp( - mu * C(eps) ).
        """
        Rval = self.R(Ne, HI)
        Cval = self.C(eps)
        return Rval * np.exp(-self.mu * Cval)

    def is_viable(self, Ne, HI, eps):
        """
        Return a boolean mask of T >= T_min.
        """
        return self.T(Ne, HI, eps) >= self.T_min

    # ---- calibration and epsilon_max ----

    def _calibrate_mu(self):
        """
        Calibrate mu from the pivot condition:

          T(Ne_star, HI_star, eps_star) = T_min.

        This implies:

          R_star * exp(-mu C_star) = T_min
          => mu = (1/C_star) * log(R_star / T_min),

        where

          R_star = R(Ne_star, HI_star),
          C_star = C(eps_star).
        """
        R_star = float(self.R(self.Ne_star, self.HI_star))
        C_star = float(self.C(self.eps_star))

        if C_star <= 0.0:
            raise ValueError("Pivot complexity C_star must be positive.")
        if R_star <= 0.0:
            raise ValueError("Pivot reality bandwidth R_star must be positive.")

        return np.log(R_star / self.T_min) / C_star

    def C_max(self, Ne, HI):
        """
        Complexity bound:

          C_max(Ne,HI) = (1/mu) * log( R(Ne,HI) / T_min ),

        defined only where R > T_min. For R <= T_min, returns 0.0
        (i.e. no allowed anisotropy).
        """
        Rval = np.asarray(self.R(Ne, HI), dtype=float)

        Cmax = np.zeros_like(Rval)
        mask = Rval > self.T_min
        Cmax[mask] = (1.0 / self.mu) * np.log(Rval[mask] / self.T_min)
        Cmax[~np.isfinite(Cmax)] = 0.0
        Cmax[Cmax < 0.0] = 0.0
        return Cmax

    def epsilon_max(self, Ne, HI):
        """
        Maximal allowed |eps| from the complexity bound:

          |eps|_max = sigma * sqrt( 2 C_max(Ne,HI) ).

        Returns zero where C_max <= 0.
        """
        Cmax = self.C_max(Ne, HI)
        Cmax = np.asarray(Cmax, dtype=float)

        eps_max = np.zeros_like(Cmax)
        mask = Cmax > 0.0
        eps_max[mask] = self.sigma * np.sqrt(2.0 * Cmax[mask])
        return eps_max
