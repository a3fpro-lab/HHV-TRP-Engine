# Impact of the TRP Anisotropy Bound on Early-Universe Models

This document explains how the Hartle–Hawking–Vireon TRP engine in this
repository constrains early-universe scenarios, using only standard
slow-roll + CMB inputs and the explicit TRP ansatz implemented in
`trp_engine.py`.

We keep all assumptions and formulas explicit so that every prediction is
reproducible from the code in this repo.

---

## 1. TRP structure

The core TRP law is multiplicative:
\[
T = R \times P.
\]

In the engine, we define:

- **Reality bandwidth** (from exit-geometry entropy):
  \[
  A(N_e, H_I) = k_A\,\frac{e^{2N_e}}{H_I^2},
  \qquad
  S_{\rm geom}(N_e,H_I) = \frac{A(N_e,H_I)}{4},
  \]
  \[
  R(N_e,H_I) = \log\!\left(\frac{S_{\rm geom}(N_e,H_I)}{S_0}\right),
  \]
  where \(S_0\) is a reference entropy scale and \(k_A\) is a geometric
  prefactor.

- **Curvature complexity** (quadratic ansatz):
  \[
  C(\epsilon) = \frac{\epsilon^2}{2\sigma^2},
  \]
  where \(\epsilon\) parametrizes exit anisotropy and \(\sigma\) is a
  scale parameter.

- **Perception gain**:
  \[
  P(\epsilon) = \exp\big(-\mu\,C(\epsilon)\big).
  \]

The total TRP capacity is:
\[
T(N_e,H_I,\epsilon) = R(N_e,H_I)\,\exp\big(-\mu\,C(\epsilon)\big).
\]

The viability condition is:
\[
T(N_e,H_I,\epsilon) \;\ge\; T_{\min}.
\]

All of these are implemented directly in `trp_engine.py` as:

- `R(Ne, HI)`
- `C(eps)`
- `P(eps)`
- `T(Ne, HI, eps) = R * P`

---

## 2. Calibration: fixing \(\mu\) from CMB anisotropy

There is a one-parameter freedom in how much “difficulty” we assign to
geometry vs perception: any rescaling
\[
R \rightarrow a\,R, 
\qquad
P \rightarrow \frac{1}{a} P
\]
leaves \(T = R P\) invariant. To make the framework predictive we fix
that freedom with a **pivot calibration** at a Planck-like point:
\[
(N_e^\*, H_I^\*, \epsilon^\*).
\]

In code, the pivot is defined in `hhv_math` via:

- `cmb_params.A_s_planck_2018` (scalar amplitude),
- `cmb_params.r_upper_approx_0p005` (an upper bound on \(r\)),
- `cmb_params.N_e_pivot` (a representative pivot \(N_e^\*\)),

and we obtain \(H_I^\*\) using:
```python
HI_star = inflation.H_from_r_As(r_star, As, Mpl=1.0)
