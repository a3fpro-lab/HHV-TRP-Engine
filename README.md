![CI](https://github.com/a3fpro-lab/HHV-TRP-Engine/actions/workflows/ci.yml/badge.svg)

# Hartle–Hawking–Vireon TRP Engine

Planck-calibrated **Total Recursive Processing (TRP)** selection factor
for the Hartle–Hawking no-boundary wavefunction.

This repo provides:

- A small, explicit Python class `TRPEngine` implementing
  \[
  T(N_e,H_I,\epsilon) = R(N_e,H_I)\,\exp\big(-\mu\,C(\epsilon)\big),
  \]
  with
  \[
  R = \log\!\frac{S_{\rm geom}}{S_0},\quad
  S_{\rm geom} = \frac{A(N_e,H_I)}{4},\quad
  A = k_A\,\frac{e^{2N_e}}{H_I^2},\quad
  C(\epsilon) = \frac{\epsilon^2}{2\sigma^2}.
  \]
- Calibration of \(\mu\) at a Planck-like pivot
  \((N_e^\*, H_I^\*, \epsilon^\*)\) consistent with CMB constraints.
- Utilities in `hhv_math/` for standard slow-roll relations
  \((H_I, r, A_s, V^{1/4})\) and Planck 2018–style parameter values.
- Examples (`examples/`) that:
  - derive \(H_I\) and \(V^{1/4}\) from \(A_s\) and \(r\),
  - compute \(|\epsilon|_{\max}(N_e,H_I)\),
  - export a CSV grid over the \((N_e,H_I)\) plane.
- Tests (`pytest`) + GitHub Actions CI to guarantee reproducibility.

The goal is not new microphysics, but a clean, test-backed bridge:
given any anisotropic Hartle–Hawking calculation, you can plug its
geometry into this TRP engine and obtain a transparent bound on the
maximal allowed anisotropy \(|\epsilon|_{\max}(N_e,H_I)\).


# HHV–TRP Engine

Hartle–Hawking–Vireon (HHV) TRP engine: a derivation-based framework that adds a **Total Recursive Processing (TRP)** selection factor on top of standard Hartle–Hawking–style exit geometries, calibrated to **Planck + BICEP/Keck** constraints on inflation and statistical isotropy.

This repo is deliberately restricted to:
- **Known, searchable physics** (standard GR, inflation, CMB constraints), and
- A **clearly stated ansatz** for the TRP factor, with all formulas explicit.

No pseudoscience, no metaphors. If a formula is here, it is either:
- Standard textbook / review-article physics, or
- An explicit, labelled model built from those ingredients.

---

## 1. Background: inflation scale and isotropy (data-backed)

### 1.1 Inflation scale from tensors and scalar amplitude

Planck 2018 + BICEP/Keck give upper bounds on the tensor-to-scalar ratio at CMB scales. A representative combined limit is (pivot \(k \sim 0.005~\mathrm{Mpc}^{-1}\)):

\[
r \lesssim 0.03 \quad (95\%~\text{C.L.})
\]

with scalar amplitude

\[
A_s \simeq 2.1 \times 10^{-9}.
\]

For single-field slow-roll, the tensor power and scalar power are

\[
P_t = \frac{2 H_I^2}{\pi^2 M_{\rm Pl}^2}, \qquad
P_s = A_s, \qquad
r = \frac{P_t}{P_s} = \frac{2 H_I^2}{\pi^2 A_s M_{\rm Pl}^2},
\]

so

\[
H_I = \frac{\pi}{\sqrt{2}} \sqrt{r A_s}\, M_{\rm Pl}.
\]

Taking a **conservative high-scale pivot**

\[
r_\star = 0.03, \qquad A_s = 2.1\times 10^{-9},
\]

gives

\[
H_{I,\star} \approx 1.8 \times 10^{-5}\,M_{\rm Pl}.
\]

We choose a standard CMB-scale exit of

\[
N_{e,\star} = 60
\]

e-folds between horizon exit and the end of inflation.

### 1.2 Statistical isotropy and quadrupole anisotropy

Planck tests statistical isotropy by constraining quadrupolar power anisotropy parameters (often denoted \(g_\*\)) and Bianchi models. For axisymmetric quadrupolar modulation of the scalar power spectrum, Planck finds (68% C.L. range, depending on exact model) roughly

\[
|g_\*| \lesssim \mathcal{O}(10^{-2}).
\]

We will use

\[
|\epsilon_\star| \approx 0.02
\]

as a **conservative upper bound** on a dimensionless exit anisotropy parameter \(\epsilon\) at our pivot point \((N_{e,\star}, H_{I,\star})\). We do **not** claim \(\epsilon\) is exactly equal to a specific \(g_\*\); we only require it be of comparable magnitude for a toy one-parameter description.

---

## 2. Geometry, entropy, and TRP ansatz

Everything here is defined explicitly; the model assumptions are separated from known physics.

### 2.1 Exit surface area and geometric entropy

We parameterize a quasi-de Sitter inflationary patch whose comoving region (our observable universe) exits the Hubble radius at \(H_I\) and undergoes \(N_e\) additional e-folds.

We define an **exit surface area** scaling as

\[
A(N_e,H_I) = k_A\,\frac{e^{2N_e}}{H_I^2},
\]

where \(k_A\) is a dimensionful prefactor capturing the detailed definition of the exit 3-surface and reheating history. For numerical work, we treat \(k_A\) as \(\mathcal{O}(10^2{-}10^3)\), which only affects overall normalization.

In Planck units with \(4G\hbar = 1\), the **geometric entropy** of the exit surface is

\[
S_{\text{geom}}(N_e,H_I)
= \frac{A(N_e,H_I)}{4}.
\]

This is a standard Bekenstein–Hawking-style area–entropy relation, used here heuristically for the quasi-local exit surface.

### 2.2 Reality bandwidth \(R(N_e,H_I)\)

We define a **reality bandwidth** \(R\) as the logarithm of geometric entropy measured against a reference scale \(S_0\):

\[
R(N_e,H_I) = \log\!\left(\frac{S_{\text{geom}}(N_e,H_I)}{S_0}\right)
           = \log\!\left(\frac{A(N_e,H_I)}{4S_0}\right).
\]

Here:

- \(S_0\) is interpreted as an entropy scale associated with the **minimal viable recursion engine** (smallest system capable of sustained self-measurement).
- Typical numerical choice: \(S_0 = 10^5\), but this is transparent and adjustable.

### 2.3 Curvature complexity and perception gain

We consider a single dimensionless anisotropy parameter \(\epsilon\) describing a global geometric distortion of the exit surface (e.g. squashing of \(S^3\)).

We define curvature complexity as

\[
C(\epsilon) = \frac{\epsilon^2}{2\sigma^2},
\]

with \(\sigma\) a normalization scale. In numerical work we set \(\sigma = 1\) unless otherwise specified.

We then define **perception gain** as an exponential penalty in curvature complexity:

\[
P(\epsilon) = \exp\!\big(-\mu\,C(\epsilon)\big)
            = \exp\!\left(-\mu \frac{\epsilon^2}{2\sigma^2}\right),
\]

with \(\mu > 0\) a dimensionless parameter controlling how fragile recursive processing is to geometric anisotropy.

This structure is an **explicit ansatz**: a simple exponential in a quadratic complexity measure.

### 2.4 Total Recursive Processing \(T(N_e,H_I,\epsilon)\)

The **Total Recursive Processing (TRP)** capacity is defined as

\[
T(N_e,H_I,\epsilon)
= R(N_e,H_I) \, P(\epsilon)
= \log\!\left(\frac{A(N_e,H_I)}{4S_0}\right)
  \exp\!\left(-\mu \frac{\epsilon^2}{2\sigma^2}\right).
\]

We postulate a **TRP viability condition**:

\[
T(N_e,H_I,\epsilon) \;\ge\; T_{\min},
\]

where \(T_{\min} > 0\) is a universal threshold required for stable recursive observers.

---

## 3. TRP complexity bound \(C_{\max}\) and \(|\epsilon|_{\max}\)

The TRP viability inequality

\[
\log\!\left(\frac{A}{4S_0}\right)\exp(-\mu C) \ge T_{\min}
\]

can be solved explicitly for the maximal allowed curvature complexity at fixed \((N_e,H_I)\).

Let

\[
R(N_e,H_I) = \log\!\left(\frac{A(N_e,H_I)}{4S_0}\right).
\]

Then the inequality is

\[
R \exp(-\mu C) \ge T_{\min}
\quad\Rightarrow\quad
\exp(-\mu C) \ge \frac{T_{\min}}{R}.
\]

Assuming \(R > T_{\min} > 0\) so that the right-hand side is in \((0,1)\), we take logarithms:

\[
-\mu C \ge \log\!\left(\frac{T_{\min}}{R}\right)
\quad\Rightarrow\quad
C \le \frac{1}{\mu} \log\!\left(\frac{R}{T_{\min}}\right).
\]

So the **TRP complexity bound** is

\[
C_{\max}(N_e,H_I)
:= \frac{1}{\mu} \log\!\left(
\frac{R(N_e,H_I)}{T_{\min}}
\right)
= \frac{1}{\mu} \log\!\left(
\frac{\log\!\big(A(N_e,H_I)/(4S_0)\big)}{T_{\min}}
\right),
\]

defined whenever \(R(N_e,H_I) > T_{\min}\) and positive there.

Using \(C(\epsilon) = \epsilon^2/(2\sigma^2)\), the maximal allowed anisotropy at fixed \((N_e,H_I)\) is

\[
|\epsilon|_{\max}(N_e,H_I)
= \sigma \sqrt{2C_{\max}(N_e,H_I)}.
\]

This is exactly what the Python engine computes.

---

## 4. Numerical calibration to data (Planck/BK + isotropy)

### 4.1. Pivot point and TRP saturation

We pick a **pivot point** \((N_{e,\star}, H_{I,\star})\) consistent with CMB data:

- \(N_{e,\star} = 60\)
- \(H_{I,\star} \approx 1.8\times 10^{-5} M_{\rm Pl}\) for \(r_\star = 0.03\).

At this pivot, we identify a conservative anisotropy scale

\[
|\epsilon_\star| = 0.02
\]

compatible with Planck’s bounds on quadrupolar statistical anisotropy parameters \(|g_\*| \lesssim \mathcal{O}(10^{-2})\).

We **calibrate** \(\mu\) by demanding that TRP **saturates** the viability condition at the pivot anisotropy:

\[
T(N_{e,\star}, H_{I,\star}, \epsilon_\star) = T_{\min}.
\]

This is a model choice that ties \(\mu\) directly to a real constraint rather than choosing it arbitrarily.

### 4.2. Solving for \(\mu\)

Let

\[
R_\star
= R(N_{e,\star},H_{I,\star})
= \log\!\left(\frac{S_{\rm geom,\star}}{S_0}\right)
= \log\!\left(\frac{A_\star}{4S_0}\right),
\]

with

\[
A_\star = A(N_{e,\star},H_{I,\star})
        = k_A \frac{e^{2 N_{e,\star}}}{H_{I,\star}^2}.
\]

The curvature complexity at the pivot anisotropy is

\[
C_\star = C(\epsilon_\star)
= \frac{\epsilon_\star^2}{2\sigma^2}.
\]

Saturation \(T_\star = T_{\min}\) gives

\[
R_\star \exp(-\mu C_\star) = T_{\min}
\quad\Rightarrow\quad
\exp(-\mu C_\star) = \frac{T_{\min}}{R_\star}.
\]

Therefore

\[
\boxed{
\mu = \frac{1}{C_\star}
      \log\!\left(\frac{R_\star}{T_{\min}}\right)
}
\]

provided \(R_\star > T_{\min} > 0\).

For typical numerical choices:

- \(k_A = 10^3\),
- \(S_0 = 10^5\),
- \(T_{\min} = 3\),
- \(\sigma = 1\),
- \(\epsilon_\star = 0.02\),
- \(N_{e,\star} = 60\), \(H_{I,\star} \approx 1.8 \times 10^{-5}\),

one finds:

- \(A_\star \sim 4 \times 10^{64}\),
- \(S_{\rm geom,\star} \sim 10^{64}\),
- \(R_\star = \log(S_{\rm geom,\star}/S_0) \sim \mathcal{O}(10^2)\),
- \(C_\star = (0.02)^2 / 2 = 2\times10^{-4}\),

and thus

\[
\mu \sim 10^4{-}10^5,
\]

with the precise value set by the exact numerical constants.

Once \(\mu\) is fixed, all **TRP predictions** for \(|\epsilon|_{\max}(N_e,H_I)\) follow from the explicit formulas above.

---

## 5. Python implementation

The file [`trp_engine.py`](trp_engine.py) implements the TRP engine:

- `area(Ne, HI)` – exit surface area \(A(N_e,H_I)\).
- `S_geom(Ne, HI)` – geometric entropy \(S_{\text{geom}}\).
- `R(Ne, HI)` – reality bandwidth \(R\).
- `C_max(Ne, HI)` – TRP complexity bound \(C_{\max}\).
- `epsilon_max(Ne, HI)` – maximal allowed anisotropy \(|\epsilon|_{\max}\).
- `T(Ne,HI,eps)` – TRP capacity \(T\).
- `is_viable(Ne,HI,eps)` – boolean viability under TRP.

On initialization, `TRPEngine`:

1. Stores \((S_0, T_{\min}, k_A, \sigma)\),
2. Takes a calibration pivot \((N_{e,\star}, H_{I,\star}, \epsilon_\star)\),
3. Solves for \(\mu\) using the formula above so that
   \[
   T(N_{e,\star},H_{I,\star},\epsilon_\star) = T_{\min}.
   \]

This ensures that the only free modeling choices are **explicit and interpretable**.

---

## 6. Usage

### 6.1 Install dependencies

The only required runtime dependency is `numpy`. For plotting (if you extend the code), add `matplotlib`.

```bash
pip install numpy

### 6.2 Minimal example

```python
from trp_engine import TRPEngine

engine = TRPEngine(
    S0=1e5,
    T_min=3.0,
    k_A=1e3,
    sigma=1.0,
    Ne_star=60.0,
    HI_star=1.8e-5,  # from r ~ 0.03, A_s ~ 2.1e-9
    eps_star=0.02    # calibrated to Planck-like anisotropy bound
)

print("Calibrated mu:", engine.mu)

# Pivot check: should return ~0.02
eps_max_star = engine.epsilon_max(60.0, 1.8e-5)
print("|eps|_max at pivot:", eps_max_star)

# Example variation: lower H_I at same N_e
for HI in [1.8e-5, 1.0e-5, 5.0e-6]:
    print(HI, engine.epsilon_max(60.0, HI))

### 6.3 Planck pivot example

A more complete example using the `hhv_math` utilities and Planck-like
parameters is provided in:

```bash
examples/pivot_from_planck.py

### 6.4 Grid scan example

For a more systematic calculation, see:

```bash
examples/scan_grid.py

## 8. References and source-backed formulas

This section lists the main external physics results used in this repo,
so every numerical input or standard relation is traceable to the
literature.

### 8.1 Slow-roll relations for r, A_s, H_I

We use the standard single-field slow-roll relations in Planck units
(\(M_{\rm Pl}=1\)):

- Tensor power:
  \[
  P_t = \frac{2 H_I^2}{\pi^2}.
  \]
- Scalar power at CMB pivot:
  \[
  P_s = A_s \simeq 2.1\times 10^{-9}.
  \]
- Tensor-to-scalar ratio:
  \[
  r = \frac{P_t}{P_s} = \frac{2 H_I^2}{\pi^2 A_s}
  \quad\Rightarrow\quad
  H_I = \frac{\pi}{\sqrt{2}}\sqrt{r A_s}.
  \]

These formulas appear in standard cosmology texts and in CMB review
papers (e.g. the Planck 2018 cosmological parameters and inflation
papers).

### 8.2 Observational inputs (example values)

1. **Scalar amplitude \(A_s\)**  
   From Planck 2018 base-\(\Lambda\)CDM fit:
   \[
   A_s \simeq 2.1\times 10^{-9}
   \]
   at pivot \(k = 0.05\,\mathrm{Mpc}^{-1}\).

2. **Tensor-to-scalar ratio \(r\)**  
   Combined Planck + BICEP/Keck analyses give upper bounds of order
   \[
   r \lesssim 0.03\quad\text{(95\% C.L., pivot \(k\sim 0.005\,\mathrm{Mpc}^{-1}\))}.
   \]
   We use \(r_\star = 0.03\) as a conservative high-scale pivot.

3. **Statistical isotropy (quadrupolar anisotropy)**  
   Planck statistical isotropy analyses constrain quadrupolar
   power-anisotropy parameters (often denoted \(g_\*\)) to
   \(|g_\*| \lesssim \mathcal{O}(10^{-2})\) at CMB scales, depending on
   the exact model. We therefore choose
   \[
   |\epsilon_\star| = 0.02
   \]
   as a conservative upper bound for a dimensionless anisotropy
   parameter comparable in magnitude to \(g_\*\).

These numbers are used **only** as calibration inputs; all code is
explicit about where they enter.

### 8.3 TRP ansatz and derived formulas

The TRP structure
\[
T(N_e,H_I,\epsilon)
= \log\!\left(\frac{S_{\rm geom}(N_e,H_I)}{S_0}\right)
  \exp\!\left(-\mu\,\frac{\epsilon^2}{2\sigma^2}\right)
\]
and the curvature-complexity ansatz
\[
C(\epsilon) = \frac{\epsilon^2}{2\sigma^2}
\]
are **model assumptions**.

Given these definitions and the TRP viability condition \(T\ge T_{\min}\),
we derive:

- Reality bandwidth:
  \[
  R(N_e,H_I)=\log\!\left(\frac{S_{\rm geom}}{S_0}\right)
            =\log\!\left(\frac{A(N_e,H_I)}{4S_0}\right).
  \]
- Complexity bound:
  \[
  C_{\max}(N_e,H_I)
  = \frac{1}{\mu}\log\!\left(\frac{R(N_e,H_I)}{T_{\min}}\right),
  \quad R > T_{\min}.
  \]
- Maximal anisotropy:
  \[
  |\epsilon|_{\max}(N_e,H_I)
  = \sigma \sqrt{2 C_{\max}(N_e,H_I)}.
  \]

The Python implementation in `trp_engine.py` is a direct translation of
these formulas.

### 8.4 Calibration of μ to data

At a chosen pivot \((N_{e,\star},H_{I,\star},\epsilon_\star)\) consistent
with CMB constraints, we impose TRP saturation:
\[
T(N_{e,\star},H_{I,\star},\epsilon_\star) = T_{\min}.
\]

Writing
\[
R_\star = R(N_{e,\star},H_{I,\star}),\quad
C_\star = C(\epsilon_\star),
\]
this is
\[
R_\star \exp(-\mu C_\star) = T_{\min}
\quad\Rightarrow\quad
\mu = \frac{1}{C_\star}\log\!\left(\frac{R_\star}{T_{\min}}\right).
\]

This is exactly what `TRPEngine._calibrate_mu()` computes.

## 9. Tests

Basic tests live in `tests/test_trp_engine.py`. They verify that:

1. **Pivot calibration**  
   At the Planck-like pivot \((N_e^\*, H_I^\*, \epsilon^\*)\),
   `epsilon_max` reproduces `eps_star` (within numerical tolerance).

2. **Monotonicity in \(H_I\)**  
   For fixed \(N_e\), lower \(H_I\) (larger exit area/entropy) allows
   larger \(|\epsilon|_{\max}\).

3. **TRP suppression for large anisotropy**  
   Choosing \(\epsilon \gg \epsilon^\*\) at the pivot drives
   \(T(N_e^\*, H_I^\*, \epsilon)\) below \(T_{\min}\), as expected.

To run the tests locally:

```bash
pip install pytest
pytest

3. Commit message: `Document tests`.
4. Commit.

---

## Step 3 – How you (or anyone) runs the tests later

On a machine with Python:

```bash
git clone https://github.com/<your-username>/HHV-TRP-Engine.git
cd HHV-TRP-Engine
pip install -r requirements.txt pytest
pytest

### TRP decomposition and gauge freedom

The TRP law is multiplicative:
\[
T = R \times P.
\]
Only \(T\) is physically constrained by the viability bound
\(T \ge T_{\min}\). Any rescaling
\[
R \rightarrow a\,R,\qquad P \rightarrow \frac{1}{a}\,P
\]
leaves \(T\) invariant. This defines a one-parameter "gauge
freedom" in how we assign difficulty to the external world (R)
versus internal processing (P).

In this repo we fix that freedom by convention:
- \(R(N_e,H_I)\) depends only on exit-geometry entropy,
- \(P(\epsilon)\) depends only on curvature complexity.

Given any two of \((T,R,P)\), the third is uniquely determined:
\[
R = \frac{T}{P},\quad P = \frac{T}{R},\quad \log T = \log R + \log P.
\]
