# Hartle–Hawking–Vireon TRP Engine: Notes for Physicists

This document summarizes the assumptions, equations, and intended use of the
HHV–TRP engine in a way aimed at cosmologists and high-energy theorists.

The goal is to make the code in this repository transparently interpretable and
easy to plug into existing Hartle–Hawking / tunneling wavefunction frameworks.

---

## 1. Setup and assumptions

We work within the standard semiclassical, slow-roll inflation framework:

- Metric: FLRW with small perturbations, plus possible homogeneous anisotropy
  at early times (e.g. Bianchi-type seeds).
- **Inflationary scale:** characterized by the Hubble parameter during inflation,
  \(H_I\), in reduced Planck units (with \(M_{\rm Pl} = 1\) in the code).
- **Observable patch:** specified by the number of e-folds between horizon
  exit of the CMB pivot scale and the end of inflation, \(N_e\).

On the quantum-cosmology side we assume:

- A **no-boundary** (Hartle–Hawking) or similar path-integral amplitude
  \(\Psi[g,\phi]\) over 4-geometries and matter fields, which induces a measure
  on exit 3-geometries \(\gamma\) (e.g. on a constant-density slice).
- For nearly FLRW universes, exit geometries can be characterized by:
  \[
  (N_e, H_I, \epsilon)
  \]
  where \(\epsilon\) is a dimensionless parameter encoding the effective
  anisotropy/inhomogeneity amplitude on the exit surface (see below).

This repository does **not** attempt to derive a specific Hartle–Hawking
wavefunction. Instead, it adds a **selection factor** on top of whatever
underlying amplitude one chooses.

---

## 2. TRP ansatz

We introduce a dimensionless quantity called **Total Recursive Processing** (TRP),
denoted \(T\), which is intended to encode the capacity of a universe to host
long-lived, self-measuring observers.

The basic ansatz is multiplicative:
\[
T = R \times P,
\]
where:

- \(R\) (“reality bandwidth”) quantifies the structural capacity of the geometry,
  and is taken to depend only on \((N_e,H_I)\).
- \(P\) (“perception gain”) quantifies how much of that capacity is usable given
  geometric complexity / anisotropy, and is taken to depend only on \(\epsilon\).

In this repository we make **one specific, explicit choice** for \(R\) and \(P\).

### 2.1 Reality bandwidth \(R(N_e,H_I)\)

We use a simple entropy-based proxy:
\[
A(N_e,H_I) = k_A\,\frac{e^{2 N_e}}{H_I^2},
\qquad
S_{\rm geom}(N_e,H_I) = \frac{A(N_e,H_I)}{4},
\]
with \(A\) the effective exit area and \(S_{\rm geom}\) the corresponding
geometric entropy in Planck units (assuming \(4 G \hbar = 1\)). The constants:

- \(k_A\): geometric prefactor (absorbing details of the exit surface choice
  and reheating history),
- \(S_0\): reference entropy scale representing a minimal “recursion engine”.

We define:
\[
R(N_e,H_I) = \log\!\left( \frac{S_{\rm geom}(N_e,H_I)}{S_0} \right)
= \log\!\left(\frac{A(N_e,H_I)}{4 S_0}\right).
\]

In the **code**, these appear as:

- `area(Ne, HI)` → \(A(N_e,H_I)\),
- `S_geom(Ne, HI)` → \(S_{\rm geom}(N_e,H_I)\),
- `R(Ne, HI)` → \(R(N_e,H_I)\).

### 2.2 Curvature complexity and perception gain

We introduce a single anisotropy parameter \(\epsilon\), which should be thought
of as a dimensionless measure of deviation from a round \(S^3\) exit geometry.
We do **not** fix its microscopic definition here; it is intended as the
effective amplitude that would be constrained by CMB isotropy measurements
(dipole/quadruple, hemispherical asymmetry, etc.).

We adopt a quadratic complexity ansatz:
\[
C(\epsilon) = \frac{\epsilon^2}{2 \sigma^2},
\]
with \(\sigma\) a scale parameter, and define the perception gain:
\[
P(\epsilon) = \exp\big(-\mu\,C(\epsilon)\big)
= \exp\!\left(-\mu\,\frac{\epsilon^2}{2 \sigma^2}\right).
\]

In the **code**, these appear as:

- `C(eps)` → \(C(\epsilon)\),
- `P(eps)` → \(P(\epsilon)\).

---

## 3. TRP and the viability condition

The total TRP capacity is:
\[
T(N_e,H_I,\epsilon) = R(N_e,H_I)\,\exp\big(-\mu\,C(\epsilon)\big).
\]

We impose a **minimal TRP requirement**:
\[
T(N_e,H_I,\epsilon) \;\ge\; T_{\min},
\]
where \(T_{\min}\) is a constant representing the minimal processing capacity
needed for long-lived, recursive observers.

This inequality is the **only** selection condition; there is no additional
anthropic tuning beyond the choice of numerical values for
\((S_0, T_{\min}, k_A, \sigma)\) and the calibration of \(\mu\).

In the **code**, the viability test is:

- `T(Ne, HI, eps)` → \(T(N_e,H_I,\epsilon)\),
- `is_viable(Ne, HI, eps)` → boolean mask for \(T \ge T_{\min}\).

---

## 4. Calibration of \(\mu\) at a Planck-like pivot

There is a one-parameter degeneracy in the decomposition \(T = R P\): rescaling
\[
R \rightarrow a R, \qquad P \rightarrow \frac{1}{a} P
\]
leaves \(T\) invariant. To make predictions, we fix this freedom by choosing a
**pivot point** in \((N_e, H_I, \epsilon)\) space.

We define a pivot:

- \(N_e^\*\): an effective number of e-folds at which CMB pivot modes exit
  (e.g. \(N_e^\*\sim 50\)–\(60\)),
- \(H_I^\*\): the corresponding inflationary Hubble rate, related to the tensor-to-scalar ratio \(r\) and scalar amplitude \(A_s\) via standard slow-roll relations,
- \(\epsilon^\*\): an anisotropy amplitude representative of (or bounded by) current CMB constraints.

Given \((N_e^\*, H_I^\*, \epsilon^\*)\), we **define** \(\mu\) by imposing:
\[
T(N_e^\*, H_I^\*, \epsilon^\*) = T_{\min}.
\]

Explicitly:
\[
R_\* = R(N_e^\*, H_I^\*),
\qquad
C_\* = C(\epsilon^\*),
\]
\[
R_\* \exp(-\mu C_\*) = T_{\min}
\quad\Rightarrow\quad
\mu = \frac{1}{C_\*}\,\log\!\left(\frac{R_\*}{T_{\min}}\right).
\]

This calibration is implemented in `TRPEngine._calibrate_mu()`. Once \(\mu\)
is fixed, the dependence of \(T\) and the resulting anisotropy bound on
\((N_e,H_I)\) is completely determined by the TRP ansatz and the chosen pivot.

The slow-roll relations and Planck-like parameters live in `hhv_math/`:

- `inflation.H_from_r_As(r, As, Mpl=1.0)` for \(H_I(r,A_s)\),
- `inflation.r_from_H_As(H, As, Mpl=1.0)` for the inverse,
- `cmb_params.*` contains representative \(A_s\), \(r\), \(N_e^\*\) values.

---

## 5. Complexity and anisotropy bounds

From the viability condition
\[
R(N_e,H_I)\,\exp\big(-\mu\,C(\epsilon)\big) \ge T_{\min},
\]
and assuming \(R(N_e,H_I) > T_{\min}\), we obtain:
\[
C(\epsilon) \le \frac{1}{\mu}
\log\!\left(\frac{R(N_e,H_I)}{T_{\min}}\right)
\;:=\; C_{\max}(N_e,H_I).
\]

With the quadratic form \(C(\epsilon) = \epsilon^2/(2\sigma^2)\), this becomes:
\[
\frac{\epsilon^2}{2 \sigma^2} \le C_{\max}(N_e,H_I)
\quad\Rightarrow\quad
|\epsilon| \le \sigma\,\sqrt{2\,C_{\max}(N_e,H_I)}
\;:=\; \epsilon_{\max}(N_e,H_I).
\]

The **code** implements:

- `C_max(Ne, HI)` → \(C_{\max}(N_e,H_I)\),
- `epsilon_max(Ne, HI)` → \(\epsilon_{\max}(N_e,H_I)\).

Given any early-universe model that predicts an effective exit anisotropy
\(\epsilon_{\rm model}(N_e,H_I)\), HHV–TRP declares it **TRP-viable** iff:
\[
|\epsilon_{\rm model}(N_e,H_I)| \le \epsilon_{\max}(N_e,H_I).
\]

---

## 6. Impact on early-universe models

### 6.1 High-scale, short-duration inflation

From
\[
A(N_e,H_I) \propto \frac{e^{2 N_e}}{H_I^2}
\]
we see that high \(H_I\) (large inflationary scale) and modest \(N_e\) imply
smaller exit area, smaller geometric entropy, and thus smaller \(R(N_e,H_I)\).

For smaller \(R\) at fixed \(T_{\min}\), the allowed complexity
\[
C_{\max}(N_e,H_I) = \frac{1}{\mu}\,\log\!\left(\frac{R(N_e,H_I)}{T_{\min}}\right)
\]
decreases, tightening the bound \(\epsilon_{\max}(N_e,H_I)\).

**Interpretation:**  
High-scale, short-duration inflationary histories (large \(r\), minimal e-folds)
are only TRP-viable if they produce extremely smooth exit geometries. Models
with significant primordial anisotropy in that regime are disfavored.

### 6.2 Low-scale, long-duration inflation

Conversely, low \(H_I\) and/or large \(N_e\) increase the exit area and
\(S_{\rm geom}\), hence increase \(R(N_e,H_I)\). This yields larger
\(C_{\max}\) and therefore a larger allowed \(\epsilon_{\max}\).

**Interpretation:**  
Low-scale, long-duration inflation can tolerate more pre-inflation geometric
structure (e.g. relic anisotropy) while still yielding enough TRP to support
observers.

### 6.3 Anisotropic (Bianchi) models

If an anisotropic Hartle–Hawking or tunneling solution yields an effective
anisotropy parameter \(\epsilon_{\rm model}(N_e,H_I)\), one can test it against
HHV–TRP in the code:

```python
from trp_engine import TRPEngine
from hhv_math import inflation, cmb_params

# 1. Build an engine from a Planck-like pivot
As = cmb_params.A_s_planck_2018
r_star = cmb_params.r_upper_approx_0p005
Ne_star = cmb_params.N_e_pivot
HI_star = inflation.H_from_r_As(r_star, As, Mpl=1.0)

engine = TRPEngine(
    S0=1e5,
    T_min=3.0,
    k_A=1e3,
    sigma=1.0,
    Ne_star=Ne_star,
    HI_star=HI_star,
    eps_star=0.02,  # reference anisotropy
)

# 2. Model prediction at some (Ne_model, HI_model)
Ne_model = 60.0
HI_model = HI_star
eps_model = 0.01  # example model's anisotropy

# 3. TRP anisotropy bound
eps_max = engine.epsilon_max(Ne_model, HI_model)
viable = abs(eps_model) <= eps_max
