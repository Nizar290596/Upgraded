# Thermophysical Coupling in mmcFoam

## Particle‑In‑Cell and Kernel Estimation, conditioned on the mixture fraction and on `phiModified`

**Prepared for committee review**
**Code base:** mmcFoam (sparse‑Lagrangian MMC for OpenFOAM, ESI/OpenCFD branch)
**Scope:** `src/lagrangian/mmc/submodels/Thermodynamic/ThermoPhysicalCouplingModel/`

---

## 1. Executive summary

mmcFoam solves turbulent reacting flows with a **hybrid Eulerian–Lagrangian**
formulation. A fluctuating scalar/thermochemical state is carried by a swarm of
stochastic **Lagrangian particles** (the Monte‑Carlo representation of the joint
composition PDF), while the mean flow and the mean thermochemical fields live on
the **Eulerian finite‑volume (FV) mesh**. For the two descriptions to remain
consistent, the particle thermochemistry must be **projected back onto the FV
grid** every time step so the density, temperature and species fields used by the
momentum/pressure solver reflect what the particles are actually doing. That
projection is the job of the **thermophysical coupling model**.

This report documents the two projection methods that are the subject of the
committee review:

| Method | Class | Idea in one line |
|---|---|---|
| **Particle‑In‑Cell (PIC)** | `ParticleInCell` | Straight **volume (super‑cell) averaging** of the particles that physically sit in each cell. |
| **Kernel Estimation (KE)** | `KernelEstimation` | **k‑nearest‑neighbour, kernel‑weighted regression** in *physical + conditioning space*, producing smooth, gap‑free target fields. |

Both methods produce the same output — Eulerian **target fields** for temperature
`TEqvETarget` and the solved species `YEqvETarget[i]` — which are then relaxed into
the flow through a source term. Their essential difference is **how they average**
and **what space they average in**. PIC averages only in physical space (cheap,
robust, but noisy and leaves holes where particles are sparse). KE additionally
conditions on a **coupling variable** — either the **mixture fraction `z`** or the
**modified progress variable `phiModified` (φ°)** — so that only particles that are
thermochemically *similar* to the target cell contribute, and it fills gaps by
kernel interpolation and extrapolation.

The report explains (i) the coupling problem and how the target enters the
governing equations, (ii) each model in detail, (iii) the two conditioning
variables (`z` and `phiModified`) and the Eulerian field `phiModEul` that makes
`phiModified`‑conditioning possible, and (iv) how both models interact with the
**second‑conditioning subset** of particles.

---

## 2. Background: why coupling is needed

### 2.1 The sparse‑Lagrangian MMC picture

In the Multiple Mapping Conditioning (MMC) / sparse‑Lagrangian approach:

- Each **particle** carries a full thermochemical state — mass `m`, statistical
  weight `wt`, temperature `T`, the composition vector `Y[i]`, and a set of
  **conditioning / reference variables** (position, mixture fraction, etc.).
  Chemistry and micro‑mixing act on the particles.
- The **FV mesh** carries the *equivalent* mean fields the pressure–velocity
  system needs. In the solver these are the **equivalent enthalpy** `hEqvE` and
  the **equivalent species** `YEqvE[i]` (and the temperature derived from them).

Because the particles are **sparse** (there may be far fewer or far more particles
than cells, and they are unevenly distributed), you cannot simply read a cell's
thermochemistry off "the particle in the cell". You must **estimate** the local
mean from the particle cloud. That estimation is the coupling model.

### 2.2 How the target enters the equations

The projected fields are not written into the mesh directly; they are used as
**relaxation targets**. From `baseMMCFoamEquations/hYEqvE_Eqn.H`:

```cpp
// once per time step (first PIMPLE iteration)
pSets.first().coupling().EqvETargetValues(Xi, YEqvETarget, TEqvETarget);
hEqvETarget = thermoEqvETarget.he(p, TEqvETarget);   // T-target -> enthalpy target

// then, for every solved species and for enthalpy:
YiEqn =
    fvm::ddt(rho, Yi) + mvConvectionYi->fvmDiv(phi, Yi)
  - fvm::laplacian(rho*DEff, Yi)
 == (rho*(YTarget - Yi)/tauRel) * Indicator;      //  <-- coupling source
```

Three features matter for this report:

1. **Relaxation source** `rho*(target − field)/tauRel`. The FV field is driven
   toward the particle‑derived target over a relaxation time `tauRel`. `tauRel`
   is itself blended from a large start value to a small target value over the
   run (`XiEqn.H`), so the coupling is loosened at start‑up and tightened as the
   solution matures.
2. The **`Indicator` field** (0/1 per cell) multiplies the source. **Where a
   coupling model could not produce a trustworthy target, it sets
   `Indicator = 0` and the source switches off** — the cell then evolves by pure
   transport/diffusion instead of being pulled toward a garbage target. This is
   the mechanism both models use to handle "no data here".
3. `EqvETargetValues(...)` is the **single virtual entry point** every coupling
   model implements. PIC and KE are two run‑time‑selectable implementations of
   it (`thermoPhysicalCouplingModel ParticleInCell;` or `... KernelEstimation;`).

So the committee should read both models as *different estimators feeding the same
relaxation source*, and the quality question is: **how faithful, how smooth, and
how well‑covered is the target field each one produces?**

---

## 3. The conditioning variables: `z` and `phiModified`

Before the models themselves, it is worth fixing the two variables the KE model
can condition on, because they are central to the request.

### 3.1 Mixture fraction `z`

The **mixture fraction** `z ∈ [0, 1]` measures the local degree of mixing between
fuel and oxidiser streams. It is a **conserved scalar** (no chemical source) and
is the classical MMC conditioning variable for **non‑premixed** combustion: at a
given `z`, the thermochemistry is strongly constrained, so *conditioning on `z`*
(averaging only particles of similar `z`) removes most of the scatter. In
`mmcVariablesDefinitions`, `z` is declared as a **passive** `couplingVar` — it is
transported on the FV grid with proper boundary conditions and simply carried
along:

```
z { name z; couplingName z; passive true; mixStep 1; MMCType couplingVar; }
```

Conditioning on `z` is the **default** mode of the Kernel Estimation model.

### 3.2 Modified progress variable `phiModified` (φ°)

For **premixed / stratified** combustion, `z` no longer discriminates the state —
reaction progress does. mmcFoam therefore introduces a **progress variable `φ`**
and a **modified progress variable `φ°` = `phiModified`**:

- **`φ` (`phi`)** is initialised from temperature,
  `φ = (T − T_u)/(T_b − T_u)` clamped to `[0,1]` (`MixingPopeCloud.C`,
  `setParticleProperties`), i.e. 0 in unburnt, 1 in fully burnt gas. It is
  advanced by a **progress‑variable reaction source**
  `W(φ) = A·(1−φ)·exp[Z·(φ−1)]` each time step (`updatePhiReaction`).
- **`φ° = φ·exp(β·ω_OU)`** is the **modified** progress variable. `ω_OU` is an
  **Ornstein–Uhlenbeck (OU)** stochastic process advanced per particle
  (`updateOUProcess`), and `β` scales how strongly the OU fluctuation perturbs
  `φ`. `φ°` can therefore **exceed 1** (typically up to ≈1.6 for `β≈0.15`),
  unlike `z` and `φ` which are bounded in `[0,1]`.

`φ°` is the physically appropriate conditioning variable for the
**second‑conditioning** subset (Section 6): conditioning on `φ°` groups particles
by *stochastically‑perturbed reaction progress*, which is what controls the
thermochemistry in a premixed flame brush.

### 3.3 Why `phiModified` needs a companion Eulerian field: `phiModEul`

Here is the subtlety that the recent code changes address, and that the committee
should understand clearly.

The Kernel Estimation model conditions in a mixed space **(x, y, z_phys, C)**,
where `C` is the coupling variable. To do that it needs **`C` on both sides**:

- the **particle** value of `C` (trivial — read it off each particle), **and**
- an **Eulerian field of `C` on the FV mesh**, because the kernel is evaluated
  *at each cell centre*, so it needs the cell's coupling‑variable coordinate to
  form the query point and to gate the flammability window `[fLow, fHigh]`.

For `z` this is free: `z` is a transported field that already lives on the mesh.
For `φ°` there is **no such field** — `φ°` is a *particle* quantity with **no
physical boundary conditions**, so it cannot simply be transported like `z`.

The solution is a dedicated **Eulerian phi‑degree field, `phiModEul`**, declared
as a **non‑passive** `couplingVar`:

```
phiModEul { name phiModEul; couplingName phiModEul; passive false; mixStep 1; MMCType couplingVar; }
```

`phiModEul` is **reconstructed from the particles** and then **relaxed**, once per
time step, inside `XiEqn.H`:

```cpp
// XiEqn.H  (phiModEul branch, first PIMPLE iteration)
XiiTarget = Xii;                          // default: no relaxation
forAllIters(pSets.first(), iter) {        // build particle-projected target
    if (filterFlagged && iter().secondCondFlag() != 1) continue;
    sumWt[celli]    += iter().m();
    sumPhiWt[celli] += iter().m() * iter().phiModified();
}
tgt[celli] = sumPhiWt[celli]/sumWt[celli];   // mass-weighted mean phi-degree
// then transported/relaxed:
XiiEqn = ddt(rho,Xii) + div(phi,Xii) - laplacian(rho*DEff,Xii)
      == rho*(XiiTarget - Xii)/tauRel;
```

Two design points to note:

- **Cells with no (subset) particle support keep `XiiTarget = Xii`** (zero
  relaxation) and are filled instead by **transport + diffusion** — this smoothly
  in‑paints `phiModEul` where particles are momentarily absent, exactly the gap
  problem that plagues raw particle projection.
- The usual `couplingVar` **upper clamp to 1 is skipped for `phiModEul`**, because
  `φ°` legitimately exceeds 1 (`Xii.min(1.0)` guarded by `if (!isPhiModEul)`).
  Only the lower clamp `Xii.max(0.0)` is applied.

So the complete `phiModified` coupling path is:

```
particle φ°  ──(mass-weighted projection in XiEqn)──►  phiModEul (Eulerian, relaxed & diffused)
     │                                                        │
     │ (particle coordinate)                                  │ (cell coordinate)
     └──────────────►  KernelEstimation kd-tree query  ◄──────┘
                              conditions on φ°
```

This is enabled with a single switch (Section 5.4): `conditionOnPhiModified true;`.

---

## 4. Model 1 — Particle‑In‑Cell (`ParticleInCell`)

**File:** `ParticleInCell.C` · **Author:** J. W. Gärtner (2022)

### 4.1 Concept

PIC is the **direct, physical‑space estimator**. For every cell it forms the
**mass/weight‑weighted average** of the particles that physically occupy the
cell's **super‑cell** and assigns that average as the cell target. A *super‑cell*
is a coarser aggregation of neighbouring FV cells (managed by
`particleNumberController`); pooling particles over a super‑cell increases the
sample size per average and reduces statistical noise relative to per‑cell
averaging.

There is **no conditioning variable** in PIC. It always conditions purely on
physical location. In the constructor it simply sets `Indicator() = 1` everywhere
— by default every cell is coupled.

### 4.2 Algorithm

From `EqvETargetValues` (`ParticleInCell.C`):

1. **Bucket particles by super‑cell** in a single `O(N_p)` pass:
   `pManager.getParticlesInSuperCellList(owner)`. (This is an explicit
   optimisation over the naive `O(N_cells × N_p)` approach of scanning every cell
   for its particles.)
2. **Zero the non‑solved species targets once** for the whole field (hoisted out
   of the per‑cell loop).
3. **For each super‑cell**, accumulate in a single pass over its particles:
   `totWt += wt`, `totT += T·wt`, and `totY[i] += Y[i]·wt` for solved species.
4. **Assign the average** to every member cell of the super‑cell:
   `TEqvETarget[cell] = totT/totWt`, `YEqvETarget[i][cell] = totY[i]/totWt`.
5. **Empty super‑cell → keep the previous time‑step value** (the `continue`
   branch). This deliberately avoids writing `0` into a cell that momentarily has
   no particles.

### 4.3 Interaction with the second‑conditioning subset

PIC is subset‑aware. When `secondCondMixingEnabled()` is true
(`filterFlagged`):

- It **starts from `Indicator = 0` everywhere** and only raises `Indicator = 1`
  on cells whose super‑cell contains at least one **flagged** particle
  (`secondCondFlag() == 1`).
- Only flagged particles contribute to the average (`if (... secondCondFlag() != 1)
  continue;`), because non‑subset particles carry stale `Y`/`T`.
- A super‑cell with **no flagged particle** (`totWt ≤ SMALL`) is **left uncoupled**
  (`Indicator` stays 0).

This is the key weakness the committee should note: **with a sparse flagged subset,
PIC leaves coupling holes** — every super‑cell without a flagged particle produces
no target and switches its source off. That is precisely the motivation for the
kernel estimator.

### 4.4 Strengths and limitations

| Strengths | Limitations |
|---|---|
| Very cheap (`O(N_p)`), no tree, no sort. | **No conditioning** on `z` or `φ°` — cannot separate thermochemically distinct states that overlap in space. |
| Local and conservative — a straight mass‑weighted mean. | **Noisy** when few particles per super‑cell. |
| Robust; no bandwidth/limit parameters to tune. | **Leaves gaps** (Indicator = 0) wherever particles/flagged particles are absent — severe for sparse or subset clouds. |
| Good baseline / reference. | Piecewise‑constant per super‑cell → blocky target fields. |

---

## 5. Model 2 — Kernel Estimation (`KernelEstimation`)

**File:** `KernelEstimation.C`

### 5.1 Concept

Kernel Estimation is a **conditioned, kernel‑weighted regression**. Instead of
"which particles are in this cell", it asks **"which particles are nearest to this
cell in the combined space of physical position and the coupling variable, and
what is their kernel‑weighted mean state?"** This delivers three things PIC cannot:

1. **Conditioning** — the coupling variable `C` (either `z` or `φ°`) is one of the
   search dimensions, so only thermochemically‑similar particles contribute.
2. **Smoothness** — a cubic (SPH‑style) smoothing kernel weights neighbours by
   distance, giving continuous target fields.
3. **Gap‑filling** — kernel evaluation at *every* cell centre, plus a
   gradient‑controlled **linear extrapolation** to neighbouring cells, keeps
   coverage high even for a sparse (flagged) subset.

### 5.2 Data structures and the mixed metric

Each particle and each cell is packed into a flat record (`densParticle =
List<scalar>`) with slots for index, `(x,y,z)`, weight, `T`, distance, the solved
species, and the coupling variables. Two lists are built each step:

- **`particleList_`** — the source data (`buildParticleList`): particle position,
  weight `m`, `T`, species `Y`, and coupling variables. In `φ°` mode the
  conditioning slot is filled from **`iter().phiModified()`**, not from a particle
  `XiC` slot; otherwise it is the particle's coupling value `iter().XiC()[j]`.
- **`LESList_`** — one record per FV cell (`buildLESParticleList`): cell centre and
  the cell's coupling value read from the **Eulerian coupling field** (`z`, or
  **`phiModEul`** in φ° mode).

A **4‑D kd‑tree** is built over the particles in dimensions
`(x, y, z_phys, C)` with per‑dimension weights `wts = {1, 1, 1, fm_}`
(`computeTargets`). The three spatial weights are unity; the **coupling‑space
weight is `fm`** — the bandwidth that sets *how much physical distance one unit of
coupling‑variable mismatch is worth*. Small `fm` ⇒ tight conditioning (particles
must match `C` closely); large `fm` ⇒ looser conditioning, wider acceptance
window (±0.5·`fm`).

### 5.3 The kernel algorithm, step by step

For each cell (walking the LES list, sorted by radial distance):

1. **Flammability gate.** If the cell's coupling value `fLES` is outside
   `[fLow, fHigh]`, skip it (`Indicator` stays 0). This avoids biased kernels at
   the edges of coupling space and saves work — e.g. for `z`, `[0.03, 0.19]`
   brackets the flammable band; for `φ°`, a band like `[0.0, 2.0]`.
2. **Recompute‑or‑extrapolate decision.** A full kernel is recomputed only when
   the cell is "far" from the last computed cell in coupling space
   (`|df| ≥ DELTAf`) or physical space (`dd ≥ DELTAd`); otherwise the target is
   obtained by cheap **linear extrapolation** from the last computed cell using
   the stored gradients `dY/dC`, `dT/dC`. `DELTAf`/`DELTAd` are derived from the
   local gradient and `C2`, and capped by `dfMax`. Setting `dfMax = 0` forces a
   full recompute everywhere (safest, most accurate, most expensive).
3. **k‑NN query.** `nNearest_` (default 20) nearest particles to the query point
   `(x_cell, y_cell, z_cell, fLES)` are retrieved from the kd‑tree.
4. **Kernel weights.** Two **cubic B‑spline smoothing kernels** (Monaghan,
   *Rep. Prog. Phys.* 68 (2005) 1703) are applied:
   - a **1‑D kernel in coupling space** (bandwidth `h = 0.25·fm`), giving weight
     `IDWf` and its derivative `dIDWf` (used to build the conditional gradient),
   - a **3‑D kernel in physical space** (bandwidth `h2 = 0.5·rMax`), giving
     `IDWd`.
   The combined particle weight is `Wt = IDWd · m · IDWf` (physical kernel ×
   particle mass × conditioning kernel).
5. **Weighted means and conditional gradients.** Accumulate `sumWt`, `sumWt·T`,
   `sumWt·Y[i]` (and the `d…` sums with `dIDWf`), then
   `T_target = ΣWt·T / ΣWt`, `Y_target[i] = ΣWt·Y[i] / ΣWt`, and the
   coupling‑gradients `dT/dC`, `dY[i]/dC`.
6. **Adaptive step control.** The largest normalised gradient sets `DELTAf`
   (`= C2·h/maxGrad`, capped by `dfMax`), which controls how far the next
   extrapolation may reach — steeper local variation ⇒ shorter reuse distance.
7. **Mark covered.** `Indicator[cell] = 1`. Inert species is set last so all
   species sum to one (`YEqvETarget[N2] = 1 − Σ others`, floored at 0).

### 5.4 Conditioning mode selection (`z` vs `phiModified`)

The mode is chosen in the constructor, **by a dedicated switch, not by
`condVariable`** (deliberately — `condVariable` must stay a real coupling variable
because the chemistry consumes it too, e.g. `ReactingPopeParticle` passes
`XiC(cVarName())` to chemistry as the mixture fraction):

```cpp
phiModEnabled_ = coeffDict().lookupOrDefault("conditionOnPhiModified", false);

if (phiModEnabled_ && !XiC().cVarInXiC().found("phiModEul"))
    FatalError << "conditionOnPhiModified requires a 'phiModEul' couplingVar ...";

word condName = phiModEnabled_ ? "phiModEul" : cVarName();   // "z" by default
condSlotXiC_  = XiC().cVarInXiC()[condName];
condSlotXi_   = XiC().cVarInXi()[condName];
```

- **`z` mode (default):** conditions on the transported mixture‑fraction field;
  particle coordinate = `iter().XiC()[z-slot]`, cell coordinate = the `z` field.
- **`φ° mode` (`conditionOnPhiModified true`):** conditions on **`phiModEul`**;
  particle coordinate = **`iter().phiModified()`**, cell coordinate = the relaxed
  Eulerian `phiModEul` field built in `XiEqn.H`. Requires the `phiModEul`
  couplingVar to exist (hard `FatalError` otherwise) and requires second
  conditioning to be enabled.

### 5.5 Robustness features (recent hardening)

The committee should note several defensive measures added for the sparse/subset
case:

- **Keep every flagged particle.** Down‑sampling to ≈0.5·N_cells is applied *only*
  to the full cloud; with the second‑conditioning subset active, **no thinning**
  is done (`subSample = (!filterFlagged) && ...`) — the subset is already sparse
  and thinning would re‑introduce holes.
- **Empty‑list guard.** With a sparse subset a processor (plus its gathered
  neighbours) can hold **zero** flagged particles. The kd‑tree constructor
  dereferences `particles_[0]`, so building it on an empty list would segfault
  that rank and surface as an **MPI/InfiniBand collective timeout** on its peers.
  The code guards this (`haveParticles`), skips the kernel (cells stay uncoupled,
  `Indicator = 0`), and **still runs the collective coverage diagnostic on every
  rank** so parallel balance is preserved.
- **Short/empty k‑NN result guard.** `result.empty()` and `nFound = min(nn,
  result.size())` handle the case where fewer than `nNearest` particles exist.
- **Halo exchange.** Before computing, each rank appends the particle lists of its
  **neighbouring processors** (`processorPolyPatch`) so kernels near processor
  boundaries are not biased by the domain cut.
- **Coverage diagnostic.** Every step prints
  `KernelEstimation coupling: X/Y cells covered (Z%), N particles in kernel list`,
  which is the primary knob‑tuning signal (raise `fm` if coverage is low).

### 5.6 Strengths and limitations

| Strengths | Limitations |
|---|---|
| **Conditions** on `z` or `φ°` — separates thermochemically distinct states. | More expensive (kd‑tree build + k‑NN per cell), though extrapolation reuse mitigates it. |
| **Smooth** target fields (cubic kernels) and **conditional gradients** for extrapolation. | Several parameters to tune (`fm`, `nNearest`, `fLow/fHigh`, `dfMax`, `C2`, `rMax`). |
| **High coverage** even for a sparse flagged subset (kernel + extrapolation fill gaps). | `φ°` mode needs the extra `phiModEul` field + `XiEqn` relaxation machinery and second conditioning. |
| Parallel‑robust (halo exchange, empty‑rank guards, coverage diagnostic). | Bandwidth `fm` must match the coupling variable's range (much larger for `φ°` than `z`). |

---

## 6. The second‑conditioning subset (context for both models)

Both models change behaviour when **second conditioning** is active, so a brief
description closes the loop.

A configurable **fraction `R`** of particles is flagged into a subset
(`secondCondFlag = 1`, assigned in `setParticleProperties`). For those particles
mmcFoam evolves the extra stochastic machinery:

- **Progress‑variable reaction** `W(φ) = A·(1−φ)·exp[Z·(φ−1)]` advances `φ`
  (`updatePhiReaction`, applied to *all* particles for consistency).
- **OU process** advances `ω_OU` and recomputes `φ° = φ·exp(β·ω_OU)` for flagged
  particles only (`updateOUProcess`). Non‑subset particles keep `ω_OU = 0`, so for
  them `φ° = φ`.

When the subset is active, **both coupling models restrict their averaging to
flagged particles** (their `Y`/`T` are the physically meaningful ones), and the
Eulerian `phiModEul` projection in `XiEqn.H` likewise uses only flagged particles.
This is where the two models diverge sharply in quality:

- **PIC** simply drops the source (`Indicator = 0`) in every super‑cell without a
  flagged particle → **sparse, hole‑ridden coupling**.
- **KE** conditions on `φ°`, interpolates by kernel, extrapolates by gradient, and
  in‑paints `phiModEul` by diffusion → **dense, smooth coupling** on the same
  sparse subset.

---

## 7. Configuration reference

### 7.1 Selecting the model

```
// in constant/cloudProperties, cloud sub-dict
thermoPhysicalCouplingModel ParticleInCell;      // or:
thermoPhysicalCouplingModel KernelEstimation;
```

`ParticleInCell` takes **no coefficients**.

### 7.2 `KernelEstimationCoeffs` — conditioning on mixture fraction `z` (default)

```
KernelEstimationCoeffs
{
    nNearest   20;      // flagged particles per cell kernel
    rMax       1e9;     // physical kernel-radius cap [m] (1e9 = uncapped)
    dfMax      0.0;     // 0 -> always recompute (no extrapolation reuse; safest)
    C2         0.1;     // extrapolation resolution control
    debug      false;

    fLow       0.03;    // lower flammability limit (in z)
    fHigh      0.19;    // upper flammability limit (in z)
    fm         0.01;    // conditioning bandwidth (acceptance window = +/- 0.5*fm)
}
```

### 7.3 `KernelEstimationCoeffs` — conditioning on `phiModified` (φ°)

```
KernelEstimationCoeffs
{
    nNearest   20;
    rMax       1e9;
    dfMax      0.0;
    C2         0.1;

    conditionOnPhiModified true;    // <-- switch to phi-degree mode
    fLow   0.0;   fHigh  2.0;        // phi-degree band (phi° can exceed 1)
    fm     0.3;                      // wider bandwidth (phi° range >> z range)
}
```

Additional requirements for φ° mode:

1. Add the **non‑passive `phiModEul` couplingVar** to `mmcVariablesDefinitions`
   (Section 3.3) so `XiEqn.H` builds and relaxes the Eulerian phi‑degree field.
2. **Enable second conditioning** (the subset that carries `ω_OU`/`φ°`).
3. **Do not** set `condVariable phiModified` — leave `condVariable` a real coupling
   variable (`z`); the chemistry also consumes it. The switch, not `condVariable`,
   selects the mode.
4. Watch the `KernelEstimation coupling: X/Y cells covered` line and raise `fm` if
   coverage is low.

> The block above is documented as a **reference configuration** in
> `tests/Cases/Case-mixing/constant/cloudProperties`.

---

## 8. Side‑by‑side comparison

| Aspect | Particle‑In‑Cell | Kernel Estimation (`z`) | Kernel Estimation (`φ°`) |
|---|---|---|---|
| Averaging space | Physical (super‑cell) | Physical + `z` | Physical + `φ°` |
| Conditioning | none | mixture fraction | modified progress variable |
| Estimator | mass‑weighted mean | cubic‑kernel k‑NN regression + gradient extrapolation | same, on `φ°` |
| Cell coupling variable source | — | transported `z` field | relaxed Eulerian `phiModEul` field |
| Particle coupling coordinate | — | `XiC()[z]` | `phiModified()` |
| Gap behaviour | holes (Indicator = 0) where sparse | kernel + extrapolation fill | + `phiModEul` diffusion in‑paints |
| Cost | lowest (`O(N_p)`) | higher (kd‑tree + k‑NN) | higher |
| Target smoothness | blocky (per super‑cell) | smooth | smooth |
| Best for | dense clouds, baseline, robustness | non‑premixed (`z` discriminates state) | premixed/stratified + second conditioning |
| Tunables | none | `fm, nNearest, fLow/fHigh, dfMax, C2, rMax` | + `conditionOnPhiModified`, `phiModEul` field |

**Recommendation summary.** Use **PIC** as a cheap, robust baseline where the
cloud is dense and no conditioning is needed. Use **Kernel Estimation on `z`** for
non‑premixed combustion where the mixture fraction discriminates the
thermochemistry. Use **Kernel Estimation on `phiModified`** for premixed /
stratified combustion with the second‑conditioning subset, where reaction progress
(not `z`) governs the state and the sparse subset would otherwise leave PIC's
coupling full of holes.

---

## 9. Key source locations

| Item | File |
|---|---|
| Coupling model base / virtual API | `.../ThermoPhysicalCouplingModel/ThermoPhysicalCouplingModel/ThermoPhysicalCouplingModel.H` |
| Particle‑In‑Cell | `.../ThermoPhysicalCouplingModel/ParticleInCell/ParticleInCell.{H,C}` |
| Kernel Estimation | `.../ThermoPhysicalCouplingModel/KernelEstimation/KernelEstimation.{H,C}` |
| Target → source term (enthalpy/species) | `applications/solvers/mmc/baseMMCFoamEquations/hYEqvE_Eqn.H` |
| Eulerian `phiModEul` relaxation | `applications/solvers/mmc/SPFoam/XiEqn.H` |
| Coupling variable declarations | `applications/solvers/mmc/SPFoam/mmcVariablesDefinitions` |
| `φ`, `φ°`, OU process, subset flag | `.../clouds/Templates/MixingPopeCloud/MixingPopeCloud.C` |
| `phiModified` accessor / storage | `.../popeParticles/Templates/MixingPopeParticle/MixingPopeParticle{.H,I.H}` |
| Reference configuration | `tests/Cases/Case-mixing/constant/cloudProperties` |

---

## 10. Glossary

- **MMC** — Multiple Mapping Conditioning: the stochastic PDF combustion model.
- **Coupling model** — projects Lagrangian particle thermochemistry onto the FV
  mesh as relaxation targets (`TEqvETarget`, `YEqvETarget`).
- **`Indicator`** — per‑cell 0/1 field; gates the relaxation source on/off.
- **Super‑cell** — coarse aggregation of FV cells used by PIC to pool particles.
- **`z`** — mixture fraction, conserved scalar, `[0,1]`; default KE conditioning
  variable.
- **`φ` (`phi`)** — progress variable, `[0,1]`, from temperature; reaction‑driven.
- **`φ°` (`phiModified`)** — modified progress variable `φ·exp(β·ω_OU)`; can exceed
  1; the φ‑mode conditioning variable.
- **`phiModEul`** — non‑passive Eulerian field of `φ°`, reconstructed from
  particles and relaxed/diffused in `XiEqn.H`; the cell‑side conditioning field.
- **`ω_OU`** — Ornstein–Uhlenbeck stochastic process on subset particles.
- **Second conditioning** — a flagged particle subset (fraction `R`) carrying the
  `φ`/`ω_OU`/`φ°` machinery for premixed/stratified modelling.
- **`fm`** — kernel bandwidth in coupling space (acceptance window ±0.5·`fm`).
