# MainModel MATLAB Walkthrough (for colleague presentation)

This guide is designed as a **talk track** so you can walk colleagues from economic intuition to code execution without getting lost in MATLAB implementation details.

---

## 1) What this project is doing at a high level

The code solves a **stationary equilibrium OTC intermediation model** with:

- customer types (high- vs low-valuation)
- dealer search frictions
- dealer inventories / chain intermediation
- bargaining in customer–dealer and dealer–dealer trades
- implied prices, markups, chain-length distribution, and welfare

The core output set you can emphasize in discussion:

- chain length distribution and average chain length
- markup levels and markup decomposition by chain position
- expected bid/ask/interdealer prices
- implied yield spread
- welfare decomposition (customers vs dealers, and gains from trade)

---

## 2) Practical note on the PDF in this environment

I attempted to parse `rdz037 (1).pdf` directly in this environment, but text extraction tools/libraries were unavailable. I could recover metadata and internal labels only (e.g., sections, propositions, equations), but not full readable body text.

Recovered metadata indicates this is TeX output with title `OP-REST190039`, containing references to sections 2–5, theorem/proposition/lemma labels, and calibration table labels (e.g., `table.2`).

So the walkthrough below is grounded in the MATLAB implementation and in-paper label structure embedded in the PDF.

---

## 3) Execution order (the one-slide overview)

When you run `MainModel_main_mod(...)`, the model executes in this order:

1. `MainModel_functions` (define reusable anonymous functions)
2. `MainModel_datamoments` (empirical targets / moments)
3. `MainModel_paramset` (baseline parameters)
4. `MainModel_backsolvedemographics` (infer demographics from moments)
5. `MainModel_distributions` (equilibrium distributions / masses)
6. `MainModel_values` (customer/dealer reservation values)
7. `MainModel_markup` (markup moments + decomposition)
8. `MainModel_display` (diagnostic prints / figures)
9. `MainModel_price_and_welfare` (prices, spreads, welfare)

That sequence is exactly how you should narrate the code in your meeting.

---

## 4) File-by-file briefing notes

## `MainModel_main_mod.m` (driver)

Use this as your entry point in live demo.

- Accepts optional parameter overrides, e.g. changing bargaining powers or valuation spreads.
- Sets numerical controls (`TolX`, `gridsize`, chain cap `K=7`).
- Runs all submodules in sequence.
- Returns a `results` struct with key parameters and outputs.
- Wrapped in `try/catch`, so failed calibration points return a failure object instead of hard-crashing batch jobs.

**Presentation line:** “This file is the reproducible pipeline from calibration inputs to equilibrium and welfare outputs.”

---

## `MainModel_datamoments.m` (empirical targets)

This is where empirical discipline enters.

Key targets/inputs:

- Yield spread target (ABX source)
- Li–Schürhoff chain-length/markup moments up to chain length 7
- Mean chain length target and inventory duration target
- Aggregate market cap, trade size, turnover, retail base, and dealer counts

Derived moment objects include:

- empirical chain length frequencies
- average markup and markup-vs-chain covariance/slope
- `M0_over_M_balanced` baseline share used by parameter set

**Presentation line:** “All structural unknowns are pinned down to hit this block of observed market moments.”

---

## `MainModel_paramset.m` (baseline calibration)

Defines baseline bargaining and flow-utility parameters:

- `theta`, `theta1` (bargaining splits)
- `yL`, `yH`
- dealer utility-flow bounds (`deltaLow`, `deltaHigh`)
- optional match heterogeneity (`epsiHigh` etc.)
- `M0_over_M` linked to balanced condition from moments

**Presentation line:** “This is the baseline from table-style calibration; we only override these in robustness runs.”

---

## `MainModel_backsolvedemographics.m` (inversion step)

Back-solves demographic/search primitives from targets and accounting identities:

- masses of asset and asset-less dealers (`M1`, `M0`, normalized `m1`, `m0`, `m`)
- search/meeting intensity objects (`lambda`, `rho`)
- customer-type shares (`piH`, `piL`)
- transition/flow masses (`muh0`, `mul1`, etc.)

**Presentation line:** “Before solving dynamic values, we infer the equilibrium population geometry implied by data.”

---

## `MainModel_functions.m` (math primitives)

Contains reusable anonymous functions used everywhere else:

- `m1EqmFunc(...)`: cubic solution (Cardano-based) for key equilibrium mass
- `Phi0/Phi1` and derivatives: dealer inventory CDF-like objects
- hazard/intensity mappings (`lambda0/lambda1/Lambda`)
- conditional distributions of initial/chain positions
- penalty function for HLW-pattern violations

**Presentation line:** “Think of this as a local mathematical library feeding all remaining blocks.”

---

## `MainModel_distributions.m` (equilibrium distribution solve)

This module uses demographics + primitives to compute:

- equilibrium masses (`m1Eqm`, `m0Eqm` etc.)
- distribution objects over dealer states
- checks and numerical objects used by value and markup modules

**Presentation line:** “This is where static equilibrium composition is pinned down before solving values/prices.”

---

## `MainModel_values.m` (value-function / reservation-value block)

Builds beta-distributed heterogeneity transforms and solves a linear system for:

- `DWH` (high-valuation customer reservation value)
- `DWL` (low-valuation customer reservation value)
- baseline dealer value level

Then computes a grid approximation `DVmat` and interpolation `DVfunc(x)` used in pricing and welfare integrals.

**Presentation line:** “This is the dynamic-value core; all price formulas downstream load this interpolation.”

---

## `MainModel_markup.m` (markup mechanics)

Computes markup outcomes and chain decompositions such as:

- expected markup by chain length
- aggregate expected markup
- decomposition matrix (`DistributionOfMarkupMat`) across chain positions
- average chain length implied by model

This is the direct comparison block to LS empirical markup moments.

**Presentation line:** “If your colleagues ask ‘does the model fit markups?’, this is the file to open.”

---

## `MainModel_display.m` (diagnostics + plots)

Formats and prints key diagnostics and often draws figures (depending on settings).

Use this for live presentation output if you want a quick sanity-check dashboard.

---

## `MainModel_price_and_welfare.m` (prices and welfare)

Computes:

- frictionless benchmark (`p_walras`) and autarky values
- OTC bid/ask/interdealer expected prices
- implied yield spread from expected interdealer price
- welfare decomposition: autarky vs frictionless vs OTC; customer vs dealer shares

**Presentation line:** “This converts structural equilibrium objects into economically interpretable outcomes.”

---

## `MainModel_batch.m` (robustness sweeps)

Runs `MainModel_main_mod` on an `ndgrid` of parameter combinations and saves `MainModel_batch_results.mat`.

Useful for sensitivity analysis slides:

- bargaining power shifts
- heterogeneity toggles (`epsiHigh`)
- low-valuation utility changes (`yL`)

---

## `MainModel_yieldspread.m`

Minimal helper that computes expected interdealer price and implied yield spread.

---

## 5) Suggested live demo flow (20–30 minutes)

1. **Start with motivation/outcomes:** “We explain markup and chain length in OTC markets.”
2. **Open `MainModel_datamoments.m`:** show the empirical anchors.
3. **Open `MainModel_paramset.m`:** show baseline calibrated knobs.
4. **Open `MainModel_main_mod.m`:** walk through call sequence.
5. **Jump to `MainModel_values.m` + `MainModel_markup.m`:** show where markups actually come from.
6. **Finish with `MainModel_price_and_welfare.m`:** show policy-relevant outcomes.
7. **If asked robustness:** show `MainModel_batch.m` and how to alter vectors.

---

## 6) Talking-point cheat sheet (for likely colleague questions)

**Q: Which parameters are most economically interpretable?**
- Bargaining (`theta`, `theta1`), valuation gap (`yH - yL`), dealer-flow bounds (`deltaLow`, `deltaHigh`), search intensities (`rho`, `lambda`).

**Q: Where do chain-length moments enter?**
- In `MainModel_datamoments.m`, then propagated through demographics/distribution and compared in markup outputs.

**Q: What creates bid-ask-interdealer wedges?**
- Search frictions + bargaining + heterogeneous reservation values (`DVfunc`, `DWH`, `DWL`) + inventory-state distributions (`Phi0`, `Phi1`).

**Q: How to run one counterfactual quickly?**
- `res = MainModel_main_mod(struct('theta',0.95,'theta1',0.6));`

**Q: How to run many?**
- Edit vectors in `MainModel_batch.m`; run once; inspect `MainModel_batch_results.mat`.

---

## 7) Minimal run commands (for your prep)

```matlab
% baseline
res0 = MainModel_main_mod(struct());

% quick bargaining-power sensitivity
res1 = MainModel_main_mod(struct('theta',0.95,'theta1',0.6));

% optional grid run
MainModel_batch
load('MainModel_batch_results.mat')
```

---

If you want, I can next produce a **10-slide speaking outline** (slide titles + what to say + which file/equation to display on each slide).
