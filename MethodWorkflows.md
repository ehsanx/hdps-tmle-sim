
# Description of Code Implementation on Methods

This document outlines the implementation details of the simulation methods described in Table 1 of the manuscript.

---

## 🔧 General Setup

- Each method was evaluated over **1,000 simulated iterations**.
- In each iteration, a unique **plasmode dataset** of size 3,000 was created by randomly sampling from the full analytic dataset.
- The dataset included:
  - Simulated treatment and outcome
  - Investigator-specified baseline covariates (demographics and labs)
  - Proxy variables (except for methods designed to exclude proxies, e.g., `TMLE.u`)

---

## 🔁 Iteration Workflows by Method

### 🟣 DC.TMLE
- Uses all covariates and proxies in the analysis.
- Applied the `DC_tmle_g1_k()` function from the `CrossFit` package.
- Performed:
  1. Double cross-fitting with 3-fold cross-validation
  2. TMLE over 25 splits
  3. Median aggregation of risk difference (RD) and standard error (SE) estimates
- Each Super Learner library (1-, 3-, and 4-learner) was evaluated separately.

---

### 🔵 TMLE Methods (`hdPS.TMLE`, `LASSO.TMLE`, `hdPS.LASSO.TMLE`, `TMLE.ks`)
1. Perform variable selection on proxy variables (see **Proxy Selection Methods** below).
2. Create analysis dataset with baseline covariates and selected proxies.
3. Specify exposure model using covariates + proxies.
4. Estimate propensity scores using Super Learner (with specific library).
5. Apply TMLE to obtain RD and SE estimates.

---

### 🟢 SL Methods (`hdPS.SL`, `LASSO.SL`, `hdPS.LASSO.SL`, `SL.ks`)
1. Perform proxy variable selection (method-specific).
2. Construct analysis dataset with covariates and selected proxies.
3. Generate propensity score weights using `create.weights.sl()` (based on corresponding TMLE method).
4. Estimate RD and SE using `extract.res()`.

---

### 🟠 Standard Methods With Proxies (`hdPS`, `LASSO`, `hdPS.LASSO`, `PS.ks`)
1. Select proxy variables via hdPS or LASSO.
2. Construct dataset with covariates and selected proxies.
3. Estimate PS using `create.weights()`.
4. Estimate RD and SE using `extract.res()`.

---

### 🔴 Standard Methods Without Proxies (`TMLE.u`, `SL.u`, `PS.u`)
- Follow the corresponding workflow above (TMLE, SL, or standard) **excluding all proxy variables**.
- Use only baseline covariates throughout the pipeline.

---

## 🔍 Proxy Selection Methods

- **hdPS**: Uses the Bross formula to prioritize variables for bias reduction. Top 100 proxies selected.
- **LASSO**: Applies L1-regularized regression to select proxies with non-zero coefficients.
- **hdPS.LASSO**: Combines both:
  - First applies hdPS to identify top 100 proxies.
  - Then applies LASSO for further selection.
- **ks (kitchen-sink)**: Uses **all proxies**; no selection performed.

