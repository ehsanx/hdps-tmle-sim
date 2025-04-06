# How Effective Are Machine Learning and Doubly Robust Estimators in Incorporating High-Dimensional Proxies to Reduce Residual Confounding?

Code and materials for a simulation study evaluating high-dimensional proxy adjustment using hdPS, machine learning (LASSO, Super Learner), and TMLE. The simulations assess bias, coverage, and the impact of model complexity in high-dimensional settings.

---

## 📁 Folder Structure

```
simulation/
├── code/              # All code files
├── data/              # Simulated datasets
├── frequent/          # Frequent exposure & outcome scenario
├── rare_exposure/     # Rare exposure scenario
├── rare_outcome/      # Rare outcome scenario
```

### 📂 Mapping of Folder Names

Each folder contains subfolders for methods. The folder names correspond to the following methods listed in **Table 1** of the manuscript. Each folder contains code for a specific method, and some folders may contain multiple files for different library sizes.

| **Folder Name**        | **Table 1 Method**                      |
|-------------------------|------------------------------------------|
| `DC.TMLE`              | DC.TMLE                                  |
| `hdPS.TMLE`            | hdPS.TMLE                                |
| `LASSO.TMLE`           | LASSO.TMLE                               |
| `hdPS.LASSO.TMLE`      | hdPS.LASSO.TMLE                          |
| `TMLE`                 | TMLE.ks and TMLE.u                       |
| `hdPS.SL`              | hdPS.SL                                  |
| `LASSO.SL`             | LASSO.SL                                 |
| `hdPS.LASSO.SL`        | hdPS.LASSO.SL                            |
| `SL.ks`                | SL.ks |
| `PS`                   | PS.ks and PS.u                           |
| `hdPS`                 | hdPS                                     |
| `LASSO`                | LASSO                                    |
| `hdPS.LASSO`           | hdPS.LASSO                               |
| `SL.u`                 | SL.u                                     |


---

## ⚙️ How to Run Simulations

### 1. Setup

- **R required** (version ≥ 4.1 recommended)
- Install dependencies:
  ```r
  install.packages(c(
  "autoCovariateSelection", # hdPS analysis
  "SuperLearner",           # Ensemble learning
  "tmle",                   # TMLE
  "glmnet",                 # LASSO and elastic-net
  "earth",                  # MARS
  "nnet",                   # Neural networks
  "kernlab",                # Kernel-based methods
  "xgboost",                # Extreme Gradient Boosting
  "WeightIt",               # PS weighting
  "cobalt",                 # Balance tables and plots
  "survey",                 # complex survey samples
  "data.table",             # Fast data manipulation
  "sandwich",               # Robust variance
  "dplyr",                  # Data manipulation
  "purrr",                  # Functional programming
  "parallel",               # Parallel computing
  "plyr"                    # Splitting/applying/combining
  ))
  ```

---

### 2. Choose a Simulation Scenario

Go to the relevant folder under `simulation/`:

- `frequent`: frequent exposure and outcome
- `rare_exposure`: rare exposure
- `rare_outcome`: rare outcome

---

### 3. Run Simulations for Each Method

#### 3.1 Methods Other Than SL

**Step 1: Navigate to the Method Folder**

- Inside each scenario folder, go to the relevant method folder.
- Folders may contain one or more files. File names ending in `_1`, `_3`, or `_4` indicate different SL library sizes:
  - `_1`: 1-learner library
  - `_3`: 3-learner library
  - `_4`: 4-learner library

**Step 2: Run the Script**

- Open and run the script to generate results.
- To control the number of simulation iterations, set:
  - `iterStart`
  - `iterStop`
- Total iterations = `iterStop - iterStart + 1`

**Step 3: View the Results**

- Individual results are saved to `simXres.Rds` files in the method folder. `X` = the simulation iteration index.
- These are aggregated into the `resultDF` object in memory.
- To save the combined output, manually add:
  ```r
  saveRDS(resultDF, file = paste0(method.under.consider, "/resultsDF.Rds"))
  ```
---

#### 3.2 SL Methods

**Step 1: Run the Corresponding TMLE Method**

- SL methods depend on output from a TMLE method.
- For example, to run `hdPS.SL` with a 1-learner library, first run the TMLE script ending in `_1`.
- For 3- or 4-learner libraries, run the TMLE script ending in `_3_4`.

**Step 2: Copy the Weight Folder**

- After running the TMLE script, a `weight` folder will be generated (e.g., `weight_hdps_tmle_3`).
- Copy this entire folder into the corresponding SL folder.
- Do **not** rename the folder.

**Step 3: Run the SL Script**

- Run the SL script to generate simulation results.
- Modify `iterStart` and `iterStop` to adjust iteration count, as above.

**Step 4: View the Results**

- Results are combined into `resultDF`, as above. Save manually if needed.

---

## 📘 Further Details

For a detailed explanation of method-specific workflows, including variable selection, model fitting, and effect estimation:

➡️ See [MethodWorkflows](MethodWorkflows.md)


---

## 📝 Notes

- The output `resultDF` is created in memory and not saved unless modified.
- Keep the file structure unchanged unless necessary.
- Folder names grouped as:
  - Double Cross-Fit TMLE (`DC.TMLE`)
  - TMLE Methods (`*.TMLE`)
  - Super Learner Methods (`*.SL` and `SL.ks`)
  - Standard Methods with Proxies (`PS.ks`, `hdPS`, `LASSO`, `hdPS.LASSO`)
  - Standard Methods without Proxies (`*.u`)
- Method names that serve multiple roles, like `TMLE` and `PS`, which are mapped to both `.ks` and `.u` versions.

  
---

## 📄 License

This project is licensed under the GPL-3.0 license.
