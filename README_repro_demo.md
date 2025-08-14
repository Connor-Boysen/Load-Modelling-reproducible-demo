# Reproducible Internship Modeling Demo

This repository contains a **fully reproducible** R script that simulates daily energy load data and compares several modeling approaches used during my internship. It mirrors a “**GETC vs GDTC**” analysis by generating two synthetic load targets and testing multiple models for accuracy and interpretability.

---

## Contents

```
reproducible_internship_demo.R   # Main script (simulates data, fits models, saves outputs)
plots/                           # Output folder created at runtime
  lm_residuals.png               # Residuals vs Fitted (LM)
  segmented_trend_example.png    # Example plot of segmented trend
  model_scorecard.csv            # Test RMSE by model
```

---

## What the demo does

- **Simulates realistic daily data** for multiple sites over ~2 years (seasonality, noise, and site effects).
- Computes **HDD69** and **CDD73**, along with **calendar dummies** (dow, holidays) and **storm/event flags** (e.g., Floods, Blizzards, Derechoes, Severe Thunderstorms).
- Generates two synthetic targets: **`GDTC_MWh`** and **`GETC_MWh`**, to emulate parallel load definitions.
- Trains and compares the following models:
  1. **Linear Model (LM)** with robust standard errors (HC3)
  2. **Linear Mixed-Effects (LME)** with random intercepts *and* temperature slopes by site (`lme4::lmer`)
  3. **Generalized Additive Model (GAM)** with smooths for temperature and trend (`mgcv::gam`)
  4. **LASSO** regularized regression with cross-validation (`glmnet::cv.glmnet`)
  5. **Segmented Trend** model to detect a breakpoint in long-term trend (`segmented::segmented`)

- Produces a **test RMSE scorecard** and several plots for diagnostics and communication.

---

## How to run

```r
# In a fresh R session, in a folder that contains the script:
install.packages(c(
  "dplyr","tidyr","tibble","ggplot2","lme4","broom.mixed",
  "glmnet","mgcv","segmented","sandwich","lmtest","MASS","readr"
))

source("reproducible_internship_demo.R")
```

- The script will create a `plots/` directory and write outputs there.
- To change simulation horizon, number of sites, or noise levels, edit the **config block** at the top of the script (e.g., `set.seed()`, date range, `n_sites`, event rates, degree-day bases).

> **Tip:** For reproducible environments, consider using `renv`:
> ```r
> install.packages("renv")
> renv::init()
> renv::install(c("dplyr","tidyr","tibble","ggplot2","lme4","broom.mixed",
>                  "glmnet","mgcv","segmented","sandwich","lmtest","MASS","readr"))
> ```

---

## Models compared

| Model | R Function | Key Features | Notes |
|---|---|---|---|
| Linear Model (LM) | `lm()` | Weather + calendar + events + trend | HC3 **robust SE** via `sandwich` + `lmtest` |
| Linear Mixed-Effects (LME) | `lme4::lmer()` | Random intercepts and temp slopes by `Site` | Captures between-site heterogeneity |
| GAM | `mgcv::gam()` | Smooths of `temp_F` and `trend` | Nonlinear response to temperature/seasonality |
| LASSO | `glmnet::cv.glmnet()` | Penalty-based variable selection | Cross-validated λ; uses model matrix of features |
| Segmented Trend | `segmented::segmented()` | Breakpoint in trend | Wrapped in `tryCatch()` for stability |

**Evaluation:** The script performs a chronological train/test split and reports **test RMSE** for each model in `plots/model_scorecard.csv`. Visual diagnostics (e.g., residuals vs fitted) are saved as PNG files.

---

## Interpreting results

- **LM (robust SE):** Coefficients show average effects; HC3 reduces false positives under heteroskedasticity. Check residual plots for curvature or variance patterns.
- **LME:** Random effects quantify site-specific baselines and temperature sensitivity. Look at variance components and BLUPs for operational insights.
- **GAM:** Smooth terms reveal nonlinear relationships (e.g., diminishing returns at high/low temps). Use `summary()`/`plot()` for partial effects.
- **LASSO:** Useful for feature selection and multicollinearity control. Inspect the selected coefficients at `lambda.min` or `lambda.1se`.
- **Segmented:** Detects potential structural breaks (e.g., policy or equipment change). Validate that inferred breakpoints align with domain context.

---

## Replacing simulated data with real data

1. Load your data (CSV or RData) early in the script in place of the simulation block.  
2. Ensure you have at least the following fields or can derive them:
   - `date` (Date), `Site` (factor), `temp_F` (numeric), `precip_in_per_day`, `avg_wind_spd`, `cloud_cvr`  
   - Optional: `Holiday` (0/1), event flags (`Floods`, `Blizzards`, `Derechoes`, `Severe_Thunderstorms`)  
3. Compute **HDD69** and **CDD73**:
   ```r
   HDD69 <- pmax(69 - temp_F, 0)
   CDD73 <- pmax(temp_F - 73, 0)
   ```
4. Create calendar dummies: day-of-week, `Summer/Fall/Winter` seasons, etc.  
5. Use **relative paths** (e.g., `data/input.csv`) rather than absolute `C:/...` paths.  
6. Re-run the models exactly as in the demo; the rest of the pipeline is unchanged.

---

## Why this approach

- Mirrors a real **operations analytics** workflow: weather-normalized load modeling, site heterogeneity, and structural change detection.
- Balances **interpretability and performance** across classic statistical models, mixed effects, smooths, and regularization.
- Fully reproducible, with explicit seeds, pinned package list, and deterministic outputs.

---

## License

You may adapt this demo for personal or academic use. For open-source distribution, consider adding an MIT license.

---

## Contact

Questions or ideas? Open an issue or reach out on LinkedIn/email.
