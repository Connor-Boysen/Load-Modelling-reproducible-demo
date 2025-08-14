# reproducible_internship_demo.R
# Purpose: Reproducible demo of internship-style analytics using simulated data.
# Features demonstrated:
#   - Data simulation for energy load vs weather
#   - Linear model (LM) with robust SE
#   - Mixed-effects model (LME) with random effects by site
#   - Generalized additive model (GAM) for nonlinearity
#   - LASSO regression with glmnet
#   - Segmented (piecewise) trend
#   - Diagnostics and quick model comparison
# Output: prints model summaries and writes plots to ./plots/

# ---- 0) Setup ----
set.seed(42)

pkg_needed <- c(
  "dplyr","tidyr","tibble","ggplot2","lme4","broom.mixed","glmnet",
  "mgcv","segmented","sandwich","lmtest","MASS","readr"
)
# Uncomment to auto-install (optional)
# install.packages(setdiff(pkg_needed, rownames(installed.packages())))
invisible(lapply(pkg_needed, require, character.only = TRUE))

# Create output folders if not present
dir.create("plots", showWarnings = FALSE)

# ---- 1) Simulate daily data for multiple sites ----
n_sites <- 8
days <- seq.Date(from = as.Date("2023-01-01"), to = as.Date("2024-12-31"), by = "day")
N <- length(days) * n_sites

sim <- tidyr::expand_grid(
  date = days,
  site = factor(paste0("Site_", seq_len(n_sites)))
) |>
  dplyr::mutate(
    dow   = weekdays(date),
    month = as.integer(format(date, "%m")),
    season = dplyr::case_when(
      month %in% c(12,1,2) ~ "Winter",
      month %in% c(3,4,5)  ~ "Spring",
      month %in% c(6,7,8)  ~ "Summer",
      TRUE                 ~ "Fall"
    ),
    season = factor(season, levels = c("Spring","Summer","Fall","Winter")),
    trend  = as.numeric(date - min(date)) # days since start
  )

# Weather drivers (simulate realistic patterns)
# Daily temperature (F): seasonal sinusoid + noise + site bias
site_bias <- rnorm(n_sites, 0, 2)
names(site_bias) <- levels(sim$site)

sim <- sim |>
  dplyr::mutate(
    doy = as.numeric(strftime(date, "%j")),
    tmp_base = 55 + 20 * sin(2*pi*(doy-30)/365),
    temp_F = tmp_base + site_bias[site] + rnorm(dplyr::n(), 0, 4),
    precip_in = pmax(0, rnorm(dplyr::n(), 0.1, 0.2)),
    cloud_cvr = pmin(100, pmax(0, rnorm(dplyr::n(), 50, 20))),
    wind_mph  = pmax(0, rnorm(dplyr::n(), 8, 3)),
    holiday   = as.integer(format(date, "%m-%d") %in% c("01-01","07-04","11-28","12-25")),
    # rare storm flags (simulate with small probabilities, vary by season)
    Floods = rbinom(dplyr::n(), 1, ifelse(season %in% c("Spring","Summer"), 0.01, 0.002)),
    Blizzards = rbinom(dplyr::n(), 1, ifelse(season == "Winter", 0.01, 0.0005)),
    Derechoes = rbinom(dplyr::n(), 1, ifelse(season == "Summer", 0.005, 0.0002)),
    Tornadoes = rbinom(dplyr::n(), 1, ifelse(season %in% c("Spring","Summer"), 0.003, 0.0001)),
    Severe_Thunderstorms = rbinom(dplyr::n(), 1, ifelse(season %in% c("Spring","Summer"), 0.01, 0.001))
  )

# Degree days (match the style HDD69, CDD73 from your notes)
base_hdd <- 69; base_cdd <- 73
sim <- sim |>
  dplyr::mutate(
    HDD69 = pmax(0, base_hdd - temp_F),
    CDD73 = pmax(0, temp_F - base_cdd)
  )

# True site random effects and coefficients (to recover with models)
re_site_intercept <- rnorm(n_sites, mean = 0, sd = 20000) # base load by site
re_site_temp      <- rnorm(n_sites, mean = 120,  sd = 20) # temp slope by site
names(re_site_intercept) <- names(re_site_temp) <- levels(sim$site)

# Construct outcome: two flavors (GDTC vs GETC) to mirror your projects
# Baseline: load driven by temp (both HDD and CDD), weather, season, trend
make_mu <- function(df) {
  500000 +                                      # system base
  re_site_intercept[as.character(df$site)] +    # site intercept
  (120*df$CDD73 +  90*df$HDD69) +               # nonlinear temp load
  - 200*df$cloud_cvr + 50*df$wind_mph +         # weather effects
  - 800*df$precip_in +
  ifelse(df$holiday==1, -15000, 0) +            # holidays slightly lower
  ifelse(df$Floods==1,  8000, 0) +              # storms push load up/down
  ifelse(df$Blizzards==1, 12000, 0) +
  ifelse(df$Derechoes==1, 6000, 0) +
  ifelse(df$Tornadoes==1,  4000, 0) +
  ifelse(df$Severe_Thunderstorms==1, 5000, 0) +
  60*df$trend                                   # mild growth trend
}

mu <- make_mu(sim)

# GDTC: "gross demand" (add site-specific temp slope)
gdtc <- mu + re_site_temp[as.character(sim$site)] * (sim$CDD73 + 0.6*sim$HDD69)

# GETC: "temperature corrected" (subtract a modeled temp component)
getc <- mu - (100*sim$CDD73 + 70*sim$HDD69)

# Add noise
sim <- sim |>
  dplyr::mutate(
    GDTC_MWh = gdtc + rnorm(dplyr::n(), 0, 12000),
    GETC_MWh = getc + rnorm(dplyr::n(), 0, 12000)
  )

# ---- 2) Train / test split ----
sim <- sim |> dplyr::arrange(date, site)
n_all <- nrow(sim)
idx_train <- seq_len(round(n_all*0.8))
train <- sim[idx_train, ]
test  <- sim[-idx_train, ]

# Helper: RMSE
rmse <- function(actual, pred) sqrt(mean((actual - pred)^2))

# ---- 3) Linear model (LM) ----
f_lm <- as.formula(GDTC_MWh ~ HDD69 + CDD73 + precip_in + cloud_cvr + wind_mph +
                     dow + holiday + season + Floods + Blizzards + Derechoes +
                     Tornadoes + Severe_Thunderstorms + trend)

lm_fit <- lm(f_lm, data = train)
cat("\n=== LM summary (GDTC) ===\n")
print(summary(lm_fit))

# Robust SE (HC3) for LM
cat("\n=== LM with robust SE (HC3) ===\n")
lm_robust <- lmtest::coeftest(lm_fit, vcov = sandwich::vcovHC(lm_fit, type = "HC3"))
print(lm_robust)

# Predictions & RMSE
pred_lm <- predict(lm_fit, newdata = test)
rmse_lm <- rmse(test$GDTC_MWh, pred_lm)
cat(sprintf("\nLM RMSE (test) = %.1f\n", rmse_lm))

# Residual plot
library(ggplot2)
p_resid <- ggplot(data.frame(fitted=fitted(lm_fit), resid=residuals(lm_fit)), aes(fitted, resid)) +
  geom_point(alpha=0.3) + geom_hline(yintercept = 0) +
  labs(title="Residuals vs Fitted — LM (GDTC)", x="Fitted", y="Residuals")
ggsave(filename = "plots/lm_residuals.png", plot = p_resid, width = 7, height = 4, dpi = 120)

# ---- 4) Mixed-effects model (LME) ----
# Random intercepts and temp slopes by site
f_lme <- as.formula(GDTC_MWh ~ HDD69 + CDD73 + precip_in + cloud_cvr + wind_mph +
                      dow + holiday + season + Floods + Blizzards + Derechoes +
                      Tornadoes + Severe_Thunderstorms + trend +
                      (1 + CDD73 + HDD69 | site))

lme_fit <- lme4::lmer(f_lme, data = train, REML = TRUE)
cat("\n=== LME summary (GDTC) ===\n")
print(summary(lme_fit))

# Tidy random effects
cat("\n=== Random effects (by site) ===\n")
print(broom.mixed::tidy(lme_fit, effects = "ran_vals") |> dplyr::slice_head(n = 10))

# Predictions & RMSE
pred_lme <- predict(lme_fit, newdata = test, allow.new.levels = TRUE)
rmse_lme <- rmse(test$GDTC_MWh, pred_lme)
cat(sprintf("\nLME RMSE (test) = %.1f\n", rmse_lme))

# ---- 5) GAM for nonlinearity ----
gam_fit <- mgcv::gam(GDTC_MWh ~ s(temp_F, k=10) + precip_in + cloud_cvr + wind_mph +
                       dow + holiday + season + s(trend, k=10),
                     data = train, method = "REML")
cat("\n=== GAM summary (GDTC) ===\n")
print(summary(gam_fit))

pred_gam <- predict(gam_fit, newdata = test)
rmse_gam <- rmse(test$GDTC_MWh, pred_gam)
cat(sprintf("\nGAM RMSE (test) = %.1f\n", rmse_gam))

# ---- 6) LASSO (glmnet) ----
# Build model matrix (no intercept); keep numeric/dummies automatically
x <- model.matrix(~ HDD69 + CDD73 + precip_in + cloud_cvr + wind_mph +
                    dow + holiday + season + Floods + Blizzards + Derechoes +
                    Tornadoes + Severe_Thunderstorms + trend,
                  data = train)[, -1]
y <- train$GDTC_MWh

cv <- glmnet::cv.glmnet(x, y, alpha = 1, family = "gaussian", nfolds = 5)
lasso <- glmnet::glmnet(x, y, alpha = 1, family = "gaussian", lambda = cv$lambda.min)
cat("\n=== LASSO nonzero coefficients (GDTC) ===\n")
print(sort(as.vector(coef(lasso)))[1:10])
pred_lasso <- predict(lasso, newx = model.matrix(~ HDD69 + CDD73 + precip_in + cloud_cvr + wind_mph +
                                                   dow + holiday + season + Floods + Blizzards + Derechoes +
                                                   Tornadoes + Severe_Thunderstorms + trend,
                                                 data = test)[, -1])
rmse_lasso <- rmse(test$GDTC_MWh, as.numeric(pred_lasso))
cat(sprintf("\nLASSO RMSE (test) = %.1f\n", rmse_lasso))

# ---- 7) Segmented regression on trend (piecewise) ----
lm_trend <- lm(GDTC_MWh ~ trend, data = train)
seg_fit <- tryCatch(
  segmented::segmented(lm_trend, seg.Z = ~ trend, psi = list(trend = median(train$trend))),
  error = function(e) { message("Segmented failed: ", e$message); NULL }
)
cat("\n=== Segmented trend break ===\n")
if (!is.null(seg_fit)) print(summary(seg_fit)) else print("No breakpoint estimated (fallback).")

# Plot segmented fit on a sample site to visualize
samp <- dplyr::filter(train, site == levels(train$site)[1]) |>
  dplyr::slice_sample(n = 400)
p_seg <- ggplot(samp, aes(trend, GDTC_MWh)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Segmented Trend (Illustrative)", x = "Trend (days)", y = "GDTC_MWh")
ggsave("plots/segmented_trend_example.png", p_seg, width = 7, height = 4, dpi = 120)

# ---- 8) Compare approaches on GDTC and also run LM on GETC ----
# Simple LM on GETC:
f_lm_getc <- update(f_lm, GETC_MWh ~ .)
lm_getc <- lm(f_lm_getc, data = train)
pred_getc <- predict(lm_getc, newdata = test)
rmse_getc <- rmse(test$GETC_MWh, pred_getc)

# Scorecard
score <- tibble::tibble(
  model = c("LM (GDTC)", "LME (GDTC)", "GAM (GDTC)", "LASSO (GDTC)", "LM (GETC)"),
  rmse  = c(rmse_lm, rmse_lme, rmse_gam, rmse_lasso, rmse_getc)
) |> dplyr::arrange(rmse)

cat("\n=== Test RMSE scorecard ===\n")
print(score)

# Save scorecard to CSV
readr::write_csv(score, "plots/model_scorecard.csv")

cat("\nAll done. Plots saved in ./plots and scorecard saved to plots/model_scorecard.csv\n")
