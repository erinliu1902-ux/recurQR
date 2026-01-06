# recurQR

`recurQR` implements the two-step quantile regression method for **recurrent episode length** data
from:

> Yi Liu, Guillermo E. Umpierrez, and Limin Peng.  
> *Exploring the Heterogeneity in Recurrent Episode Lengths Based On Quantile Regression.*

The method fits the model

\[
Q_{X_{ij}}(\tau \mid Z_i, L_{ij}) = \beta_1(\tau)^\top Z_i + \beta_2(\tau)^\top L_{ij}
\]

and is designed for recurrent episode data with dependent truncation/censoring and informative
cluster size.

## Install (GitHub)

```r
# install.packages("remotes")
remotes::install_github("YOUR_GITHUB_USERNAME/recurQR")
```

## Quick start

```r
library(recurQR)

dat <- read_example_sim_data()

fit <- rqrecur(
  epi_length ~ z1 + z2 + z3,
  data = dat,
  id = "id",
  td_vars = "z3",
  tau = c(0.25, 0.5, 0.75),
  # In real use, provide cohort censoring/follow-up times here:
  censor_time = NULL
)

fit
summary(fit)
coef(fit)

newdat <- data.frame(z1 = 1.5, z2 = 0, z3 = 3)
predict(fit, newdata = newdat)

# Bootstrap (optional) for pointwise confidence intervals
fit <- rqrecur(
  epi_length ~ z1 + z2 + z3,
  data = dat,
  id = "id",
  td_vars = "z3",
  tau = seq(0.1, 0.9, by = 0.1),
  censor_time = NULL,
  boot = TRUE,
  B = 50,   # increase for real analyses
  seed = 1
)

# Plot coefficient functions with pointwise CIs (if bootstrapped)
plot(fit, parm = c("(Intercept)", "z1", "z2", "z3"))
```

## Expected input shape (based on the simulation data)

The fitting function expects a **long** data.frame with one row per episode and at minimum:

- an ID column (`id`)
- an episode length column (the response in `formula`)
- covariate columns used in `formula`
- an episode ordering column (recommended), e.g. `indx` in the simulations

Optional but recommended for truncation correction:

- `censor_time`: optional censoring/follow-up times. If `NULL`, `rqrecur()` will use a `censor_time` column in `data` (one unique value per subject) when available; otherwise it sets `G(t)=1`. You can also pass a numeric vector of cohort follow-up times, or a column name in `data`.
- `w1`: a column containing `W_{i1}` (time from entry to start of first episode); the simulations use `T1`

To handle datasets that include an incomplete (censored) final episode, use the `complete` flag:

- If `complete = TRUE` (default): all rows are assumed to be complete episodes.
- If `complete = FALSE`: provide an episode censoring indicator via `status` (a column name in `data`),
  and only complete episodes are kept (typically `status == 1`).

Example:

```r
fit <- rqrecur(
  epi_length ~ z1 + z2 + z3,
  data = dat,
  id = "id",
  td_vars = "z3",
  tau = c(0.25, 0.5, 0.75),
  complete = FALSE,
  status = "delta",         # 1=complete, 0=censored
  status_complete = 1
)
```

## Simulation code

The original simulation scripts included by the authors are shipped under:

```r
system.file("simulations", package = "recurQR")
```

They are **not** run automatically by the package (they are heavy by design), but are included for
reproducibility and as worked examples.

## Notes / limitations (current version)

- `rqrecur()` currently supports **additive** covariate terms only (no interactions or transformations),
  because the step-1 "delta" construction requires an unambiguous definition of the time-dependent covariate changes.
- Time-dependent covariates specified in `td_vars` must be numeric.

## License

MIT.