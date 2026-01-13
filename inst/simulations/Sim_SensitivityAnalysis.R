# ─── Load Required Packages ─────────────────────────────────────────────────────
# quantreg: quantile regression functions (rq)
# tidyr:    data frame utilities (expand_grid)
library(tidyr)

# ─── Data Generation ────────────────────────────────────────────────────────────
# Generate recurrent‐event data up to subject‐specific censoring times.
# Args:
#   n     – number of subjects
#   max.C – numeric vector of length n giving each subject’s maximum follow-up time
# Returns:
#   data.frame with one row per observed episode, including start/end times,
#   episode length, covariates z1, z2, time-dependent count z3, and derived fields.
data.genrt <- function(n, max.C) {
  z1   <- runif(n) + 1
  z2   <- rbinom(n, 1, 0.5)
  d.int <- rnorm(n, mean = 0, sd = sqrt(1/2))
  all  <- NULL
  
  for (k in 1:n) {
    start <- 0
    end   <- 0
    obs   <- NULL
    indx  <- 0
    
    while (start < max.C[k]) {
      if (indx == 0) ht <- 0
      if (indx > 0)  ht <- obs[indx, 4]
      
      dW  <- rexp(1, 0.2 * exp(ht * 0.2))
      end <- start + dW
      curr.z3 <- indx
      Xt <- 10 + (d.int - z1 + z2)[k] +
        0.4 * curr.z3 +
        rnorm(1, mean = 0, sd = sqrt(((z1 + z2)[k] + 0.2 * curr.z3 + 1)^2 / 4 - 1/2))
      
      if (end + Xt < max.C[k]) {
        indx <- indx + 1
        obs  <- rbind(obs,
                      c(k, start, end, Xt, z1[k], z2[k], curr.z3, indx))
      }
      
      start <- end + Xt
    }
    
    if (indx >= 1) {
      n_epi <- indx
      T1    <- obs[1, 3]
      T2    <- obs[1, 3] + obs[1, 4]
      first <- rep(0, indx); first[1] <- 1
      last  <- rep(0, indx); last[indx] <- 1
      d.Xt  <- obs[, 4] - obs[1, 4]
      d.z3  <- obs[, 7] - obs[1, 7]
      
      obs.1 <- cbind(obs, first, last)
      obs.2 <- cbind(obs.1, n_epi)
      obs.3 <- cbind(obs.2, d.Xt, d.z3, T1, T2)
      all   <- rbind(all, obs.3)
    }
  }
  
  da <- data.frame(all)
  names(da) <- c(
    "id", "start", "end", "epi_length", "z1", "z2", "z3", "indx",
    "first", "last", "n_episode", "chgn_epi_length", "chgn_z3", "T1", "T2"
  )
  return(da)
}

# ─── Empirical Survival Function ────────────────────────────────────────────────
# Compute S_C(x) = P(C ≥ x) from censoring vector fix.C.
# Args:
#   x     – time point
#   fix.C – numeric vector of censoring times
# Returns:
#   Proportion of fix.C values ≥ x
S.C <- function(x, fix.C) {
  return(length(which(fix.C >= x)) / length(fix.C))
}

# ─── Censoring Time Generator ───────────────────────────────────────────────────
# Generate n normal draws truncated below at min.
# Args:
#   n   – number of draws
#   mu  – mean of normal
#   sd  – standard deviation
#   min – minimum allowed value
# Returns:
#   Numeric vector of length n
get.C <- function(n, mu, sd, min) {
  res <- rnorm(n, mu, sd)
  res[which(res < min)] <- min
  return(res)
}

# ─── Bootstrap Quantile Regression ──────────────────────────────────────────────
# Perform B bootstrap replicates of two-stage QR estimator.
# Args:
#   da    – data.frame from data.genrt
#   B     – number of bootstrap samples
#   fix.C – censoring vector for weights
# Returns:
#   List of: mean coefs, sd coefs, lower 95% CI, upper 95% CI
boot.qr <- function(da, B, fix.C) {
  coef.w3 <- NULL
  taus    <- c(0.25, 0.5, 0.75)
  
  for (i_B in 1:B) {
    id_B <- sample(unique(da$id), length(unique(da$id)), replace = TRUE)
    da_B <- NULL
    for (i_id in id_B) {
      da_B <- rbind(da_B, da[da$id == i_id, ])
    }
    
    fit.w3.2 <- quantreg::rq(chgn_epi_length ~ chgn_z3,
                   data = da_B[da_B$first == 0, ], tau = taus)
    
    da_B$n_epi_length_1 <- da_B$epi_length - fit.w3.2$coef[2, 1] * da_B$z3
    da_B$n_epi_length_2 <- da_B$epi_length - fit.w3.2$coef[2, 2] * da_B$z3
    da_B$n_epi_length_3 <- da_B$epi_length - fit.w3.2$coef[2, 3] * da_B$z3
    
    wc1 <- 1 / sapply(da_B$T1 + da_B$n_epi_length_1, S.C, fix.C = fix.C)
    wc2 <- 1 / sapply(da_B$T1 + da_B$n_epi_length_2, S.C, fix.C = fix.C)
    wc3 <- 1 / sapply(da_B$T1 + da_B$n_epi_length_3, S.C, fix.C = fix.C)
    
    fit1 <- quantreg::rq(n_epi_length_1 ~ z1 + z2,
               data = da_B, tau = taus[1],
               weights = 1 / n_episode * wc1)
    fit2 <- quantreg::rq(n_epi_length_2 ~ z1 + z2,
               data = da_B, tau = taus[2],
               weights = 1 / n_episode * wc2)
    fit3 <- quantreg::rq(n_epi_length_3 ~ z1 + z2,
               data = da_B, tau = taus[3],
               weights = 1 / n_episode * wc3)
    
    coef.w3 <- rbind(
      coef.w3,
      c(
        rbind(
          cbind(fit1$coef, fit2$coef, fit3$coef),
          fit.w3.2$coef[2, ]
        )
      )
    )
  }
  
  res_mean <- apply(coef.w3, 2, mean)
  res_sd   <- apply(coef.w3, 2, sd)
  res_low  <- apply(coef.w3, 2, quantile, probs = 0.025)
  res_high <- apply(coef.w3, 2, quantile, probs = 0.975)
  return(list(res_mean, res_sd, res_low, res_high))
}

# ─── Coverage Rate Calculators ─────────────────────────────────────────────────
# Gaussian approx. coverage
coverage.rate1 <- function(mean, sd, true) {
  res <- as.numeric(true <= mean + 1.96 * sd &
                      true >= mean - 1.96 * sd)
  return(sum(res) / length(res))
}
# Empirical interval coverage
coverage.rate2 <- function(low, high, true) {
  res <- as.numeric(true <= high & true >= low)
  return(sum(res) / length(res))
}

# ─── Main Simulation Routine ───────────────────────────────────────────────────
# Run 1000 simulations under two sample sizes and censoring levels,
# compare various QR estimators and return bias/sd table plus event rates.
# Args:
#   i_sim – index (1–4) selecting N = {250,500} and muC = {23,50}
# Returns:
#   List(sim.res: 12×12 matrix, sim.freq, sim.trun)
sim.RED <- function(i_sim) {
  N          <- c(250, 500)
  muC        <- c(23, 50)
  sim.setups <- expand_grid(N, muC)
  
  taus      <- c(0.25, 0.5, 0.75)
  true.coef <- t(cbind(
    qnorm(taus) / 2 + 10,
    qnorm(taus) / 2 - 1,
    qnorm(taus) / 2 + 1,
    0.4 + 0.2 * qnorm(taus) / 2
  ))
  
  size  <- as.numeric(sim.setups[i_sim, "N"])
  fix.C <- get.C(size, as.numeric(sim.setups[i_sim, "muC"]), 10, 10)
  
  coef.uw <- coef.w4 <- coef.w5 <- coef.w6 <- NULL
  sd.w6   <- mean.w6 <- low.w6 <- high.w6 <- NULL
  freq    <- trun <- NULL
  
  ctime <- Sys.time()
  for (sim in 1:1000) {
    da   <- data.genrt(n = size, max.C = fix.C)
    freq <- c(freq, nrow(da) / size)
    trun <- c(trun, 1 - length(unique(da$id)) / size)
    
    fit.uw <- quantreg::rq(epi_length ~ z1 + z2 + z3, data = da, tau = taus)
    coef.uw <- rbind(coef.uw, c(fit.uw$coefficients))
    
    fit.w4 <- quantreg::rq(epi_length ~ z1 + z2 + z3,
                 data = da, tau = taus,
                 weights = 1 / da$n_episode)
    coef.w4 <- rbind(coef.w4, c(fit.w4$coefficients))
    
    fit.w5.2 <- quantreg::rq(chgn_epi_length ~ chgn_z3,
                   data = da[da$first == 0, ], tau = taus)
    da$n_epi_length_1 <- da$epi_length - fit.w5.2$coef[2, 1] * da$z3
    da$n_epi_length_2 <- da$epi_length - fit.w5.2$coef[2, 2] * da$z3
    da$n_epi_length_3 <- da$epi_length - fit.w5.2$coef[2, 3] * da$z3
    
    f1 <- quantreg::rq(n_epi_length_1 ~ z1 + z2,
             data = da, tau = taus[1],
             weights = 1 / n_episode)
    f2 <- quantreg::rq(n_epi_length_2 ~ z1 + z2,
             data = da, tau = taus[2],
             weights = 1 / n_episode)
    f3 <- quantreg::rq(n_epi_length_3 ~ z1 + z2,
             data = da, tau = taus[3],
             weights = 1 / n_episode)
    coef.w5 <- rbind(coef.w5, c(
      rbind(cbind(f1$coef, f2$coef, f3$coef), fit.w5.2$coef[2, ])
    ))
    
    fit.w6.2 <- quantreg::rq(chgn_epi_length ~ chgn_z3,
                   data = da[da$first == 0, ], tau = taus)
    da$n_epi_length_1 <- da$epi_length - fit.w6.2$coef[2, 1] * da$z3
    da$n_epi_length_2 <- da$epi_length - fit.w6.2$coef[2, 2] * da$z3
    da$n_epi_length_3 <- da$epi_length - fit.w6.2$coef[2, 3] * da$z3
    
    wc1 <- 1 / sapply(da$T1 + da$n_epi_length_1, S.C, fix.C = fix.C)
    wc2 <- 1 / sapply(da$T1 + da$n_epi_length_2, S.C, fix.C = fix.C)
    wc3 <- 1 / sapply(da$T1 + da$n_epi_length_3, S.C, fix.C = fix.C)
    
    w1 <- quantreg::rq(n_epi_length_1 ~ z1 + z2,
             data = da, tau = taus[1],
             weights = 1 / n_episode * wc1)
    w2 <- quantreg::rq(n_epi_length_2 ~ z1 + z2,
             data = da, tau = taus[2],
             weights = 1 / n_episode * wc2)
    w3 <- quantreg::rq(n_epi_length_3 ~ z1 + z2,
             data = da, tau = taus[3],
             weights = 1 / n_episode * wc3)
    coef.w6 <- rbind(coef.w6, c(
      rbind(cbind(w1$coef, w2$coef, w3$coef), fit.w6.2$coef[2, ])
    ))
    
    br <- boot.qr(da, 200, fix.C)
    sd.w6   <- rbind(sd.w6, br[[2]])
    mean.w6 <- rbind(mean.w6, br[[1]])
    low.w6  <- rbind(low.w6, br[[3]])
    high.w6 <- rbind(high.w6, br[[4]])
  }
  Sys.time() - ctime
  
  tab1 <- matrix(nrow = 12, ncol = 12)
  tab1[, 1]  <- as.vector(true.coef)
  tab1[, 2]  <- as.vector(array(apply(coef.uw, 2, mean), c(4, 3)) - true.coef)
  tab1[, 3]  <- as.vector(array(apply(coef.uw, 2, sd), c(4, 3)))
  tab1[, 4]  <- as.vector(array(apply(coef.w4, 2, mean), c(4, 3)) - true.coef)
  tab1[, 5]  <- as.vector(array(apply(coef.w4, 2, sd), c(4, 3)))
  tab1[, 6]  <- as.vector(array(apply(coef.w5, 2, mean), c(4, 3)) - true.coef)
  tab1[, 7]  <- as.vector(array(apply(coef.w5, 2, sd), c(4, 3)))
  tab1[, 8]  <- as.vector(array(apply(coef.w6, 2, mean), c(4, 3)) - true.coef)
  tab1[, 9]  <- as.vector(array(apply(coef.w6, 2, sd), c(4, 3)))
  tab1[, 10] <- as.vector(array(apply(mean.w6, 2, mean), c(4, 3)) - true.coef)
  tab1[, 11] <- as.vector(array(apply(sd.w6, 2, mean), c(4, 3)))
  cr <- NULL
  for (i_col in 1:12) {
    cr <- c(cr, coverage.rate2(
      low  = low.w6[, i_col],
      high = high.w6[, i_col],
      true = true.coef[i_col]
    ))
  }
  tab1[, 12] <- cr
  colnames(tab1) <- c(
    "True", rep(c("bias", "sd"), 4),
    "mean bias", "estimated var", "coverage rate"
  )
  
  sim.res  <- tab1
  sim.freq <- mean(freq)
  sim.trun <- mean(trun)
  return(list(sim.res, sim.freq, sim.trun))
}
