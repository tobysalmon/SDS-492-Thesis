# Refactored base Kinsa nowcast model as a callable function.
#
# Usage:
#   source("R/fit_base_kinsa.R")
#   preds <- fit_base_kinsa(data_with_pop, as_of_date = "2025-02-23")
#
# `data_with_pop` must already contain the columns:
#   time, flu_count, rsv_count, covid_count, percent_ill
# (i.e. the respnet + kinsa merge that Base_Model_Kinsa.Rmd produces).
#
# Returns a long-format tibble with one row per (disease x held-out week):
#   model, as_of, target_date, horizon, disease, truth,
#   mean, median, q025, q25, q75, q975, crps

suppressPackageStartupMessages({
  library(rjags)
  library(tidyverse)
  library(coda)
  library(scoringRules)
})

# --- JAGS model string (identical to Base_Model_Kinsa.Rmd) -------------------
.base_kinsa_model_string <- "
model {
  for(i in 1:tmax) {
    flu[i]   ~ dnegbin(p1[i], r[1])
    rsv[i]   ~ dnegbin(p2[i], r[2])
    covid[i] ~ dnegbin(p3[i], r[3])
    kinsa[i] ~ dlnorm(phi[i,4], tau_kinsa)

    p1[i] <- r[1] / (r[1] + lambda1[i])
    p2[i] <- r[2] / (r[2] + lambda2[i])
    p3[i] <- r[3] / (r[3] + lambda3[i])

    lambda1[i] <- exp(phi[i,1])
    lambda2[i] <- exp(phi[i,2])
    lambda3[i] <- exp(phi[i,3])
    lambda4[i] <- exp(phi[i,4])
  }

  for(k in 1:4) {
    mu0[k] ~ dnorm(0, 1e-4)
    beta_sin[k] ~ dnorm(0, 1e-4)
    beta_cos[k] ~ dnorm(0, 1e-4)
    prec.phi[k] ~ dgamma(0.01, 0.01)
    rho1[k] ~ dunif(0, 1)

    mu_t[1,k] <- mu0[k] + beta_sin[k] * sin52[1] + beta_cos[k] * cos52[1]
    phi[1,k] ~ dnorm(mu_t[1,k], prec.phi[k] * (1 - pow(rho1[k], 2)))

    for(i in 2:tmax) {
      mu_t[i,k] <- mu0[k] + beta_sin[k] * sin52[i] + beta_cos[k] * cos52[i]
      phi[i,k] ~ dnorm(mu_t[i,k] + rho1[k] * (phi[i-1,k] - mu_t[i-1,k]),
                        prec.phi[k])
    }
  }

  r[1] ~ dgamma(0.1, 0.1)
  r[2] ~ dgamma(0.1, 0.1)
  r[3] ~ dgamma(0.1, 0.1)
  tau_kinsa ~ dgamma(0.01, 0.01)
}
"

# Write model file once on source (in tempdir so it doesn't clutter the repo)
.base_kinsa_model_file <- file.path(tempdir(), "base_model_kinsa.jags")
writeLines(.base_kinsa_model_string, .base_kinsa_model_file)


fit_base_kinsa <- function(data_with_pop,
                           as_of_date,
                           n_holdout_specific = 3,
                           n_adapt   = 1000,
                           n_burnin  = 10000,
                           n_iter    = 20000,
                           n_thin    = 10,
                           verbose   = FALSE) {

  as_of_date <- as.Date(as_of_date)

  # Truncate data to simulate "today = as_of_date"
  df <- data_with_pop %>%
    filter(time <= as_of_date) %>%
    arrange(time)

  tmax <- nrow(df)
  if (tmax < 52) {
    stop("Not enough data before as_of = ", as_of_date,
         " (need at least 52 weeks, have ", tmax, ")")
  }

  # Harmonics for seasonality
  t_vec <- 1:tmax
  sin52 <- sin(2 * pi * t_vec / 52)
  cos52 <- cos(2 * pi * t_vec / 52)

  # Hold out the last n_holdout_specific weeks of respnet
  actuals_holdout <- df %>%
    tail(n_holdout_specific) %>%
    select(time, flu_count, rsv_count, covid_count)

  flu_data   <- df$flu_count
  rsv_data   <- df$rsv_count
  covid_data <- df$covid_count
  kinsa_vec  <- df$percent_ill

  flu_data[(tmax - n_holdout_specific + 1):tmax]   <- NA
  rsv_data[(tmax - n_holdout_specific + 1):tmax]   <- NA
  covid_data[(tmax - n_holdout_specific + 1):tmax] <- NA

  jags_data <- list(
    flu   = flu_data,
    rsv   = rsv_data,
    covid = covid_data,
    kinsa = kinsa_vec,
    tmax  = tmax,
    sin52 = sin52,
    cos52 = cos52
  )

  inits <- list(
    list(".RNG.seed" = 123, ".RNG.name" = "base::Wichmann-Hill"),
    list(".RNG.seed" = 456, ".RNG.name" = "base::Wichmann-Hill"),
    list(".RNG.seed" = 789, ".RNG.name" = "base::Wichmann-Hill")
  )

  params <- c("flu", "rsv", "covid")

  run_jags <- function() {
    jm <- jags.model(
      file     = .base_kinsa_model_file,
      data     = jags_data,
      n.chains = 3,
      n.adapt  = n_adapt,
      inits    = inits,
      quiet    = !verbose
    )
    update(jm, n.iter = n_burnin, progress.bar = if (verbose) "text" else "none")
    coda.samples(
      model          = jm,
      variable.names = params,
      n.iter         = n_iter,
      thin           = n_thin,
      progress.bar   = if (verbose) "text" else "none"
    )
  }

  samples <- if (verbose) run_jags() else suppressMessages(suppressWarnings(run_jags()))
  samples_combined <- do.call(rbind, samples)

  # Build long-format prediction tibble
  disease_info <- tibble::tribble(
    ~disease, ~truth_col,
    "flu",    "flu_count",
    "rsv",    "rsv_count",
    "covid",  "covid_count"
  )

  out <- map_dfr(seq_len(n_holdout_specific), function(h) {
    idx <- (tmax - n_holdout_specific) + h
    target_date <- actuals_holdout$time[h]

    map_dfr(seq_len(nrow(disease_info)), function(r) {
      dname <- disease_info$disease[r]
      truth <- actuals_holdout[[disease_info$truth_col[r]]][h]
      s <- samples_combined[, paste0(dname, "[", idx, "]")]
      qs <- quantile(s, c(0.025, 0.25, 0.5, 0.75, 0.975))
      tibble(
        model       = "base_kinsa",
        as_of       = as_of_date,
        target_date = target_date,
        horizon     = as.integer(round(as.numeric(difftime(target_date, as_of_date, units = "weeks")))),
        disease     = dname,
        truth       = truth,
        mean        = mean(s),
        median      = qs[["50%"]],
        q025        = qs[["2.5%"]],
        q25         = qs[["25%"]],
        q75         = qs[["75%"]],
        q975        = qs[["97.5%"]],
        crps        = crps_sample(y = truth, dat = s)
      )
    })
  })

  out
}
