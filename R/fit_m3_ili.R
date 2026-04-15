# Model 3 (ILI): MVN AR(1) with shared rho
# Joint multivariate normal innovations, single shared autocorrelation.
# Uses shared rho to match the Kinsa M3 for fair comparison.

suppressPackageStartupMessages({
  library(rjags)
  library(tidyverse)
  library(coda)
  library(scoringRules)
})

.m3_ili_model_string <- "
model {
  for(i in 1:tmax) {
    flu[i]   ~ dnegbin(p1[i], r[1])
    rsv[i]   ~ dnegbin(p2[i], r[2])
    covid[i] ~ dnegbin(p3[i], r[3])
    ili[i]   ~ dnegbin(p4[i], r[4])

    p1[i] <- r[1] / (r[1] + lambda1[i])
    p2[i] <- r[2] / (r[2] + lambda2[i])
    p3[i] <- r[3] / (r[3] + lambda3[i])
    p4[i] <- r[4] / (r[4] + lambda4[i])

    lambda1[i] <- exp(phi[i,1])
    lambda2[i] <- exp(phi[i,2])
    lambda3[i] <- exp(phi[i,3])
    lambda4[i] <- exp(phi[i,4]) * num_patients[i]
  }

  # Seasonal means
  for(k in 1:4) {
    mu0[k] ~ dnorm(0, 1e-4)
    beta_sin[k] ~ dnorm(0, 1e-4)
    beta_cos[k] ~ dnorm(0, 1e-4)

    for(i in 1:tmax) {
      mu_t[i,k] <- mu0[k] + beta_sin[k] * sin52[i] + beta_cos[k] * cos52[i]
    }
  }

  # Single shared AR(1) autocorrelation
  rho_shared ~ dunif(0, 1)

  # First time point
  for(k in 1:4) {
    phi[1,k] ~ dnorm(mu_t[1,k], 0.01)
  }

  # Multivariate AR(1)
  for(i in 2:tmax) {
    for(k in 1:4) {
      mu_phi[i,k] <- mu_t[i,k] + rho_shared * (phi[i-1,k] - mu_t[i-1,k])
    }
    phi[i,1:4] ~ dmnorm(mu_phi[i,1:4], Omega[1:4,1:4])
  }

  Omega[1:4,1:4] ~ dwish(R[1:4,1:4], 5)
  Sigma[1:4,1:4] <- inverse(Omega[1:4,1:4])

  r[1] ~ dgamma(0.1, 0.1)
  r[2] ~ dgamma(0.1, 0.1)
  r[3] ~ dgamma(0.1, 0.1)
  r[4] ~ dgamma(0.1, 0.1)
}
"

.m3_ili_model_file <- file.path(tempdir(), "model3_mvn_ar1_ili.jags")
writeLines(.m3_ili_model_string, .m3_ili_model_file)


fit_m3_ili <- function(data_with_pop,
                       as_of_date,
                       n_holdout_specific = 3,
                       n_holdout_ili = 1,
                       n_adapt   = 2000,
                       n_burnin  = 20000,
                       n_iter    = 30000,
                       n_thin    = 15,
                       verbose   = FALSE) {

  as_of_date <- as.Date(as_of_date)

  df <- data_with_pop %>%
    filter(time <= as_of_date) %>%
    arrange(time)

  tmax <- nrow(df)
  if (tmax < 52) stop("Not enough data before as_of = ", as_of_date)

  t_vec <- 1:tmax
  sin52 <- sin(2 * pi * t_vec / 52)
  cos52 <- cos(2 * pi * t_vec / 52)

  actuals_holdout <- df %>%
    tail(n_holdout_specific) %>%
    select(time, flu_count, rsv_count, covid_count)

  flu_data     <- df$flu_count
  rsv_data     <- df$rsv_count
  covid_data   <- df$covid_count
  ili_data_vec <- df$ili_count

  flu_data[(tmax - n_holdout_specific + 1):tmax]   <- NA
  rsv_data[(tmax - n_holdout_specific + 1):tmax]   <- NA
  covid_data[(tmax - n_holdout_specific + 1):tmax] <- NA
  ili_data_vec[(tmax - n_holdout_ili + 1):tmax]    <- NA

  jags_data <- list(
    flu = flu_data, rsv = rsv_data, covid = covid_data,
    ili = ili_data_vec,
    num_patients = df$ili_number_patients_tested,
    tmax = tmax, sin52 = sin52, cos52 = cos52,
    R = diag(4)
  )

  inits <- list(
    list(".RNG.seed" = 123, ".RNG.name" = "base::Wichmann-Hill"),
    list(".RNG.seed" = 456, ".RNG.name" = "base::Wichmann-Hill"),
    list(".RNG.seed" = 789, ".RNG.name" = "base::Wichmann-Hill")
  )

  params <- c("flu", "rsv", "covid")

  run_jags <- function() {
    jm <- jags.model(
      file = .m3_ili_model_file, data = jags_data,
      n.chains = 3, n.adapt = n_adapt, inits = inits, quiet = !verbose
    )
    update(jm, n.iter = n_burnin, progress.bar = if (verbose) "text" else "none")
    coda.samples(jm, variable.names = params, n.iter = n_iter,
                 thin = n_thin, progress.bar = if (verbose) "text" else "none")
  }

  samples <- if (verbose) run_jags() else suppressMessages(suppressWarnings(run_jags()))
  samples_combined <- do.call(rbind, samples)

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
        model = "m3_ili", as_of = as_of_date, target_date = target_date,
        horizon = as.integer(round(as.numeric(difftime(target_date, as_of_date, units = "weeks")))),
        disease = dname, truth = truth, mean = mean(s), median = qs[["50%"]],
        q025 = qs[["2.5%"]], q25 = qs[["25%"]], q75 = qs[["75%"]], q975 = qs[["97.5%"]],
        crps = crps_sample(y = truth, dat = s)
      )
    })
  })

  out
}
