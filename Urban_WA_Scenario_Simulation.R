# ============================================================
# Vegetation-induced changes in urban water availability (WA)
# using a conditional Vine Copula framework
#
# Description
# ------------------------------------------------------------
# This script quantifies vegetation-induced changes in
# precipitation (P), evapotranspiration (ET), and urban
# water availability (WA = P − ET) by fitting a conditional
# Vine Copula model.
#
# The model characterizes the dependence structure among:
#
#   1. Precipitation (P)
#   2. Evapotranspiration (ET)
#   3. Urban vegetation (LAIur)
#   4. Peri-urban vegetation (LAIperi)
#
# Four vegetation scenarios are evaluated:
#
#   • Both_Dynamic
#       Observed urban and peri-urban LAI trajectories.
#
#   • Urban_Change
#       Urban LAI fixed at +1 SD relative to the baseline,
#       peri-urban LAI fixed at baseline.
#
#   • Periurban_Change
#       Peri-urban LAI fixed at +1 SD relative to the
#       baseline, urban LAI fixed at baseline.
#
#   • Both_Fixed
#       Urban and peri-urban LAI fixed at baseline.
#
# Outputs
# ------------------------------------------------------------
# Daily simulations of:
#
#   • Precipitation (P)
#   • Evapotranspiration (ET)
#   • Water availability (WA = P − ET)
#
# for all cities and all scenarios.
#
# ============================================================
#Variables order:
#
#1. Precipitation (P)
#2. Evapotranspiration (ET)
#3. Urban LAI (LAIu)
#4. Peri-urban LAI (LAIr)
#
#Conditional variables:
#LAIu and LAIr
#
#Response variables:
#P and ET
rm(list = ls())
gc()

# ============================================================
# 1. Load required packages
# ============================================================

library(VineCopula)
library(CDVineCopulaConditional)
library(openxlsx)
library(readxl)

# ============================================================
# 2. Define study period
# ============================================================

START_DATE <- as.Date("1981-01-01")
END_DATE   <- as.Date("2023-12-31")

dates <- seq.Date(
  from = START_DATE,
  to   = END_DATE,
  by   = "day"
)

EXPECTED_N_DAYS <- 15705

if(length(dates) != EXPECTED_N_DAYS){
  stop(
    sprintf(
      "Unexpected number of days: %s (expected %s).",
      length(dates),
      EXPECTED_N_DAYS
    )
  )
}

# ============================================================
# 3. Harmonic reconstruction of daily LAI
# ============================================================
#
# Description:
# Monthly LAI observations are converted to continuous
# daily trajectories using a harmonic regression approach.
#
# A single-harmonic model is fitted when fewer than
# 12 monthly observations are available; otherwise
# a double-harmonic model is used.
#
# Inputs:
#   lai_daily_raw : LAI time series
#   dates         : Date vector
#
# Output:
#   Daily LAI series
#
# ============================================================

fit_harmonic_lai <- function(lai_daily_raw, dates) {
  
  if(length(lai_daily_raw) != length(dates)){
    stop("Length of lai_daily_raw does not match dates.")
  }
  
# ------------------------------------------------------------
# Convert leap-year day-of-year values to a 365-day calendar.
#
# Feb 29 is mapped to day 59 and all subsequent dates are
# shifted backward by one day.
# ------------------------------------------------------------
  doy365 <- function(d){
    doy <- as.integer(format(d, "%j"))
    is_leap <- as.integer(format(d, "%Y")) %% 4 == 0 &
      (as.integer(format(d, "%Y")) %% 100 != 0 | as.integer(format(d, "%Y")) %% 400 == 0)
    # after Feb 28 in leap year, shift doy back by 1 (so Feb 29 -> 59, Mar 1 -> 60, ..., Dec 31 -> 365)
    doy[is_leap & doy > 59] <- doy[is_leap & doy > 59] - 1L
    doy
  }
  
  month_index_check <- format(dates, "%Y-%m")
  unique_per_month <- tapply(lai_daily_raw, month_index_check, function(x) length(unique(x[!is.na(x)])))
  is_monthly_data <- mean(unique_per_month == 1, na.rm=TRUE) > 0.8
  
  month_index <- format(dates, "%Y-%m")
  lai_month <- tapply(lai_daily_raw, month_index, mean, na.rm=TRUE)
  
  month_mid_dates <- as.Date(paste0(names(lai_month), "-15"))
  doy_month <- doy365(month_mid_dates)
  
  lai_month <- as.numeric(lai_month)
  valid <- !is.na(lai_month)
  if(sum(valid) < 6){
    return(lai_daily_raw)
  }
  
  doy_month_valid <- doy_month[valid]
  use_double_harmonic <- sum(valid) >= 12
  
  if(use_double_harmonic){
    df_model <- data.frame(
      lai = lai_month[valid],
      cos1 = cos(2*pi*doy_month_valid/365),
      sin1 = sin(2*pi*doy_month_valid/365),
      cos2 = cos(4*pi*doy_month_valid/365),
      sin2 = sin(4*pi*doy_month_valid/365)
    )
    harmonic_model <- lm(lai ~ cos1 + sin1 + cos2 + sin2, data = df_model)
    
  } else {
    df_model <- data.frame(
      lai = lai_month[valid],
      cos1 = cos(2*pi*doy_month_valid/365),
      sin1 = sin(2*pi*doy_month_valid/365)
    )
    harmonic_model <- lm(lai ~ cos1 + sin1, data = df_model)
  }
  
  doy_daily <- doy365(dates)
  
  if(use_double_harmonic){
    df_predict <- data.frame(
      cos1 = cos(2*pi*doy_daily/365),
      sin1 = sin(2*pi*doy_daily/365),
      cos2 = cos(4*pi*doy_daily/365),
      sin2 = sin(4*pi*doy_daily/365)
    )
  } else {
    df_predict <- data.frame(
      cos1 = cos(2*pi*doy_daily/365),
      sin1 = sin(2*pi*doy_daily/365)
    )
  }
  
  lai_daily <- predict(harmonic_model, newdata = df_predict)
# Enforce non-negative LAI values
  lai_daily <- pmax(lai_daily, 0)

  
  attr(lai_daily, "is_monthly_data") <- is_monthly_data
  return(lai_daily)
}

# ============================================================
# 4. Empirical CDF transformation
# ============================================================
#
# The Vine Copula model requires all variables to be
# transformed into the uniform space [0,1].
#
# Precipitation is treated as a mixed distribution because
# of the occurrence of zero-precipitation days.
#
# Continuous variables (ET and LAI) are transformed using
# standard rank-based empirical probabilities.
#
# ============================================================
# ------------------------------------------------------------
# pobs_zero()
#
# Convert precipitation to pseudo-observations while
# preserving the probability mass at zero precipitation.
#
# Input:
#   x : precipitation time series
#
# Output:
#   Uniform pseudo-observations in [0,1]
#
# ------------------------------------------------------------
pobs_zero <- function(x){
  x_clean <- x[!is.na(x)]
  zero_prop <- mean(x_clean == 0)
  out <- rep(NA, length(x))

# Randomize the zero-precipitation component  
  if(zero_prop > 0){
    idx_zero <- which(x == 0 & !is.na(x))
    idx_pos  <- which(x > 0 & !is.na(x))
    if(length(idx_zero) > 0) out[idx_zero] <- runif(length(idx_zero), 0, zero_prop*0.99)
    if(length(idx_pos)  > 0){
      posvals <- x[idx_pos]
      r <- rank(posvals)/(length(posvals)+1)
      out[idx_pos] <- zero_prop + (1-zero_prop)*r
    }
  } else {
    r <- rank(x_clean)/(length(x_clean)+1)
    out[!is.na(x)] <- r
  }
  return(out)
}
# ------------------------------------------------------------
# pobs_std()
#
# Rank-based empirical CDF transformation for continuous
# variables (e.g., ET and LAI).
#
# ------------------------------------------------------------
pobs_std <- function(x){
  x_clean <- x[!is.na(x)]
  r_clean <- rank(x_clean)/(length(x_clean)+1)
  out <- rep(NA, length(x))
  out[!is.na(x)] <- r_clean
  return(out)
}

# ============================================================
# 5. Back-transformation to physical space
# ============================================================
# ------------------------------------------------------------
# inv_precip()
#
# Transform simulated copula probabilities back to
# precipitation values using the empirical distribution
# of historical observations.
#
# Zero precipitation probability is preserved.
#
# ------------------------------------------------------------
inv_precip <- function(u, p_hist){
# Replace missing probabilities
  u[is.na(u)] <- 0.5 
  u <- pmin(pmax(u, 1e-6), 0.999999)
  
  x_pos <- p_hist[p_hist>0 & !is.na(p_hist)]
  if(length(x_pos) == 0){
    return(rep(0, length(u)))
  }
  
  zero_prop <- mean(p_hist==0, na.rm=TRUE)
  out <- rep(0, length(u))
  idx <- which(u > zero_prop)
  if(length(idx)>0){
    uu <- (u[idx] - zero_prop) / (1 - zero_prop)
    uu <- pmin(pmax(uu, 1e-6), 0.999999)
    out[idx] <- quantile(x_pos, uu, type=8, na.rm=TRUE)
  }
  return(out)
}
# ------------------------------------------------------------
# inv_std()
#
# Transform simulated copula probabilities back to
# the original variable space using empirical quantiles.
#
# ------------------------------------------------------------
inv_std <- function(u, x_hist){
# Replace missing probabilities
  u[is.na(u)] <- 0.5  
  u <- pmin(pmax(u, 1e-6), 0.999999)
  
  x_hist_clean <- x_hist[!is.na(x_hist)]
  if(length(x_hist_clean) == 0){
    return(rep(NA, length(u)))
  }
  
  quantile(x_hist_clean, u, type=8, na.rm=TRUE, names=FALSE)
}

# ======================================================
# 6. Input data
# ======================================================

base_dir <- "G:/WA/"

p_file  <- file.path(base_dir, "P19812023.xlsx")
et_file <- file.path(base_dir, "ET19812023.xlsx")
lai_ur_file <- file.path(base_dir, "LAI19812023ur.xlsx")
lai_peri_file <- file.path(base_dir, "LAI19812023peri.xlsx")
# ------------------------------------------------------------
# read_all()
#
# Read and standardize input data files.
#
# Expected dimensions:
#   15705 rows (daily data from 1981–2023)
#
# ------------------------------------------------------------
read_all <- function(path){
  df_raw <- read_excel(path)
  df_raw <- df_raw[, -1]
  df_raw <- df_raw[-1, ]
  df <- as.data.frame(lapply(df_raw, as.numeric))
  df <- df[rowSums(!is.na(df))>0, ]
  if(nrow(df)==15704){ df[15705,] <- NA }
  if(nrow(df)>15705) df <- df[1:15705,]
  if(nrow(df)<15705) stop("Input file contains fewer than 15,705 rows.")
  return(df)
}

P_all  <- read_all(p_file)
ET_all <- read_all(et_file)
LAIur_all <- read_all(lai_ur_file)
LAIperi_all <- read_all(lai_peri_file)

city_ids <- colnames(P_all)
n_city <- length(city_ids)

outdir <- file.path(base_dir, "Vine_Final_SD1")
dir.create(outdir, showWarnings = FALSE)

scenarios <- c(
  "Both_Dynamic",
  "Urban_Change",
  "Periurban_Change",
  "Both_Fixed"
)
# Initialize output files for each simulation scenario
for(s in scenarios){
  write.xlsx(data.frame(Date=dates), file.path(outdir,paste0(s,"_P.xlsx")), overwrite=TRUE)
  write.xlsx(data.frame(Date=dates), file.path(outdir,paste0(s,"_ET.xlsx")), overwrite=TRUE)
  write.xlsx(data.frame(Date=dates), file.path(outdir,paste0(s,"_WA.xlsx")), overwrite=TRUE)
}

# ============================================================
# 7. City-level D-vine copula simulation
# ============================================================
# ------------------------------------------------------------
# simulate_city()
#
# Fit a conditional Vine Copula model for a single city
# and simulate precipitation (P), evapotranspiration (ET),
# and urban water availability (WA = P − ET) under four
# vegetation scenarios:
#
#   1. Both Dynamic
#   2. Urban Change
#   3. Peri-urban Change
#   4. Both Fixed
#
# ------------------------------------------------------------
simulate_city <- function(i){
  
  city <- city_ids[i]
  cat("\nProcessing city:", city, "\n")

  P  <- P_all[[i]]
  ET <- ET_all[[i]]
  LAIur_raw <- LAIur_all[[i]]
  LAIperi_raw <- LAIperi_all[[i]]
  
  LAIur <- fit_harmonic_lai(LAIur_raw, dates)
  LAIperi <- fit_harmonic_lai(LAIperi_raw, dates)

  # Standardize LAI relative to the 1981–1990 baseline period
  baseline <- format(dates,"%Y") %in% as.character(1981:1990)
  LAIur_mean <- mean(LAIur[baseline], na.rm=TRUE)
  LAIur_sd   <- sd(LAIur[baseline], na.rm=TRUE)
  LAIperi_mean <- mean(LAIperi[baseline], na.rm=TRUE)
  LAIperi_sd   <- sd(LAIperi[baseline], na.rm=TRUE)

  LAIur_std <- (LAIur - LAIur_mean) / LAIur_sd
  LAIperi_std <- (LAIperi - LAIperi_mean) / LAIperi_sd
 #Both_Dynamic
 #Observed urban and peri-urban LAI trajectories.
 #
 #Urban_Change
 #Urban LAI fixed at +1 SD relative to baseline;
 #peri-urban LAI fixed at baseline.
 #
 #Periurban_Change
 #Peri-urban LAI fixed at +1 SD relative to baseline;
 #urban LAI fixed at baseline.
 #
 #Both_Fixed
 #Urban and peri-urban LAI fixed at baseline.
  delta <- 1  # +1 SD
  
  # Define vegetation scenarios
  LAIur_dynamic <- LAIur_std
  LAIperi_dynamic <- LAIperi_std
  
  LAIur_fixed <- rep(0, length(dates))
  LAIperi_fixed <- rep(0, length(dates))
  
  LAIur_urban <- LAIur_fixed + delta
  LAIperi_urban <- LAIperi_fixed
  
  LAIur_peri <- LAIur_fixed
  LAIperi_peri <- LAIperi_fixed + delta

  # Transform variables to pseudo-observations
  u_LAIur_dynamic <- pobs_std(LAIur_dynamic)
  u_LAIperi_dynamic <- pobs_std(LAIperi_dynamic)
  
  u_LAIur_fixed <- pobs_std(LAIur_fixed)
  u_LAIperi_fixed <- pobs_std(LAIperi_fixed)
  
  u_LAIur_urban <- pobs_std(LAIur_urban)
  u_LAIperi_urban <- pobs_std(LAIperi_urban)
  
  u_LAIur_peri <- pobs_std(LAIur_peri)
  u_LAIperi_peri <- pobs_std(LAIperi_peri)
  
  # Fit conditional Vine Copula
  u_P  <- pobs_zero(P)
  u_ET <- pobs_std(ET)

  U  <- cbind(u_P, u_ET, u_LAIur_dynamic, u_LAIperi_dynamic)
  Uc <- U[complete.cases(U), ]
  if(nrow(Uc)<300){
    cat("Insufficient observations. Skipping city.\n")
    return(NULL)
  }
  
  vc <- CDVineCondFit(Uc, Nx=2, type="DVine", treecrit="BIC")
  
  cond <- list(
  Both_Dynamic = cbind(
    u_LAIur_dynamic,
    u_LAIperi_dynamic
  ),

  Urban_Change = cbind(
    u_LAIur_urban,
    u_LAIperi_fixed
  ),

  Periurban_Change = cbind(
    u_LAIur_fixed,
    u_LAIperi_peri
  ),

  Both_Fixed = cbind(
    u_LAIur_fixed,
    u_LAIperi_fixed
  )
)
  
  results <- list()
  # Conditional simulation
  for(s in scenarios){
    cat(" → Scenario:",s,"\n")
    U_sim <- CDVineCondSim(vc, cond[[s]])
    U_sim <- U_sim[,1:2]
    
    P_sim  <- inv_precip(U_sim[,1], P)
    ET_sim <- inv_std(U_sim[,2], ET)
    
    results[[s]] <- list(
      P  = P_sim,
      ET = ET_sim,
      WA = P_sim - ET_sim
    )
    
  }
  
 # Export simulation results
  for(s in scenarios){
    fP  <- file.path(outdir, paste0(s,"_P.xlsx"))
    fET <- file.path(outdir, paste0(s,"_ET.xlsx"))
    fWA <- file.path(outdir, paste0(s,"_WA.xlsx"))
    
    dP  <- read_excel(fP)
    dET <- read_excel(fET)
    dWA <- read_excel(fWA)
    
    dP[[city]]  <- results[[s]]$P
    dET[[city]] <- results[[s]]$ET
    dWA[[city]] <- results[[s]]$WA
    
    write.xlsx(dP,  fP, overwrite=TRUE)
    write.xlsx(dET, fET, overwrite=TRUE)
    write.xlsx(dWA, fWA, overwrite=TRUE)
  }
  
  cat("Completed:", city, "\n")
}

# ============================================================
# Run simulations for all cities
# ============================================================
for(i in 1:n_city){
  simulate_city(i)
  gc()
}
cat(
  "\n====================================================\n",
  "All simulations completed successfully.\n",
  "====================================================\n"
)
