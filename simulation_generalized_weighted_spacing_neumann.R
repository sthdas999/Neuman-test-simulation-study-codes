# ============================================================
# SIMULATION RESULTS
# Generalized Weighted Spacing Neumann Test
# Comparison with Neumann, KS and CvM
# ============================================================

rm(list = ls())

set.seed(20260905)

# ------------------------------------------------------------
# Basic settings
# ------------------------------------------------------------

alpha <- 0.05
Nrep  <- 5000

sample_sizes <- c(50, 100, 200, 500, 1000)

# ------------------------------------------------------------
# 1. Classical von Neumann statistic
# ------------------------------------------------------------

neumann_stat <- function(x) {
  
  n <- length(x)
  
  num <- sum(diff(x)^2)
  
  den <- sum((x - mean(x))^2)
  
  if (den <= 0) return(NA_real_)
  
  (n - 1) * num / (2 * den)
}


# ------------------------------------------------------------
# 2. Proposed weighted spacing Neumann statistic
# ------------------------------------------------------------

spacing_neumann <- function(x, r = 1, weight = "Uniform") {
  
  n <- length(x)
  
  if (r >= n - 1) return(NA_real_)
  
  # Order statistics
  xs <- sort(x)
  
  # First-order spacings
  D <- diff(xs)
  
  # r-th order spacings
  Dr <- xs[(r + 1):n] - xs[1:(n - r)]
  
  # Successive differences of r-th order spacings
  delta <- diff(Dr)
  
  m <- length(delta)
  
  # Locations at which weights are evaluated
  u <- seq(0, 1, length.out = m)
  
  # ----------------------------------------------------------
  # Weight functions
  # ----------------------------------------------------------
  
  if (weight == "Uniform") {
    
    w <- rep(1, m)
    
  } else if (weight == "Central") {
    
    w <- 1 - (2 * u - 1)^2
    
  } else if (weight == "Tail") {
    
    w <- (2 * u - 1)^2
    
  } else {
    
    stop("Unknown weighting scheme.")
  }
  
  # Weighted quadratic numerator
  numerator <- sum(w * delta^2)
  
  # Centering denominator
  denominator <- sum((delta - mean(delta))^2)
  
  if (denominator <= 0) return(NA_real_)
  
  numerator / denominator
}


# ------------------------------------------------------------
# 3. KS statistic
# ------------------------------------------------------------

ks_stat <- function(x) {
  
  ks.test(x, "punif")$statistic
}


# ------------------------------------------------------------
# 4. Cramer-von Mises statistic
# ------------------------------------------------------------

cvm_stat <- function(x) {
  
  # Cramer-von Mises statistic for U(0,1)
  
  n <- length(x)
  
  xs <- sort(x)
  
  i <- seq_len(n)
  
  1 / (12 * n) +
    sum((xs - (2 * i - 1) / (2 * n))^2)
}


# ============================================================
# NULL CALIBRATION
# ============================================================

# Critical values are calibrated separately for each statistic.
# This is important because the proposed statistic depends on r
# and on the weighting scheme.
# ============================================================

Bcal <- 20000


# ------------------------------------------------------------
# Critical value: classical Neumann
# ------------------------------------------------------------

calibrate_neumann <- function(n, B = Bcal, alpha = 0.05) {
  
  vals <- numeric(B)
  
  for (b in seq_len(B)) {
    
    x <- runif(n)
    
    vals[b] <- neumann_stat(x)
  }
  
  quantile(vals, probs = 1 - alpha, na.rm = TRUE)
}


# ------------------------------------------------------------
# Critical value: proposed spacing statistic
# ------------------------------------------------------------

calibrate_spacing <- function(
    n,
    r,
    weight = "Uniform",
    B = Bcal,
    alpha = 0.05) {
  
  vals <- numeric(B)
  
  for (b in seq_len(B)) {
    
    x <- runif(n)
    
    vals[b] <- spacing_neumann(
      x,
      r = r,
      weight = weight
    )
  }
  
  quantile(vals, probs = 1 - alpha, na.rm = TRUE)
}


# ------------------------------------------------------------
# Critical value: KS
# ------------------------------------------------------------

calibrate_ks <- function(n, B = Bcal, alpha = 0.05) {
  
  vals <- numeric(B)
  
  for (b in seq_len(B)) {
    
    x <- runif(n)
    
    vals[b] <- as.numeric(ks_stat(x))
  }
  
  quantile(vals, probs = 1 - alpha)
}


# ------------------------------------------------------------
# Critical value: CvM
# ------------------------------------------------------------

calibrate_cvm <- function(n, B = Bcal, alpha = 0.05) {
  
  vals <- numeric(B)
  
  for (b in seq_len(B)) {
    
    x <- runif(n)
    
    vals[b] <- cvm_stat(x)
  }
  
  quantile(vals, probs = 1 - alpha)
}


# ============================================================
# CONTIGUOUS / CLUSTERED ALTERNATIVE
# ============================================================

# A simple contiguous clustered alternative on (0,1).
#
# The perturbation is concentrated around the centre.
# c controls the strength of the departure.
# ============================================================

generate_clustered <- function(n, c = 1) {
  
  # Local perturbation around x = 0.5
  #
  # f_n(x) proportional to
  # 1 + c/sqrt(n) * h(x)
  
  h <- function(x) {
    exp(-((x - 0.5)^2) / (2 * 0.08^2))
  }
  
  # Rejection sampling
  out <- numeric(n)
  
  i <- 1
  
  # Upper bound for rejection probability
  grid <- seq(0, 1, length.out = 5000)
  
  hmax <- max(h(grid))
  
  M <- 1 + abs(c) / sqrt(n) * hmax
  
  while (i <= n) {
    
    y <- runif(1)
    
    fy <- 1 + c / sqrt(n) * h(y)
    
    if (runif(1) <= fy / M) {
      
      out[i] <- y
      
      i <- i + 1
    }
  }
  
  out
}


# ============================================================
# MIXTURE ALTERNATIVE
# ============================================================

# Bimodal Gaussian mixture.
# The sample is transformed through the normal CDF so that the
# resulting observations lie on (0,1).
# ============================================================

generate_mixture <- function(
    n,
    pi = 0.20,
    mu = 1.0) {
  
  z <- numeric(n)
  
  component <- runif(n) < pi
  
  z[!component] <- rnorm(
    sum(!component),
    mean = 0,
    sd = 1
  )
  
  z[component] <- rnorm(
    sum(component),
    mean = mu,
    sd = 1
  )
  
  # Probability integral transform
  pnorm(z)
}


# ============================================================
# MONTE CARLO STANDARD ERROR
# ============================================================

mcse <- function(p, N = Nrep) {
  
  sqrt(p * (1 - p) / N)
}


# ============================================================
# SIMULATION 1:
# EMPIRICAL SIZE
# ============================================================

cat("\n============================================\n")
cat("SIMULATION: EMPIRICAL SIZE\n")
cat("============================================\n")


size_results <- data.frame()


for (n in sample_sizes) {
  
  cat("n =", n, "\n")
  
  # ----------------------------------------------------------
  # Calibrate critical values
  # ----------------------------------------------------------
  
  cv_neu <- calibrate_neumann(
    n,
    B = Bcal,
    alpha = alpha
  )
  
  cv_ks <- calibrate_ks(
    n,
    B = Bcal,
    alpha = alpha
  )
  
  cv_cvm <- calibrate_cvm(
    n,
    B = Bcal,
    alpha = alpha
  )
  
  # Use r = 1 and uniform weighting for the basic proposed test
  cv_prop <- calibrate_spacing(
    n,
    r = 1,
    weight = "Uniform",
    B = Bcal,
    alpha = alpha
  )
  
  reject_prop <- logical(Nrep)
  reject_neu  <- logical(Nrep)
  reject_ks   <- logical(Nrep)
  reject_cvm  <- logical(Nrep)
  
  for (b in seq_len(Nrep)) {
    
    x <- runif(n)
    
    Tprop <- spacing_neumann(
      x,
      r = 1,
      weight = "Uniform"
    )
    
    Tneu <- neumann_stat(x)
    
    Tks <- as.numeric(ks_stat(x))
    
    Tcvm <- cvm_stat(x)
    
    reject_prop[b] <- Tprop > cv_prop
    reject_neu[b]  <- Tneu  > cv_neu
    reject_ks[b]   <- Tks   > cv_ks
    reject_cvm[b]  <- Tcvm  > cv_cvm
  }
  
  p_prop <- mean(reject_prop, na.rm = TRUE)
  p_neu  <- mean(reject_neu, na.rm = TRUE)
  p_ks   <- mean(reject_ks, na.rm = TRUE)
  p_cvm  <- mean(reject_cvm, na.rm = TRUE)
  
  size_results <- rbind(
    size_results,
    data.frame(
      n = n,
      Proposed = 100 * p_prop,
      Proposed_MCSE = 100 * mcse(p_prop),
      Neumann = 100 * p_neu,
      Neumann_MCSE = 100 * mcse(p_neu),
      KS = 100 * p_ks,
      KS_MCSE = 100 * mcse(p_ks),
      CvM = 100 * p_cvm,
      CvM_MCSE = 100 * mcse(p_cvm)
    )
  )
}


print(size_results)


# ============================================================
# SIMULATION 2:
# POWER UNDER CONTIGUOUS CLUSTERED ALTERNATIVES
# ============================================================

cat("\n============================================\n")
cat("SIMULATION: CONTIGUOUS CLUSTERED ALTERNATIVE\n")
cat("============================================\n")


power_local <- data.frame()

c_local <- 1


for (n in sample_sizes) {
  
  cat("n =", n, "\n")
  
  # ----------------------------------------------------------
  # Critical values from the null
  # ----------------------------------------------------------
  
  cv_neu <- calibrate_neumann(
    n,
    B = Bcal,
    alpha = alpha
  )
  
  cv_ks <- calibrate_ks(
    n,
    B = Bcal,
    alpha = alpha
  )
  
  cv_cvm <- calibrate_cvm(
    n,
    B = Bcal,
    alpha = alpha
  )
  
  cv_prop <- calibrate_spacing(
    n,
    r = 1,
    weight = "Uniform",
    B = Bcal,
    alpha = alpha
  )
  
  reject_prop <- logical(Nrep)
  reject_neu  <- logical(Nrep)
  reject_ks   <- logical(Nrep)
  reject_cvm  <- logical(Nrep)
  
  for (b in seq_len(Nrep)) {
    
    x <- generate_clustered(
      n,
      c = c_local
    )
    
    Tprop <- spacing_neumann(
      x,
      r = 1,
      weight = "Uniform"
    )
    
    Tneu <- neumann_stat(x)
    
    Tks <- as.numeric(ks_stat(x))
    
    Tcvm <- cvm_stat(x)
    
    reject_prop[b] <- Tprop > cv_prop
    reject_neu[b]  <- Tneu  > cv_neu
    reject_ks[b]   <- Tks   > cv_ks
    reject_cvm[b]  <- Tcvm  > cv_cvm
  }
  
  p_prop <- mean(reject_prop, na.rm = TRUE)
  p_neu  <- mean(reject_neu, na.rm = TRUE)
  p_ks   <- mean(reject_ks, na.rm = TRUE)
  p_cvm  <- mean(reject_cvm, na.rm = TRUE)
  
  power_local <- rbind(
    power_local,
    data.frame(
      n = n,
      Proposed = 100 * p_prop,
      Proposed_MCSE = 100 * mcse(p_prop),
      Neumann = 100 * p_neu,
      Neumann_MCSE = 100 * mcse(p_neu),
      KS = 100 * p_ks,
      KS_MCSE = 100 * mcse(p_ks),
      CvM = 100 * p_cvm,
      CvM_MCSE = 100 * mcse(p_cvm)
    )
  )
}


print(power_local)


# ============================================================
# SIMULATION 3:
# POWER UNDER MIXTURE ALTERNATIVES
# ============================================================

cat("\n============================================\n")
cat("SIMULATION: MIXTURE ALTERNATIVE\n")
cat("============================================\n")


power_mixture <- data.frame()


for (n in sample_sizes) {
  
  cat("n =", n, "\n")
  
  # Null critical values
  
  cv_neu <- calibrate_neumann(
    n,
    B = Bcal,
    alpha = alpha
  )
  
  cv_ks <- calibrate_ks(
    n,
    B = Bcal,
    alpha = alpha
  )
  
  cv_cvm <- calibrate_cvm(
    n,
    B = Bcal,
    alpha = alpha
  )
  
  cv_prop <- calibrate_spacing(
    n,
    r = 1,
    weight = "Uniform",
    B = Bcal,
    alpha = alpha
  )
  
  reject_prop <- logical(Nrep)
  reject_neu  <- logical(Nrep)
  reject_ks   <- logical(Nrep)
  reject_cvm  <- logical(Nrep)
  
  for (b in seq_len(Nrep)) {
    
    x <- generate_mixture(
      n,
      pi = 0.20,
      mu = 1.0
    )
    
    Tprop <- spacing_neumann(
      x,
      r = 1,
      weight = "Uniform"
    )
    
    Tneu <- neumann_stat(x)
    
    Tks <- as.numeric(ks_stat(x))
    
    Tcvm <- cvm_stat(x)
    
    reject_prop[b] <- Tprop > cv_prop
    reject_neu[b]  <- Tneu  > cv_neu
    reject_ks[b]   <- Tks   > cv_ks
    reject_cvm[b]  <- Tcvm  > cv_cvm
  }
  
  p_prop <- mean(reject_prop, na.rm = TRUE)
  p_neu  <- mean(reject_neu, na.rm = TRUE)
  p_ks   <- mean(reject_ks, na.rm = TRUE)
  p_cvm  <- mean(reject_cvm, na.rm = TRUE)
  
  power_mixture <- rbind(
    power_mixture,
    data.frame(
      n = n,
      Proposed = 100 * p_prop,
      Proposed_MCSE = 100 * mcse(p_prop),
      Neumann = 100 * p_neu,
      Neumann_MCSE = 100 * mcse(p_neu),
      KS = 100 * p_ks,
      KS_MCSE = 100 * mcse(p_ks),
      CvM = 100 * p_cvm,
      CvM_MCSE = 100 * mcse(p_cvm)
    )
  )
}


print(power_mixture)


# ============================================================
# SIMULATION 4:
# EFFECT OF HIGHER-ORDER SPACING
# ============================================================

cat("\n============================================\n")
cat("SIMULATION: HIGHER-ORDER SPACINGS\n")
cat("============================================\n")


r_values <- c(1, 2, 3, 5, 8)


r_results <- data.frame()


for (n in sample_sizes) {
  
  cat("\nn =", n, "\n")
  
  for (r in r_values) {
    
    if (r >= n - 1) next
    
    cat("  r =", r, "\n")
    
    # --------------------------------------------------------
    # Null calibration
    # --------------------------------------------------------
    
    cv <- calibrate_spacing(
      n,
      r = r,
      weight = "Uniform",
      B = Bcal,
      alpha = alpha
    )
    
    reject <- logical(Nrep)
    
    for (b in seq_len(Nrep)) {
      
      # Mixture alternative
      x <- generate_mixture(
        n,
        pi = 0.20,
        mu = 1.0
      )
      
      T <- spacing_neumann(
        x,
        r = r,
        weight = "Uniform"
      )
      
      reject[b] <- T > cv
    }
    
    p <- mean(reject, na.rm = TRUE)
    
    r_results <- rbind(
      r_results,
      data.frame(
        n = n,
        r = r,
        Power = 100 * p,
        MCSE = 100 * mcse(p)
      )
    )
  }
}


print(r_results)


# ============================================================
# SIMULATION 5:
# EFFECT OF WEIGHTING SCHEMES
# ============================================================

cat("\n============================================\n")
cat("SIMULATION: WEIGHTING SCHEMES\n")
cat("============================================\n")


weights <- c(
  "Uniform",
  "Central",
  "Tail"
)


weight_results <- data.frame()


# Use r = 1 for direct weighting comparison

r_weight <- 1


for (n in sample_sizes) {
  
  cat("\nn =", n, "\n")
  
  for (wt in weights) {
    
    cat("  Weight =", wt, "\n")
    
    # --------------------------------------------------------
    # Null critical value
    # --------------------------------------------------------
    
    cv <- calibrate_spacing(
      n,
      r = r_weight,
      weight = wt,
      B = Bcal,
      alpha = alpha
    )
    
    reject <- logical(Nrep)
    
    for (b in seq_len(Nrep)) {
      
      x <- generate_clustered(
        n,
        c = c_local
      )
      
      T <- spacing_neumann(
        x,
        r = r_weight,
        weight = wt
      )
      
      reject[b] <- T > cv
    }
    
    p <- mean(reject, na.rm = TRUE)
    
    weight_results <- rbind(
      weight_results,
      data.frame(
        n = n,
        Weight = wt,
        Power = 100 * p,
        MCSE = 100 * mcse(p)
      )
    )
  }
}


print(weight_results)


# ============================================================
# SAVE RESULTS
# ============================================================

write.csv(
  size_results,
  "simulation_empirical_size.csv",
  row.names = FALSE
)

write.csv(
  power_local,
  "simulation_contiguous_power.csv",
  row.names = FALSE
)

write.csv(
  power_mixture,
  "simulation_mixture_power.csv",
  row.names = FALSE
)

write.csv(
  r_results,
  "simulation_higher_order_spacing.csv",
  row.names = FALSE
)

write.csv(
  weight_results,
  "simulation_weighting_schemes.csv",
  row.names = FALSE
)


# ============================================================
# DISPLAY FINAL RESULTS
# ============================================================

cat("\n\n============================================\n")
cat("FINAL SIZE RESULTS\n")
cat("============================================\n")

print(size_results)


cat("\n\n============================================\n")
cat("FINAL CONTIGUOUS POWER RESULTS\n")
cat("============================================\n")

print(power_local)


cat("\n\n============================================\n")
cat("FINAL MIXTURE POWER RESULTS\n")
cat("============================================\n")

print(power_mixture)


cat("\n\n============================================\n")
cat("FINAL HIGHER-ORDER RESULTS\n")
cat("============================================\n")

print(r_results)


cat("\n\n============================================\n")
cat("FINAL WEIGHTING RESULTS\n")
cat("============================================\n")

print(weight_results)