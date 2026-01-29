# Load required libraries
library(MASS)      # for negative binomial if needed
library(broom)     # tidy summaries (optional but clean)
library(dplyr)
library(faraway)


############################################################
############################################################

#  Part 1: GLM, Modeling Micro-Credit Uptake in Rural Communities 

############################################################
############################################################


# Load data
data <- read.csv("microcredit_impact.csv")

##################   Model: M0   ######################################

m0 <- glm(
  applications ~ 1 + offset(log(meetings_held)),
  family = poisson(link = "log"),
  data = data
)

summary(m0)

##################   Model: M1   ######################################

m1 <- glm(
  applications ~ literacy_rate + offset(log(meetings_held)),
  family = poisson(link = "log"),
  data = data
)

summary(m1)

##################   Model: M2   ######################################


m2 <- glm(
  applications ~ literacy_rate + infrastructure + offset(log(meetings_held)),
  family = poisson(link = "log"),
  data = data
)

summary(m2)

##################   Model: M3   ######################################


m3 <- glm(
  applications ~ literacy_rate + infrastructure + avg_income +
    offset(log(meetings_held)),
  family = poisson(link = "log"),
  data = data
)

summary(m3)

##################   Summary table   ######################################

model_table <- tibble(
  Model = c("M0", "M1", "M2", "M3"),
  
  # Information criterion
  AIC = c(AIC(m0), AIC(m1), AIC(m2), AIC(m3)),
  
  # Residual deviance and degrees of freedom
  Residual_Deviance = c(
    deviance(m0),
    deviance(m1),
    deviance(m2),
    deviance(m3)
  ),
  
  Residual_DF = c(
    df.residual(m0),
    df.residual(m1),
    df.residual(m2),
    df.residual(m3)
  ),
  
  # Deviance-based GOF
  Deviance_DF_Ratio = Residual_Deviance / Residual_DF,
  
  # Likelihood ratio test vs previous model
  Delta_Deviance = c(
    NA,
    deviance(m0) - deviance(m1),
    deviance(m1) - deviance(m2),
    deviance(m2) - deviance(m3)
  ),
  
  LR_p_value = c(
    NA,
    pchisq(deviance(m0) - deviance(m1),
           df.residual(m0) - df.residual(m1),
           lower.tail = FALSE),
    pchisq(deviance(m1) - deviance(m2),
           df.residual(m1) - df.residual(m2),
           lower.tail = FALSE),
    pchisq(deviance(m2) - deviance(m3),
           df.residual(m2) - df.residual(m3),
           lower.tail = FALSE)
  ),
  
  # Pearson Chi-square / DF (overdispersion diagnostic)
  Pearson_ChiSq_DF = c(
    sum(residuals(m0, type = "pearson")^2) / df.residual(m0),
    sum(residuals(m1, type = "pearson")^2) / df.residual(m1),
    sum(residuals(m2, type = "pearson")^2) / df.residual(m2),
    sum(residuals(m3, type = "pearson")^2) / df.residual(m3)
  )
)

model_table




############################################################
############################################################

# Part 2: Computational Statistics 

############################################################
############################################################

##################   Poisson IRLS   ######################################

poisson_irls <- function(X, y, tol = 1e-6, maxit = 100, quasi = FALSE) {
  # ===============================================================
  # IRLS for Poisson Regression
  # Input:
  #   X: Design matrix (n x p), including intercept
  #   y: Response vector of observed counts
  #   tol: Convergence tolerance
  #   maxit: Maximum number of iterations
  #   quasi: If TRUE, compute quasi-Poisson standard errors (overdispersion)
  # Output:
  #   beta: Estimated coefficients
  #   se: Standard errors of coefficients
  #   iterations: Number of IRLS iterations performed
  # ===============================================================
  
  n <- nrow(X)
  p <- ncol(X)
  
  # 1. Initialization
  k <- 0                     # iteration counter
  beta <- rep(0, p)          # starting value for beta (all zeros)
  
  # 2. Iteration Loop: repeat until convergence
  for (k in 1:maxit) {
    
    # ---------------------------------------------------------------
    # a) Calculate the Linear Predictor: eta^(k) = X %*% beta^(k)
    eta <- as.vector(X %*% beta)
    
    # Cap eta to avoid numerical overflow in exp()
    eta[eta > 20] <- 20
    eta[eta < -20] <- -20
    
    # ---------------------------------------------------------------
    # b) Calculate the Fitted Values: mu^(k) = exp(eta^(k))
    mu <- exp(eta)
    
    # Ensure mu is not too small or too large to avoid zero or Inf weights
    mu[mu < 1e-8] <- 1e-8
    mu[mu > 1e8] <- 1e8
    
    # ---------------------------------------------------------------
    # c) Construct the Weight Matrix W^(k)
    # w_ii = mu_i^(k)
    # We store only the sqrt(mu) vector to scale rows for QR decomposition
    W_sqrt <- sqrt(mu)
    
    # ---------------------------------------------------------------
    # d) Calculate the Adjusted Dependent Variable (z^(k))
    z <- eta + (y - mu) / mu
    
    # ---------------------------------------------------------------
    # e) Update beta^(k+1) using Weighted Least Squares via QR decomposition
    # Scale each row of X and z by sqrt(weight)
    WX <- sweep(X, 1, W_sqrt, FUN = "*")  # weighted X
    Wz <- z * W_sqrt                       # weighted z
    
    qr_fit <- qr(WX)
    beta_new <- qr.coef(qr_fit, Wz)
    
    # Check convergence
    if (max(abs(beta_new - beta), na.rm = TRUE) < tol) break
    
    beta <- beta_new
  }
  
  # 3. Output: final beta and standard errors
  # ---------------------------------------------------------------
  
  # Compute final weights
  eta <- X %*% beta
  eta[eta > 20] <- 20
  eta[eta < -20] <- -20
  mu <- exp(eta)
  mu[mu < 1e-8] <- 1e-8
  mu[mu > 1e8] <- 1e8
  W_sqrt <- sqrt(mu)
  
  # Construct weighted design matrix for SE computation
  WX <- sweep(X, 1, W_sqrt, FUN = "*")
  qr_fit <- qr(WX)
  R <- qr.R(qr_fit)
  
  # Variance-covariance matrix
  vcov <- tryCatch(solve(t(R) %*% R), error = function(e) MASS::ginv(t(R) %*% R))
  
  # Quasi-Poisson adjustment for overdispersion, if requested
  if (quasi) {
    phi_hat <- sum((y - mu)^2 / mu) / (n - p)
    vcov <- phi_hat * vcov
  }
  
  # Standard errors
  se <- sqrt(diag(vcov))
  
  # Return results
  list(beta = beta, se = se, iterations = k)
}



################## Comparison with GLM on Galapagos dataset    ####################################


data(gala)

# Design matrix including intercept
X <- model.matrix(Species ~ ., data = gala)
y <- gala$Species

# Poisson IRLS
fit_poisson <- poisson_irls(X, y)

# Poisson IRLS with quasi-Poisson adjustment
fit_quasi <- poisson_irls(X, y, quasi = TRUE)

# Fit Glm
fit_glm <- glm(Species ~ ., family = poisson, data = gala)
fit_glm_qp <- glm(Species ~ ., family = quasipoisson, data = gala)

# Collect coefficients and SEs
coef_manual_poisson <- fit_poisson$beta
se_manual_poisson <- fit_poisson$se

coef_manual_quasi <- fit_quasi$beta
se_manual_quasi <- fit_quasi$se

coef_glm_poisson <- coef(fit_glm)
se_glm_poisson <- sqrt(diag(vcov(fit_glm)))

coef_glm_quasi <- coef(fit_glm_qp)
se_glm_quasi <- sqrt(diag(vcov(fit_glm_qp)))

# Create comparison table
variables <- names(coef_glm_poisson)  # variable names from glm
comparison_table <- data.frame(
  Variable = variables,
  Manual_Poisson = coef_manual_poisson,
  GLM_Poisson = coef_glm_poisson,
  SE_Manual_Poisson = se_manual_poisson,
  SE_GLM_Poisson = se_glm_poisson,
  Manual_Quasi = coef_manual_quasi,
  GLM_Quasi = coef_glm_quasi,
  SE_Manual_Quasi = se_manual_quasi,
  SE_GLM_Quasi = se_glm_quasi
)

# Display table
print(comparison_table, digits = 6)


############################################################
############################################################

# Part 3: GLM and Computational Statistics 

############################################################
############################################################


# GLM Model 

#Fit Poisson GLM
model <- glm(Species ~ Area + Elevation + Nearest + Scruz + Adjacent,
             family = poisson(link = "log"),
             data = gala)


#################### Non-parametric Bootstrap #########

#Setup
set.seed(123)
B <- 1000
n <- nrow(gala)

# Storage for bootstrap coefficients
beta_area <- numeric(B)

#Resampling and ftting
for (b in 1:B) {
  
  # Resample rows with replacement
  idx <- sample(1:n, size = n, replace = TRUE)
  gala_boot <- gala[idx, ]
  
  # Fit Poisson GLM to bootstrap sample
  boot_model <- glm(Species ~ Area + Elevation + Nearest + Scruz + Adjacent,
                    family = poisson(link = "log"),
                    data = gala_boot)
  
  # Store coefficient of interest (Area)
  beta_area[b] <- coef(boot_model)["Area"]
}


################ Comparison ###################

# Info from bootstrap

# Mean and median
boot_mean   <- mean(beta_area)
boot_median <- median(beta_area)

# Mode estimation using kernel density
dens <- density(beta_area)
boot_mode <- dens$x[which.max(dens$y)]

# Print results
c(Mean = boot_mean,
  Median = boot_median,
  Mode = boot_mode)

# Bootstrap Percentile CI
boot_ci <- quantile(beta_area, probs = c(0.025, 0.975))
boot_ci

# Info from GLM model
original_estimate <- coef(model)["Area"]
c(original_estimate)

# Built-in asymptotic CI
asym_ci <- confint(model, parm = "Area")
asym_ci


#Side by side comparison
ci_comparison <- rbind(
  Bootstrap = boot_ci,
  Asymptotic = asym_ci
)

#plot
ci_comparison

hist(beta_area, breaks = 40, col = "lightgray",
     main = "Bootstrap Distribution of Area Coefficient",
     xlab = "Coefficient value")

abline(v = asym_ci, col = "red", lwd = 2)
abline(v = boot_ci, col = "blue", lwd = 2, lty = 2)

legend("topright",
       legend = c("Asymptotic CI", "Bootstrap CI"),
       col = c("red", "blue"),
       lwd = 2,
       lty = c(1, 2))


