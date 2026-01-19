############################################################
## Shared utilities: KL-ball robustification by entropic tilting
############################################################

# Robust risk over a KL ball around a discrete baseline distribution
# losses: numeric vector L_s
# w: baseline weights w_s (will be normalised)
# C: KL radius
# returns sup_{q: KL(q||w) <= C} sum_s q_s * L_s
robust_kl_risk_entropic <- function(losses, w, C, lambda_max = 100) {
  losses <- as.numeric(losses)
  w <- as.numeric(w)
  w <- w / sum(w)
  base_risk <- sum(w * losses)
  
  if (C <= 0) return(base_risk)
  if (max(losses) - min(losses) < 1e-12) return(base_risk)  # constant loss
  
  # Dual objective ψ(λ) = (C + log Σ w_s exp(λ L_s)) / λ
  psi <- function(lambda) {
    if (lambda <= 0) return(Inf)
    z <- lambda * losses
    z_max <- max(z)
    logZ <- z_max + log(sum(w * exp(z - z_max)))  # log-sum-exp
    (C + logZ) / lambda
  }
  
  opt <- optimize(psi, interval = c(1e-8, lambda_max))
  lambda_star <- opt$minimum
  
  # Least-favourable (tilted) weights: q_s ∝ w_s exp(λ* L_s)
  z <- lambda_star * losses
  z_max <- max(z)
  w_tilt_unnorm <- w * exp(z - z_max)
  w_tilt <- w_tilt_unnorm / sum(w_tilt_unnorm)
  
  sum(w_tilt * losses)
}




############################################################
## Experiment A: ER vs 2-block SBM (known labels)
############################################################

## Re-use your make_ER_SBM_setup(), simulate_ER_graph(),
## simulate_SBM_graph(), compute_loglik_ER(), compute_loglik_SBM()
## and single_rep_ER_SBM() as they are.

## Small-C normalisation for 0–1 loss
run_ER_SBM_experiment_A <- function(
    n, c, lambda,
    n_rep   = 1000,
    C_grid  = 10^seq(-6, -3, length.out = 12),
    model_sampling = "mix"
) {
  setup <- make_ER_SBM_setup(n)
  
  e0_vec    <- numeric(n_rep)
  e_rob_mat <- matrix(NA_real_, nrow = n_rep, ncol = length(C_grid))
  
  for (r in seq_len(n_rep)) {
    rep_res <- single_rep_ER_SBM(
      n = n, c = c, lambda = lambda,
      setup = setup,
      model_sampling = model_sampling,
      C_grid = C_grid
    )
    e0_vec[r]      <- rep_res$e0
    e_rob_mat[r, ] <- rep_res$e_rob
  }
  
  R0   <- mean(e0_vec)
  Rrob <- colMeans(e_rob_mat)
  
  ## Small-radius expansion for Bernoulli loss:
  ## e_rob(G) - e0(G) ≈ sqrt{2 e0(G) (1 - e0(G))} * sqrt(C)
  diff_mat <- sweep(e_rob_mat, 1, e0_vec, FUN = "-")
  denom_row <- sqrt(2 * e0_vec * (1 - e0_vec))
  norm_mat  <- sweep(diff_mat, 1, denom_row, FUN = "/")
  norm_mat  <- sweep(norm_mat, 2, sqrt(C_grid), FUN = "/")
  norm_curve <- colMeans(norm_mat)
  
  list(
    n = n, c = c, lambda = lambda,
    C_grid = C_grid,
    R0 = R0, Rrob = Rrob,
    e0_all = e0_vec,
    e_rob_all = e_rob_mat,
    norm_curve = norm_curve
  )
}

## Example: n = 400, c = 3, lambda = 0.4 (near detection threshold)
set.seed(1)
C_grid_A <- 10^seq(-6, -3, length.out = 12)

expA_400_lam04 <- run_ER_SBM_experiment_A(
  n = 400, c = 3, lambda = 0.4,
  n_rep  = 1000,
  C_grid = C_grid_A,
  model_sampling = "mix"
)

## Plot robustness sensitivity curve
plot(
  sqrt(expA_400_lam04$C_grid), expA_400_lam04$norm_curve,
  type = "b", pch = 19,
  xlab = expression(sqrt(C)),
  ylab = expression(
    (e[rob] - e[0]) /
      (sqrt(2 * e[0] * (1 - e[0])) * sqrt(C))
  ),
  main = "Synthetic experiment A: ER vs SBM (n=400, c=3, lambda=0.4)"
)
abline(h = 1, lty = 2, col = "red")


############################################################
## Experiment B: configuration model 
############################################################

## Keep your rgamma_trunc() and single_rep_config_Poisson() as they are.

run_config_experiment <- function(
    n, Delta,
    n_rep   = 200,
    C_grid  = 10^seq(-4, -1, length.out = 10),
    alpha0  = 1, beta0 = 1,
    S_post  = 2000
) {
  rho0_vec    <- numeric(n_rep)
  rho_rob_mat <- matrix(NA_real_, nrow = n_rep, ncol = length(C_grid))
  
  for (r in seq_len(n_rep)) {
    rep_res <- single_rep_config_Poisson(
      n = n,
      Delta = Delta,
      alpha0 = alpha0,
      beta0  = beta0,
      S_post = S_post,
      C_grid = C_grid
    )
    rho0_vec[r]     <- rep_res$rho0
    rho_rob_mat[r,] <- rep_res$rho_rob
  }
  
  R0   <- mean(rho0_vec)
  Rrob <- colMeans(rho_rob_mat)
  
  ## Normalised robustness sensitivity curve:
  ## (rho_rob - rho_0)/(rho_0 sqrt(C)) → 2
  diff_mat <- sweep(rho_rob_mat, 1, rho0_vec, FUN = "-")
  denom    <- outer(rho0_vec, sqrt(C_grid))
  norm_mat <- diff_mat / denom
  norm_curve <- colMeans(norm_mat)
  
  list(
    n = n,
    Delta = Delta,
    C_grid = C_grid,
    R0 = R0,
    Rrob = Rrob,
    norm_curve = norm_curve,
    rho0_all = rho0_vec,
    rho_rob_all = rho_rob_mat
  )
}

## Example: larger n
set.seed(2)
C_grid_B <- 10^seq(-4, -1, length.out = 12)

expB_n5000_Delta02 <- run_config_experiment(
  n      = 5000,
  Delta  = 0.2,
  n_rep  = 200,
  C_grid = C_grid_B,
  S_post = 2000
)

plot(
  sqrt(expB_n5000_Delta02$C_grid), expB_n5000_Delta02$norm_curve,
  type = "b", pch = 19,
  xlab = expression(sqrt(C)),
  ylab = expression((rho[rob] - rho[0]) / (rho[0] * sqrt(C))),
  main = "Experiment B: configuration model (n=5000, Delta=0.2)"
)
abline(h = 2, lty = 2, col = "red")




run_config_grid_and_slopes <- function(
    n,
    Delta_vec = c(0.4, 0.3, 0.25, 0.2, 0.17, 0.15),
    n_rep     = 200,
    C_grid    = 10^seq(-4, -1, length.out = 10),
    alpha0    = 1, beta0 = 1,
    S_post    = 2000,
    C_fixed   = 1e-3
) {
  res_list <- vector("list", length(Delta_vec))
  names(res_list) <- paste0("Delta_", Delta_vec)
  
  for (i in seq_along(Delta_vec)) {
    cat("Running config experiment for Delta =", Delta_vec[i], "\n")
    res_list[[i]] <- run_config_experiment(
      n      = n,
      Delta  = Delta_vec[i],
      n_rep  = n_rep,
      C_grid = C_grid,
      alpha0 = alpha0,
      beta0  = beta0,
      S_post = S_post
    )
  }
  
  R0_vec    <- sapply(res_list, `[[`, "R0")
  Delta_vec <- sapply(res_list, `[[`, "Delta")
  n_val     <- res_list[[1]]$n
  C_grid0   <- res_list[[1]]$C_grid
  
  ## Choose C closest to C_fixed
  C_idx  <- which.min(abs(C_grid0 - C_fixed))
  C_used <- C_grid0[C_idx]
  Rrob_C <- sapply(res_list, function(res) res$Rrob[C_idx])
  
  x  <- -log(Delta_vec)
  y1 <- log(n_val * R0_vec)
  y2 <- log(n_val * (Rrob_C - R0_vec) / sqrt(C_used))
  
  slope_base   <- coef(lm(y1 ~ x))[2]
  slope_robust <- coef(lm(y2 ~ x))[2]
  
  list(
    n = n_val,
    Delta_vec = Delta_vec,
    C_used = C_used,
    R0_vec = R0_vec,
    Rrob_C_vec = Rrob_C,
    x = x,
    y_base = y1,
    y_robust = y2,
    slope_baseline = slope_base,
    slope_robust = slope_robust,
    results = res_list
  )
}

## Run grid and plot
set.seed(3)
slope_res <- run_config_grid_and_slopes(
  n         = 5000,  ## <-- main change from your script
  Delta_vec = c(0.4, 0.3, 0.25, 0.2, 0.17, 0.15),
  n_rep     = 200,
  C_grid    = C_grid_B,
  S_post    = 2000,
  C_fixed   = 1e-3
)

cat("Estimated slope for log(n * rho_0) vs -log Delta  ≈",
    slope_res$slope_baseline, "\n")
cat("Estimated slope for robust excess vs -log Delta ≈",
    slope_res$slope_robust, "\n")

par(mfrow = c(1, 2))
with(slope_res, {
  plot(
    x, y_base,
    xlab = "-log Delta", ylab = "log(n * rho_0)",
    main = "Baseline MSE scaling"
  )
  abline(lm(y_base ~ x), col = "red", lty = 2)
  
  plot(
    x, y_robust,
    xlab = "-log Delta",
    ylab = "log(n * (rho[rob] - rho[0]) / sqrt(C_used))",
    main = "Robust excess scaling"
  )
  abline(lm(y_robust ~ x), col = "red", lty = 2)
})
par(mfrow = c(1, 1))





## --- R: rerun experiments with final settings ---

set.seed(1)
C_grid_A <- 10^seq(-4, -2, length.out = 12)

expA_400_04 <- run_ER_SBM_experiment(
  n       = 400,
  c       = 3,
  lambda  = 0.4,
  n_rep   = 1000,
  C_grid  = C_grid_A,
  model_sampling = "mix"
)

set.seed(2)
C_grid_B <- 10^seq(-4, -1, length.out = 12)

expB_n5000_Delta02 <- run_config_experiment(
  n       = 5000,
  Delta   = 0.2,
  n_rep   = 200,
  C_grid  = C_grid_B,
  S_post  = 2000
)

set.seed(3)
slope_res <- run_config_grid_and_slopes(
  n         = 5000,
  Delta_vec = c(0.40, 0.30, 0.25, 0.20, 0.17, 0.15),
  n_rep     = 200,
  C_grid    = C_grid_B,
  S_post    = 2000,
  C_fixed   = 1e-3
)




pdf("exp_small_radius.pdf", width = 7, height = 3.5)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1) + 0.1)

## (a) ER vs SBM
R0A <- expA_400_04$R0
plot(
  sqrt(expA_400_04$C_grid), expA_400_04$norm_curve,
  type = "b", pch = 19,
  xlab = expression(sqrt(C)),
  ylab = expression(
    (R[rob] - R[0]) /
      (sqrt(2 * R[0] * (1 - R[0])) * sqrt(C))
  ),
  main = "(a) ER vs SBM (n = 400, c = 3, λ = 0.4)"
)
abline(h = 1, lty = 2)

## (b) configuration model
plot(
  sqrt(expB_n5000_Delta02$C_grid), expB_n5000_Delta02$norm_curve,
  type = "b", pch = 19,
  xlab = expression(sqrt(C)),
  ylab = expression(
    (rho[rob] - rho[0]) / (rho[0] * sqrt(C))
  ),
  main = "(b) Config. model (n = 5000, Δ = 0.2)"
)
abline(h = 2, lty = 2)

dev.off()


pdf("exp_config_slopes.pdf", width = 7, height = 3.5)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1) + 0.1)

with(slope_res, {
  ## Baseline MSE scaling
  plot(
    x, y_base, pch = 19,
    xlab = "-log Delta", ylab = "log(n * rho_0)",
    main = "Baseline MSE scaling"
  )
  abline(lm(y_base ~ x), lty = 2)
  
  ## Robust excess scaling
  plot(
    x, y_robust, pch = 19,
    xlab = "-log Delta",
    ylab = expression(
      log(n * (rho[rob] - rho[0]) / sqrt(C[used]))
    ),
    main = "Robust excess scaling"
  )
  abline(lm(y_robust ~ x), lty = 2)
})

dev.off()


tab_slopes <- data.frame(
  n                = slope_res$n,
  slope_baseline   = slope_res$slope_baseline,
  slope_robust     = slope_res$slope_robust
)

tab_slopes
#     n slope_baseline slope_robust
# 1 5000      4.48703     4.652716

# LaTeX via knitr::kable (or xtable)
knitr::kable(
  tab_slopes,
  digits = c(0, 2, 2),
  caption = "Estimated scaling exponents in the configuration-model experiment."
)










####### RUN ######


## --- R: rerun experiments with final settings ---

set.seed(1)
C_grid_A <- 10^seq(-4, -2, length.out = 12)

expA_400_04 <- run_ER_SBM_experiment(
  n       = 400,
  c       = 3,
  lambda  = 0.4,
  n_rep   = 1000,
  C_grid  = C_grid_A,
  model_sampling = "mix"
)

set.seed(2)
C_grid_B <- 10^seq(-4, -1, length.out = 12)

expB_n5000_Delta02 <- run_config_experiment(
  n       = 5000,
  Delta   = 0.2,
  n_rep   = 200,
  C_grid  = C_grid_B,
  S_post  = 2000
)

set.seed(3)
slope_res <- run_config_grid_and_slopes(
  n         = 5000,
  Delta_vec = c(0.40, 0.30, 0.25, 0.20, 0.17, 0.15),
  n_rep     = 200,
  C_grid    = C_grid_B,
  S_post    = 2000,
  C_fixed   = 1e-3
)

par(mfrow = c(1, 2), mar = c(4, 4, 2, 1) + 0.1)

## (a) ER vs SBM
R0A <- expA_400_04$R0
plot(
  sqrt(expA_400_04$C_grid), expA_400_04$norm_curve,
  type = "b", pch = 19,
  xlab = expression(sqrt(C)),
  ylab = expression(
    (R[rob] - R[0]) /
      (sqrt(2 * R[0] * (1 - R[0])) * sqrt(C))
  ),
  main = expression("(a) ER vs SBM (n = 400, c = 3, " * lambda * " = 0.4)")
)
abline(h = 1, lty = 2)

## (b) configuration model
plot(
  sqrt(expB_n5000_Delta02$C_grid), expB_n5000_Delta02$norm_curve,
  type = "b", pch = 19,
  xlab = expression(sqrt(C)),
  ylab = expression((rho[rob] - rho[0]) / (rho[0] * sqrt(C))),
  main = expression("(b) Config. model (n = 5000, " * Delta * " = 0.2)")
)
abline(h = 2, lty = 2)





par(mfrow = c(1, 2), mar = c(4, 4, 2, 1) + 0.1)

with(slope_res, {
  ## Baseline MSE scaling
  plot(
    x, y_base, pch = 19,
    xlab = "-log Delta", ylab = "log(n * rho_0)",
    main = "Baseline MSE scaling"
  )
  abline(lm(y_base ~ x), lty = 2)
  
  ## Robust excess scaling
  plot(
    x, y_robust, pch = 19,
    xlab = "-log Delta",
    ylab = expression(
      log(n * (rho[rob] - rho[0]) / sqrt(C[used]))
    ),
    main = "Robust excess scaling"
  )
  abline(lm(y_robust ~ x), lty = 2)
})




