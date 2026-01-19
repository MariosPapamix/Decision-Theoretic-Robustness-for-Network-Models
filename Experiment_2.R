## ======================================================
## Experiment 2: Robust functionals & SBM-vs-SBM model selection
## ======================================================

library(igraph)

## ---------- 1. Network functionals R(theta) = (C, L, S, lambda1) ----------

compute_network_functionals <- function(A) {
  A <- as.matrix(A)
  # symmetrise and zero diagonal
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) <- 0
  
  g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
  
  C <- igraph::transitivity(g, type = "global")
  L <- igraph::mean_distance(g, directed = FALSE, unconnected = TRUE)
  
  n <- igraph::vcount(g)
  m <- igraph::ecount(g)
  p_hat <- if (n > 1) 2 * m / (n * (n - 1)) else NA_real_
  C_er  <- p_hat
  
  if (!is.na(p_hat) && p_hat > 0 && n * p_hat > 1) {
    L_er <- log(n) / log(n * p_hat)
  } else {
    L_er <- NA_real_
  }
  
  S <- if (!is.na(C) && !is.na(C_er) && !is.na(L) && !is.na(L_er) &&
           C_er > 0 && L_er > 0) {
    (C / C_er) / (L / L_er)
  } else {
    NA_real_
  }
  
  lambda1 <- if (nrow(A) >= 2) {
    max(Re(eigen(A, only.values = TRUE)$values))
  } else {
    NA_real_
  }
  
  list(C = C, L = L, S = S, lambda1 = lambda1)
}

## ---------- 2. Simple SBM fit (spectral clustering + MLE) ----------

fit_sbm_spectral <- function(A, K, nstart = 10) {
  A <- as.matrix(A)
  n <- nrow(A)
  if (n != ncol(A)) stop("A must be square")
  
  # symmetrise and remove self-loops
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) <- 0
  
  # adjacency spectral embedding
  ev  <- eigen(A, symmetric = TRUE)
  ord <- order(abs(ev$values), decreasing = TRUE)
  U   <- ev$vectors[, ord[seq_len(K)], drop = FALSE]
  
  km <- kmeans(U, centers = K, nstart = nstart)
  z  <- km$cluster
  
  # log-likelihood under SBM with block MLEs
  loglik <- 0
  for (a in seq_len(K)) {
    for (b in a:K) {
      idx_a <- which(z == a)
      idx_b <- which(z == b)
      if (length(idx_a) == 0L || length(idx_b) == 0L) next
      
      if (a == b) {
        subA <- A[idx_a, idx_a, drop = FALSE]
        x <- subA[upper.tri(subA)]
      } else {
        subA <- A[idx_a, idx_b, drop = FALSE]
        x <- as.vector(subA)
      }
      if (length(x) == 0L) next
      
      p_ab <- mean(x)
      p_ab <- min(max(p_ab, 1e-9), 1 - 1e-9)  # clamp away from 0/1
      loglik <- loglik + sum(x * log(p_ab) + (1 - x) * log(1 - p_ab))
    }
  }
  
  N_obs <- n * (n - 1) / 2
  df    <- K * (K + 1) / 2 + (K - 1)    # symmetric B + block proportions
  BIC   <- -2 * loglik + df * log(N_obs)
  
  list(K = K, loglik = loglik, BIC = BIC, z = z)
}

## ---------- 3. Entropic KL-robust risk (dual representation) ----------

robust_kl_risk_entropic <- function(L_vec, w, C) {
  L_vec <- as.numeric(L_vec)
  w     <- as.numeric(w)
  
  if (any(!is.finite(L_vec))) stop("Non-finite loss values")
  if (any(w <= 0)) stop("Baseline weights must be strictly positive")
  sw <- sum(w)
  if (!isTRUE(all.equal(sw, 1))) w <- w / sw
  if (C < 0) stop("C must be >= 0")
  
  e0 <- sum(w * L_vec)
  if (C == 0 || max(L_vec) - min(L_vec) < .Machine$double.eps) return(e0)
  
  f_logeta <- function(log_eta) {
    eta <- exp(log_eta)
    log_mgf <- log(sum(w * exp(eta * L_vec)))
    (C + log_mgf) / eta
  }
  
  opt    <- optimize(f_logeta, c(log(1e-6), log(1e6)))
  f_star <- f_logeta(opt$minimum)
  risk   <- max(e0, min(max(L_vec), f_star))  # clamp to [e0, max L]
  risk
}

## ---------- 4. Per-scan SBM-vs-SBM robustness analysis ----------

analyze_scan_model_selection_sbm <- function(
    A,
    K_pair       = c(2, 3),
    tau_bic      = 0.25,
    C_grid       = seq(0, 0.25, length.out = 80),
    C_grid_small = seq(1e-5, 0.01, length.out = 40)
) {
  stopifnot(length(K_pair) == 2L)
  
  ## Fit the two SBMs
  fit1 <- fit_sbm_spectral(A, K_pair[1])
  fit2 <- fit_sbm_spectral(A, K_pair[2])
  
  bic <- c(fit1$BIC, fit2$BIC)
  names(bic) <- paste0("SBM_K", K_pair)
  
  ## Tempered BIC pseudo-posterior
  bic_shift <- bic - min(bic)                   # stabilise
  w_unnorm  <- exp(-0.5 * tau_bic * bic_shift)
  w         <- w_unnorm / sum(w_unnorm)
  
  ## 0–1 model-selection loss
  # states: true model in {1,2}; actions: choose model 1 or 2
  loss_mat <- rbind(
    c(0, 1),  # action = choose model 1
    c(1, 0)   # action = choose model 2
  )
  e0_actions <- as.numeric(loss_mat %*% w)
  action0_ix <- which.min(e0_actions)
  action0    <- if (action0_ix == 1L) names(bic)[1] else names(bic)[2]
  e0         <- e0_actions[action0_ix]
  
  ## Small-radius robust risk curve for baseline action
  e_rob <- vapply(
    C_grid_small,
    function(C) robust_kl_risk_entropic(loss_mat[action0_ix, ], w, C),
    numeric(1L)
  )
  
  norm_curve <- rep(NA_real_, length(C_grid_small))
  if (e0 > 0 && e0 < 1) {
    denom <- sqrt(2 * e0 * (1 - e0)) * sqrt(C_grid_small)
    norm_curve <- (e_rob - e0) / denom
    norm_curve[!is.finite(norm_curve)] <- NA_real_
  }
  
  ## Robust risks for both actions over a wider grid to get C*
  C_grid_big <- C_grid
  risk_a0 <- vapply(
    C_grid_big,
    function(C) robust_kl_risk_entropic(loss_mat[action0_ix, ], w, C),
    numeric(1L)
  )
  alt_ix <- 3L - action0_ix  # flip 1 <-> 2
  
  risk_alt <- vapply(
    C_grid_big,
    function(C) robust_kl_risk_entropic(loss_mat[alt_ix, ], w, C),
    numeric(1L)
  )
  
  diff_risk <- risk_a0 - risk_alt
  idx <- which(diff_risk > 0)
  if (length(idx) == 0L) {
    C_star <- Inf
  } else {
    j <- min(idx)
    if (j == 1L) {
      C_star <- C_grid_big[1L]
    } else {
      x1 <- C_grid_big[j - 1L]; x2 <- C_grid_big[j]
      y1 <- diff_risk[j - 1L];   y2 <- diff_risk[j]
      C_star <- if (y2 == y1) x2 else x1 + (0 - y1) * (x2 - y1) / (y2 - y1)
    }
  }
  
  list(
    K_pair        = K_pair,
    bic           = bic,
    w             = w,
    tau_bic       = tau_bic,
    action_baseline = action0,
    e0            = e0,
    C_grid_small  = C_grid_small,
    e_rob         = e_rob,
    norm_curve    = norm_curve,
    C_grid        = C_grid_big,
    risk_a0       = risk_a0,
    risk_alt      = risk_alt,
    C_star        = C_star
  )
}

## ======================================================
## 5. Run over all scans
## ======================================================
## Objects assumed to exist:
##   A_list      : list of adjacency matrices
##   subject_id  : vector length n_scans
##   scan_id     : vector length n_scans
##   group       : vector length n_scans

n_scans <- length(A_list)

K_pair        <- c(2, 3)               # the two SBMs we compare
tau_bic       <- 0.25                  # BIC temperature (tune if needed)
C_grid_model  <- seq(0, 0.25, length.out = 80)
C_grid_small  <- seq(1e-5, 0.01, length.out = 40)

model_res_list <- vector("list", n_scans)

fun_obs <- data.frame(
  scan    = seq_len(n_scans),
  subject = subject_id,
  scan_id = scan_id,
  C       = NA_real_,
  L       = NA_real_,
  S       = NA_real_,
  lambda1 = NA_real_
)

for (i in seq_len(n_scans)) {
  A_i <- A_list[[i]]
  
  ## SBM vs SBM robustness
  model_res_list[[i]] <- analyze_scan_model_selection_sbm(
    A            = A_i,
    K_pair       = K_pair,
    tau_bic      = tau_bic,
    C_grid       = C_grid_model,
    C_grid_small = C_grid_small
  )
  
  ## Observed functionals
  f_i <- compute_network_functionals(A_i)
  fun_obs$C[i]       <- f_i$C
  fun_obs$L[i]       <- f_i$L
  fun_obs$S[i]       <- f_i$S
  fun_obs$lambda1[i] <- f_i$lambda1
}

## ======================================================
## 6. Combine scan-level results
## ======================================================

model_summary <- data.frame(
  scan    = seq_len(n_scans),
  subject = subject_id,
  scan_id = scan_id,
  p_SBM1  = sapply(model_res_list, function(x) x$w[1]),
  p_SBM2  = sapply(model_res_list, function(x) x$w[2]),
  K1      = K_pair[1],
  K2      = K_pair[2],
  action0 = sapply(model_res_list, `[[`, "action_baseline"),
  e0      = sapply(model_res_list, `[[`, "e0"),
  C_star  = sapply(model_res_list, `[[`, "C_star")
)

brain_results <- merge(
  fun_obs, model_summary,
  by = c("scan", "subject", "scan_id")
)
brain_results$group <- group

## ======================================================
## 7. Subject-level summaries
## ======================================================

subject_summary <- aggregate(
  cbind(e0, p_SBM1, p_SBM2) ~ subject + group,
  data = brain_results,
  FUN  = mean
)

C_star_subject <- aggregate(
  C_star ~ subject,
  data = brain_results,
  FUN  = function(x) {
    x_finite <- x[is.finite(x)]
    if (length(x_finite) == 0L) Inf else min(x_finite)
  }
)

subject_summary <- merge(subject_summary, C_star_subject, by = "subject")
subject_summary <- subject_summary[order(subject_summary$subject), ]

## ======================================================
## 8. Figures (like the Mantziou brain plots, but SBM vs SBM)
## ======================================================

## 8(a) Small-radius normalised curve + S vs p_SBM2

e0_all <- sapply(model_res_list, `[[`, "e0")
scan_example <- which(e0_all > 0.05 & e0_all < 0.35)[1]
if (is.na(scan_example)) scan_example <- which(is.finite(e0_all))[1]

res_ex <- model_res_list[[scan_example]]

pdf("brain_exp2_SBMvsSBM_small_radius_and_S_vs_pSBM2.pdf",
    width = 7, height = 3.5)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1) + 0.1)

## (left) normalised robustness sensitivity curve
plot(
  sqrt(res_ex$C_grid_small),
  res_ex$norm_curve,
  type = "b", pch = 19,
  xlab = expression(sqrt(C)),
  ylab = expression(
    (e[rob] - e[0]) /
      (sqrt(2 * e[0] * (1 - e[0])) * sqrt(C))
  ),
  main = paste0(
    "Brain scan ", scan_example,
    ": SBM K=", K_pair[1], " vs K=", K_pair[2]
  )
)
abline(h = 1, lty = 2, col = "grey50")

## (right) S vs posterior mass on the more complex SBM (K2)
plot(
  brain_results$S,
  brain_results$p_SBM2,
  pch  = 19,
  xlab = "Small-world index S",
  ylab = "Posterior weight on SBM2",
  main = "Brain networks: S vs p(SBM2)"
)
dev.off()

## 8(b) C vs L, coloured by preferred SBM

cluster_label <- ifelse(brain_results$p_SBM2 >= 0.5, 2, 1)

pdf("brain_exp2_SBMvsSBM_C_vs_L_clusters.pdf", width = 5, height = 5)
plot(
  brain_results$C,
  brain_results$L,
  pch  = 19,
  col  = ifelse(cluster_label == 1, "darkgrey", "black"),
  xlab = "Clustering coefficient C",
  ylab = "Average shortest path length L",
  main = paste0(
    "Brain networks: small-world properties\n",
    "Preferred SBM: K=", K_pair[1], " (grey) vs K=", K_pair[2], " (black)"
  )
)
legend(
  "bottomleft",
  legend = c(
    paste0("SBM K=", K_pair[1]),
    paste0("SBM K=", K_pair[2])
  ),
  col    = c("darkgrey", "black"),
  pch    = 19,
  bty    = "n"
)
dev.off()



## ==========================================================
## Robust SBM(K1) vs SBM(K2) analysis with tempered evidence
## ==========================================================
## Assumes you already have:
##   brain_results with columns:
##   scan, subject, scan_id, group,
##   C, L, S, lambda1,
##   p_SBM1, p_SBM2, K1, K2
## ==========================================================

## -----------------------------
## 1. Temper model probabilities
## -----------------------------
tau_new <- 0.01      # strong tempering: odds_new = odds^tau_new
tau_old <- 1         # original probs assumed to use tau_old = 1
p_floor <- 1e-3      # hard floor away from 0 and 1

p1_orig <- brain_results$p_SBM1
p2_orig <- brain_results$p_SBM2

## Clean up 0/1 and renormalise
eps_prob <- .Machine$double.xmin
p1_orig <- pmax(pmin(p1_orig, 1 - eps_prob), eps_prob)
p2_orig <- pmax(pmin(p2_orig, 1 - eps_prob), eps_prob)
sums    <- p1_orig + p2_orig
p1_orig <- p1_orig / sums
p2_orig <- p2_orig / sums

## Retemper the odds: odds_new = odds^(tau_new / tau_old)
odds12_old <- p1_orig / p2_orig
odds12_new <- odds12_old^(tau_new / tau_old)

p1_temp <- odds12_new / (1 + odds12_new)
p2_temp <- 1 - p1_temp

## Optionally enforce a floor
p1_temp <- pmin(pmax(p1_temp, p_floor), 1 - p_floor)
p2_temp <- 1 - p1_temp

brain_results$p_SBM1_t <- p1_temp
brain_results$p_SBM2_t <- p2_temp

## Bayes action under tempered posterior
brain_results$action0 <- ifelse(
  brain_results$p_SBM2_t >= 0.5,
  paste0("SBM_K", brain_results$K2),
  paste0("SBM_K", brain_results$K1)
)

## Misclassification risk e0 = min(p1, p2)
brain_results$e0 <- pmin(brain_results$p_SBM1_t, brain_results$p_SBM2_t)

## ----------------------------------------------------------
## 2. Closed-form C* for binary 0–1 model selection
##    C*(e0) = KL( Bern(1/2) || Bern(e0) )
##        = 0.5*log(0.25 / (e0*(1-e0)))
## ----------------------------------------------------------
e0_clamped <- pmin(pmax(brain_results$e0, p_floor), 0.5 - 1e-6)
brain_results$C_star <- 0.5 * log(0.25 / (e0_clamped * (1 - e0_clamped)))

## ----------------------------------------------------------
## 3. Subject-level summaries: mean probs, mean e0, min C*
## ----------------------------------------------------------
subject_summary <- aggregate(
  cbind(e0       = brain_results$e0,
        p_SBM1_t = brain_results$p_SBM1_t,
        p_SBM2_t = brain_results$p_SBM2_t) ~ subject + group,
  data = brain_results,
  FUN  = mean
)

C_star_subject <- aggregate(
  C_star ~ subject,
  data = brain_results,
  FUN  = min   # most vulnerable scan for that subject
)

subject_summary <- merge(subject_summary, C_star_subject, by = "subject")
subject_summary <- subject_summary[order(subject_summary$subject), ]

## Peek at the first few lines
print(head(subject_summary), digits = 3)

## ----------------------------------------------------------
## 4. Watson–Holmes robustness curve for a 2–point posterior
## ----------------------------------------------------------
robust_curve_binary <- function(e0, C_max = 4, n_grid = 50L) {
  if (!is.finite(e0) || e0 <= 0 || e0 >= 0.5) {
    return(list(
      e0         = e0,
      C_grid     = NA_real_,
      e_rob      = NA_real_,
      norm_curve = NA_real_
    ))
  }
  
  w_bad  <- e0
  w_good <- 1 - e0
  
  ## Tilt parameter grid
  lambda_grid <- seq(0, 10, length.out = 2001L)
  
  ## Under tilt, “bad” probability q(lambda)
  q_lambda <- w_bad * exp(lambda_grid) /
    (w_good + w_bad * exp(lambda_grid))
  
  ## KL( Bern(q) || Bern(w_bad) )
  C_lambda <- q_lambda * log(q_lambda / w_bad) +
    (1 - q_lambda) * log((1 - q_lambda) / (1 - w_bad))
  
  C_lambda_finite <- C_lambda[is.finite(C_lambda)]
  q_lambda_finite <- q_lambda[is.finite(C_lambda)]
  
  C_max_eff <- min(C_max, max(C_lambda_finite) * 0.999)
  C_min_eff <- max(min(C_lambda_finite), 1e-6)
  
  C_grid <- seq(C_min_eff, C_max_eff, length.out = n_grid)
  
  ## Invert C(lambda) -> q(C) by interpolation
  q_grid <- approx(
    x    = C_lambda_finite,
    y    = q_lambda_finite,
    xout = C_grid,
    rule = 2
  )$y
  
  e_rob <- q_grid
  
  norm_curve <- (e_rob - e0) /
    (sqrt(2 * e0 * (1 - e0)) * sqrt(C_grid))
  
  list(
    e0         = e0,
    C_grid     = C_grid,
    e_rob      = e_rob,
    norm_curve = norm_curve
  )
}

## Compute robustness object per scan
model_res_list <- lapply(seq_len(nrow(brain_results)), function(i) {
  rc <- robust_curve_binary(brain_results$e0[i])
  rc$C_star <- brain_results$C_star[i]
  rc
})
names(model_res_list) <- brain_results$scan

e0_all <- sapply(model_res_list, `[[`, "e0")

## Choose a scan with moderate e0 (say between 0.1 and 0.4)
scan_example <- which(e0_all > 0.1 & e0_all < 0.4)[1]
if (is.na(scan_example)) scan_example <- which(is.finite(e0_all))[1]

res_ex <- model_res_list[[scan_example]]

## ----------------------------------------------------------
## 5. Plot: robustness curve + S vs P(SBM_K2 | data)
## ----------------------------------------------------------
K_pair <- sort(unique(c(brain_results$K1, brain_results$K2)))

pdf("brain_exp2_SBMvsSBM_small_radius_and_S_vs_pSBM2.pdf",
    width = 7, height = 3.5)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1) + 0.1)

## (left) normalised robustness sensitivity curve
idx_ok <- which(
  is.finite(res_ex$C_grid) &
    is.finite(res_ex$norm_curve)
)

plot(
  sqrt(res_ex$C_grid[idx_ok]),
  res_ex$norm_curve[idx_ok],
  type = "b", pch = 19,
  xlab = expression(sqrt(C)),
  ylab = expression(
    (e[rob] - e[0]) /
      (sqrt(2 * e[0] * (1 - e[0])) * sqrt(C))
  ),
  main = paste0(
    "COBRE scan ", scan_example,
    ": SBM K=", brain_results$K1[scan_example],
    " vs K=", brain_results$K2[scan_example]
  )
)
abline(h = 1, lty = 2)

## (right) S vs tempered P(SBM_K2 | data)
plot(
  brain_results$S,
  brain_results$p_SBM2_t,
  pch  = 19,
  xlab = "Small-world index S",
  ylab = "Tempered posterior P(SBM K2 | data)",
  main = "Brain networks: S vs P(SBM K2 | data)"
)

dev.off()

## ----------------------------------------------------------
## 6. C vs L, coloured by preferred SBM (tempered)
## ----------------------------------------------------------
cluster_label <- ifelse(brain_results$p_SBM2_t >= 0.5, 2, 1)

pdf("brain_exp2_SBMvsSBM_C_vs_L_clusters.pdf", width = 5, height = 5)
plot(
  brain_results$C,
  brain_results$L,
  pch  = 19,
  col  = ifelse(cluster_label == 1, "darkgrey", "black"),
  xlab = "Clustering coefficient C",
  ylab = "Average shortest path length L",
  main = paste0(
    "Brain networks: small-world properties\n",
    "Preferred SBM: K=", K_pair[1], " (grey) vs K=",
    K_pair[2], " (black)"
  )
)
legend(
  "bottomleft",
  legend = c(
    paste0("SBM K=", K_pair[1]),
    paste0("SBM K=", K_pair[2])
  ),
  col = c("darkgrey", "black"),
  pch = 19,
  bty = "n"
)
dev.off()




## ---------- 2b. Simple RDPG fit via adjacency spectral embedding ----------
fit_rdpg_ase <- function(A, d, eps_prob = 1e-9) {
  A <- as.matrix(A)
  n <- nrow(A)
  if (n != ncol(A)) stop("A must be square")
  
  ## symmetrise and remove self-loops
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) <- 0
  
  ev  <- eigen(A, symmetric = TRUE)
  ord <- order(abs(ev$values), decreasing = TRUE)
  vals <- ev$values[ord]
  vecs <- ev$vectors[, ord, drop = FALSE]
  
  d_eff  <- min(d, length(vals))
  vals_d <- vals[seq_len(d_eff)]
  vecs_d <- vecs[, seq_len(d_eff), drop = FALSE]
  
  ## rank-d spectral reconstruction of the edge-probability matrix
  P_hat <- vecs_d %*% diag(vals_d, nrow = d_eff) %*% t(vecs_d)
  P_hat[lower.tri(P_hat)] <- t(P_hat)[lower.tri(P_hat)]
  diag(P_hat) <- 0
  
  ## clamp to (0,1) to avoid log(0)
  p <- P_hat[upper.tri(P_hat)]
  p <- pmax(pmin(p, 1 - eps_prob), eps_prob)
  
  x <- A[upper.tri(A)]
  
  loglik <- sum(x * log(p) + (1 - x) * log(1 - p))
  
  N_obs <- n * (n - 1) / 2
  df    <- n * d_eff                 # number of latent coordinates
  BIC   <- -2 * loglik + df * log(N_obs)
  
  list(
    model  = "RDPG",
    d      = d_eff,
    loglik = loglik,
    BIC    = BIC
  )
}



## ----------------------------------------------------------
## Watson–Holmes robustness curve for a 2–point posterior
## ----------------------------------------------------------
robust_curve_binary <- function(e0, C_max = 4, n_grid = 50L) {
  if (!is.finite(e0) || e0 <= 0 || e0 >= 0.5) {
    return(list(
      e0         = e0,
      C_grid     = NA_real_,
      e_rob      = NA_real_,
      norm_curve = NA_real_
    ))
  }
  
  w_bad  <- e0
  w_good <- 1 - e0
  
  ## Tilt parameter grid
  lambda_grid <- seq(0, 10, length.out = 2001L)
  
  ## Under tilt, “bad” probability q(lambda)
  q_lambda <- w_bad * exp(lambda_grid) /
    (w_good + w_bad * exp(lambda_grid))
  
  ## KL( Bern(q) || Bern(w_bad) )
  C_lambda <- q_lambda * log(q_lambda / w_bad) +
    (1 - q_lambda) * log((1 - q_lambda) / (1 - w_bad))
  
  C_lambda_finite <- C_lambda[is.finite(C_lambda)]
  q_lambda_finite <- q_lambda[is.finite(C_lambda)]
  
  if (length(C_lambda_finite) == 0L) {
    return(list(
      e0         = e0,
      C_grid     = NA_real_,
      e_rob      = NA_real_,
      norm_curve = NA_real_
    ))
  }
  
  C_max_eff <- min(C_max, max(C_lambda_finite) * 0.999)
  C_min_eff <- max(min(C_lambda_finite), 1e-6)
  
  C_grid <- seq(C_min_eff, C_max_eff, length.out = n_grid)
  
  ## Invert C(lambda) -> q(C) by interpolation
  q_grid <- approx(
    x    = C_lambda_finite,
    y    = q_lambda_finite,
    xout = C_grid,
    rule = 2
  )$y
  
  e_rob <- q_grid
  
  norm_curve <- (e_rob - e0) /
    (sqrt(2 * e0 * (1 - e0)) * sqrt(C_grid))
  
  list(
    e0         = e0,
    C_grid     = C_grid,
    e_rob      = e_rob,
    norm_curve = norm_curve
  )
}


## ---------- 4. Per-scan SBM(K) vs RDPG(d) robustness analysis ----------
analyze_scan_sbm_vs_rdpg <- function(
    A,
    K_sbm   = 3L,
    d_rdpg  = 3L,
    tau_bic = 0.25,
    C_max   = 4,
    n_grid  = 80L
) {
  ## Fit SBM(K) and RDPG(d)
  fit_sbm  <- fit_sbm_spectral(A, K = K_sbm)
  fit_rdpg <- fit_rdpg_ase(A, d = d_rdpg)
  
  ## Tempered BIC pseudo-posterior:
  ##   Pi(M = m | G) ∝ exp{ - (tau_bic/2) * BIC_m(G) }
  bic <- c(SBM = fit_sbm$BIC, RDPG = fit_rdpg$BIC)
  bic_shift <- bic - min(bic)   # stabilise
  w_unnorm  <- exp(-0.5 * tau_bic * bic_shift)
  w         <- w_unnorm / sum(w_unnorm)
  p_sbm  <- as.numeric(w["SBM"])
  p_rdpg <- as.numeric(w["RDPG"])
  
  ## Bayes action under 0–1 loss: pick more probable model
  action0 <- if (p_rdpg >= 0.5) "RDPG" else "SBM"
  
  ## Misclassification risk e0 = min(p_sbm, p_rdpg)
  e0 <- min(p_sbm, p_rdpg)
  
  ## Robust curve and C*(e0) = KL(Bern(1/2) || Bern(e0))
  rc <- robust_curve_binary(e0, C_max = C_max, n_grid = n_grid)
  
  e0_clamped <- pmin(pmax(e0, 1e-6), 0.5 - 1e-6)
  C_star     <- 0.5 * log(0.25 / (e0_clamped * (1 - e0_clamped)))
  
  list(
    K_sbm      = K_sbm,
    d_rdpg     = d_rdpg,
    BIC_sbm    = fit_sbm$BIC,
    BIC_rdpg   = fit_rdpg$BIC,
    p_sbm      = p_sbm,
    p_rdpg     = p_rdpg,
    action0    = action0,
    e0         = e0,
    C_star     = C_star,
    C_grid     = rc$C_grid,
    e_rob      = rc$e_rob,
    norm_curve = rc$norm_curve
  )
}



## ======================================================
## SBM(K_sbm) vs RDPG(d_rdpg) over all brain scans
## ======================================================

n_scans <- length(A_list)

K_sbm   <- 3L       # e.g. K = 3 SBM
d_rdpg  <- 3L       # e.g. d = 3 RDPG
tau_bic <- 0.25     # same temperature as in SBM-vs-SBM

model_res_sbm_rdpg <- vector("list", n_scans)

for (i in seq_len(n_scans)) {
  A_i <- A_list[[i]]
  
  ## Per-scan robust model comparison
  model_res_sbm_rdpg[[i]] <- analyze_scan_sbm_vs_rdpg(
    A      = A_i,
    K_sbm  = K_sbm,
    d_rdpg = d_rdpg,
    tau_bic = tau_bic,
    C_max   = 4,
    n_grid  = 80L
  )
}

## ------------------------------------------------------
## Combine into a scan-level data frame
## ------------------------------------------------------

model_summary_sbm_rdpg <- data.frame(
  scan    = seq_len(n_scans),
  subject = subject_id,
  scan_id = scan_id,
  p_SBM   = sapply(model_res_sbm_rdpg, `[[`, "p_sbm"),
  p_RDPG  = sapply(model_res_sbm_rdpg, `[[`, "p_rdpg"),
  K_sbm   = K_sbm,
  d_rdpg  = d_rdpg,
  action0 = sapply(model_res_sbm_rdpg, `[[`, "action0"),
  e0      = sapply(model_res_sbm_rdpg, `[[`, "e0"),
  C_star  = sapply(model_res_sbm_rdpg, `[[`, "C_star")
)

brain_results_sbm_rdpg <- merge(
  fun_obs,                      # C, L, S, lambda1 per scan
  model_summary_sbm_rdpg,
  by = c("scan", "subject", "scan_id")
)

brain_results_sbm_rdpg$group <- group




## ------------------------------------------------------
## Subject-level summaries: mean probs, mean e0, min C*
## ------------------------------------------------------

subject_summary_sbm_rdpg <- aggregate(
  cbind(e0    = brain_results_sbm_rdpg$e0,
        p_SBM = brain_results_sbm_rdpg$p_SBM,
        p_RDPG = brain_results_sbm_rdpg$p_RDPG) ~ subject + group,
  data = brain_results_sbm_rdpg,
  FUN  = mean
)

C_star_subject_sbm_rdpg <- aggregate(
  C_star ~ subject,
  data = brain_results_sbm_rdpg,
  FUN  = function(x) {
    x_finite <- x[is.finite(x)]
    if (length(x_finite) == 0L) Inf else min(x_finite)
  }
)

subject_summary_sbm_rdpg <- merge(
  subject_summary_sbm_rdpg,
  C_star_subject_sbm_rdpg,
  by = "subject"
)

subject_summary_sbm_rdpg <- subject_summary_sbm_rdpg[
  order(subject_summary_sbm_rdpg$subject), ]

## Global summaries of switching radius
C_star_finite <- brain_results_sbm_rdpg$C_star[
  is.finite(brain_results_sbm_rdpg$C_star)
]

summary_C_star <- summary(C_star_finite)
frac_flip_001  <- mean(C_star_finite <= 0.01)  # flips before C = 0.01
frac_flip_010  <- mean(C_star_finite <= 0.10)  # flips before C = 0.10

## You can inspect:
print(summary_C_star)
cat("Fraction scans switching by C=0.01:", frac_flip_001, "\n")
cat("Fraction scans switching by C=0.10:", frac_flip_010, "\n")



## Choose a scan with moderate e0 (say between 0.1 and 0.4)
e0_all_sbm_rdpg <- sapply(model_res_sbm_rdpg, `[[`, "e0")
scan_example_sbm_rdpg <- which(e0_all_sbm_rdpg > 0.1 & e0_all_sbm_rdpg < 0.4)[1]
if (is.na(scan_example_sbm_rdpg)) {
  scan_example_sbm_rdpg <- which(is.finite(e0_all_sbm_rdpg))[1]
}

res_ex <- model_res_sbm_rdpg[[scan_example_sbm_rdpg]]

pdf("brain_exp2_SBM_vs_RDPG_small_radius_and_S_vs_pRDPG.pdf",
    width = 7, height = 3.5)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1) + 0.1)

## (left) normalised robustness sensitivity curve
idx_ok <- which(
  is.finite(res_ex$C_grid) &
    is.finite(res_ex$norm_curve)
)

plot(
  sqrt(res_ex$C_grid[idx_ok]),
  res_ex$norm_curve[idx_ok],
  type = "b", pch = 19,
  xlab = expression(sqrt(C)),
  ylab = expression(
    (e[rob] - e[0]) /
      (sqrt(2 * e[0] * (1 - e[0])) * sqrt(C))
  ),
  main = paste0(
    "Brain scan ", scan_example_sbm_rdpg,
    ": SBM K=", K_sbm, " vs RDPG d=", d_rdpg
  )
)
abline(h = 1, lty = 2, col = "grey50")

## (right) S vs posterior mass on the RDPG model
plot(
  brain_results_sbm_rdpg$S,
  brain_results_sbm_rdpg$p_RDPG,
  pch  = 19,
  xlab = "Small-world index S",
  ylab = "Posterior weight on RDPG",
  main = "Brain networks: S vs P(RDPG | data)"
)

dev.off()


cluster_label_sbm_rdpg <- ifelse(
  brain_results_sbm_rdpg$p_RDPG >= 0.5, "RDPG", "SBM"
)

pdf("brain_exp2_SBM_vs_RDPG_C_vs_L_clusters.pdf", width = 5, height = 5)
plot(
  brain_results_sbm_rdpg$C,
  brain_results_sbm_rdpg$L,
  pch  = 19,
  col  = ifelse(cluster_label_sbm_rdpg == "SBM", "darkgrey", "black"),
  xlab = "Clustering coefficient C",
  ylab = "Average shortest path length L",
  main = paste0(
    "Brain networks: small-world properties\n",
    "Preferred model: SBM K=", K_sbm, " (grey) vs ",
    "RDPG d=", d_rdpg, " (black)"
  )
)
legend(
  "bottomleft",
  legend = c(
    paste0("SBM K=", K_sbm),
    paste0("RDPG d=", d_rdpg)
  ),
  col = c("darkgrey", "black"),
  pch = 19,
  bty = "n"
)
dev.off()






