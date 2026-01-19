## ==========================================================
## Experiment 2 on COBRE data: robust ER vs SBM + small-world
## ==========================================================

## 0. Packages and data -------------------------------------

library(graphclass)
library(igraph)

## COBRE brain networks (Relión et al., via graphclass)
data("COBRE.data")

n_scans    <- nrow(COBRE.data$X.cobre)
subject_id <- COBRE.data$subject.label
scan_id    <- rep(1L, n_scans)
group      <- COBRE.data$Y.cobre   # 0 = control, 1 = schizophrenia

## 1. Build binary adjacencies from Fisher-z correlations ----

threshold_to_binary <- function(A, tau = 0.3) {
  A <- as.matrix(A)
  A <- (abs(A) > tau)
  diag(A) <- 0
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  storage.mode(A) <- "integer"
  A
}

A_list <- vector("list", n_scans)
for (i in seq_len(n_scans)) {
  mat_i <- graphclass::get_matrix(COBRE.data$X.cobre[i, ])
  A_list[[i]] <- threshold_to_binary(mat_i, tau = 0.3)
}

p_nodes <- nrow(A_list[[1]])

## 2. KL-ball robust risk via entropic tilting ----------------

sanitize_weights <- function(w, eps = 1e-6) {
  w <- as.numeric(w)
  if (any(w < 0)) stop("Weights must be non-negative.")
  sw <- sum(w)
  if (sw <= 0) stop("Sum of weights must be positive.")
  w <- w / sw
  w <- pmax(w, eps)
  w / sum(w)
}

# robust risk  sup_{q:KL(q||w)<=C} sum_s q_s * L_s
robust_kl_risk_entropic <- function(losses, w, C,
                                    tol_L = 1e-12,
                                    eps_w = 1e-6,
                                    maxit_bracket = 60) {
  L <- as.numeric(losses)
  w <- sanitize_weights(w, eps = eps_w)
  L_range  <- max(L) - min(L)
  base_risk <- sum(w * L)
  
  # trivial cases
  if (C <= 0 || L_range < tol_L) {
    return(base_risk)
  }
  
  # helper: for a given lambda, compute KL(lambda)-C and rho(lambda)
  g_and_rho <- function(lambda) {
    if (lambda <= 0) {
      rho <- base_risk
      kl  <- 0
    } else {
      z    <- lambda * L
      zmax <- max(z)
      expz <- exp(z - zmax)
      Z    <- sum(w * expz)
      q    <- w * expz / Z
      rho  <- sum(q * L)
      logZ <- log(Z) + zmax
      kl   <- lambda * rho - logZ
    }
    list(g = kl - C, kl = kl, rho = rho)
  }
  
  # bracket a root in lambda >= 0
  lambda_lo <- 0
  out_lo    <- g_and_rho(lambda_lo)
  kl_lo     <- out_lo$kl
  
  lambda_hi <- 1
  out_hi    <- g_and_rho(lambda_hi)
  kl_hi     <- out_hi$kl
  it        <- 0
  
  while (is.finite(kl_hi) && kl_hi < C && it < maxit_bracket) {
    lambda_hi <- lambda_hi * 2
    out_hi    <- g_and_rho(lambda_hi)
    kl_hi     <- out_hi$kl
    it        <- it + 1
  }
  
  # if we never reach KL >= C in a stable way, just return the last rho
  if (!is.finite(kl_hi) || kl_hi <= 0) {
    return(out_hi$rho)
  }
  
  g_only <- function(lambda) g_and_rho(lambda)$g
  
  root <- tryCatch(
    uniroot(g_only, interval = c(lambda_lo, lambda_hi), tol = 1e-10),
    error = function(e) NULL
  )
  
  if (is.null(root)) {
    return(out_hi$rho)
  }
  
  lambda_star <- root$root
  out_star    <- g_and_rho(lambda_star)
  out_star$rho
}

## 3. Network functionals: C, L, S, lambda1 ------------------

compute_network_functionals <- function(A) {
  g <- graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
  
  # global clustering coefficient
  C <- suppressWarnings(transitivity(g, type = "global"))
  
  # average shortest path on giant component
  comps    <- components(g)
  giant_id <- which.max(comps$csize)
  v_giant  <- which(comps$membership == giant_id)
  g_giant  <- induced_subgraph(g, v_giant)
  L <- suppressWarnings(average.path.length(g_giant, directed = FALSE))
  
  # ER baseline with same density
  n      <- vcount(g)
  m      <- gsize(g)
  n_pair <- n * (n - 1) / 2
  p_hat  <- if (n_pair > 0) m / n_pair else NA_real_
  
  C_ER <- p_hat
  if (!is.na(p_hat) && p_hat > 0 && n * p_hat > 1) {
    L_ER <- log(n) / log(n * p_hat)
  } else {
    L_ER <- NA_real_
  }
  
  S <- if (!any(is.na(c(C, L, C_ER, L_ER)))) {
    (C / C_ER) / (L / L_ER)
  } else {
    NA_real_
  }
  
  eig     <- eigen(A, symmetric = TRUE, only.values = TRUE)
  lambda1 <- Re(eig$values[1])
  
  list(
    C       = C,
    L       = L,
    S       = S,
    lambda1 = lambda1,
    p_hat   = p_hat,
    C_ER    = C_ER,
    L_ER    = L_ER
  )
}

## 4. Working models: ER(p) and a simple K-block SBM ----------

loglik_ER <- function(A, p) {
  upper <- upper.tri(A)
  x     <- A[upper]
  m1    <- sum(x)
  m     <- length(x)
  if (p <= 0 || p >= 1) return(-Inf)
  m1 * log(p) + (m - m1) * log(1 - p)
}

fit_ER_model <- function(A, alpha0 = 1, beta0 = 1) {
  upper <- upper.tri(A)
  x     <- A[upper]
  e     <- sum(x)
  m     <- length(x)
  
  p_hat      <- e / m
  loglik_hat <- loglik_ER(A, p_hat)
  k_par      <- 1
  BIC        <- -2 * loglik_hat + k_par * log(m)
  
  alpha_post <- alpha0 + e
  beta_post  <- beta0 + m - e
  
  list(
    p_hat      = p_hat,
    alpha_post = alpha_post,
    beta_post  = beta_post,
    BIC        = BIC,
    loglik_hat = loglik_hat,
    m          = m
  )
}

fit_SBM_model <- function(A, K = 2) {
  n <- nrow(A)
  if (nrow(A) != ncol(A)) stop("A must be square.")
  
  eig <- eigen(A, symmetric = TRUE)
  U   <- eig$vectors[, 1:K, drop = FALSE]
  cl  <- kmeans(U, centers = K, nstart = 10)
  z   <- cl$cluster
  
  B        <- matrix(0, K, K)
  n_pairs  <- matrix(0, K, K)
  e_pairs  <- matrix(0, K, K)
  
  for (i in seq_len(K)) {
    idx_i <- which(z == i)
    for (j in i:K) {
      idx_j <- which(z == j)
      if (i == j) {
        if (length(idx_i) >= 2) {
          pairs <- combn(idx_i, 2)
          edges <- A[cbind(pairs[1, ], pairs[2, ])]
        } else {
          edges <- integer(0)
        }
      } else {
        if (length(idx_i) > 0 && length(idx_j) > 0) {
          edges <- A[idx_i, idx_j, drop = FALSE]
        } else {
          edges <- integer(0)
        }
      }
      m_ij <- length(edges)
      e_ij <- sum(edges)
      n_pairs[i, j] <- n_pairs[j, i] <- m_ij
      e_pairs[i, j] <- e_pairs[j, i] <- e_ij
      p_ij <- if (m_ij > 0) e_ij / m_ij else 0
      B[i, j] <- B[j, i] <- p_ij
    }
  }
  
  loglik <- 0
  for (i in seq_len(K)) {
    for (j in i:K) {
      m_ij <- n_pairs[i, j]
      if (m_ij == 0) next
      e_ij <- e_pairs[i, j]
      p_ij <- B[i, j]
      p_ij <- pmin(pmax(p_ij, 1e-6), 1 - 1e-6)
      loglik <- loglik +
        e_ij * log(p_ij) + (m_ij - e_ij) * log(1 - p_ij)
    }
  }
  
  m_total <- sum(n_pairs[upper.tri(n_pairs, diag = TRUE)])
  k_par   <- K * (K + 1) / 2
  BIC     <- -2 * loglik + k_par * log(m_total)
  
  list(
    K       = K,
    z       = z,
    B       = B,
    loglik  = loglik,
    BIC     = BIC,
    m_total = m_total
  )
}

## 5. Robust ER vs SBM model selection for one scan --------------

analyze_scan_model_selection <- function(A,
                                         K_range = 2:4,
                                         C_grid  = 10^seq(-4, -1,
                                                          length.out = 12)) {
  er_fit <- fit_ER_model(A)
  
  sbm_list <- lapply(K_range, function(K) fit_SBM_model(A, K = K))
  BIC_sbm  <- sapply(sbm_list, `[[`, "BIC")
  idx_best <- which.min(BIC_sbm)
  sbm_best <- sbm_list[[idx_best]]
  K_best   <- K_range[idx_best]
  
  # BIC ~ -2 log marginal likelihood
  log_marg_ER    <- -0.5 * er_fit$BIC
  log_marg_sbm_K <- -0.5 * BIC_sbm
  log_marg <- c(
    ER = log_marg_ER,
    setNames(log_marg_sbm_K, paste0("SBM", K_range))
  )
  
  log_shift <- log_marg - max(log_marg)
  w_models  <- exp(log_shift)
  w_models  <- w_models / sum(w_models)
  
  p_ER  <- w_models["ER"]
  p_SBM <- 1 - p_ER
  
  # two-point posterior (ER vs "some SBM")
  w2 <- sanitize_weights(c(p_ER, p_SBM), eps = 1e-6)
  names(w2) <- c("ER", "SBM")
  
  # baseline Bayes action under 0-1 loss
  action_idx  <- which.max(w2)  # 1 = ER, 2 = SBM
  action_name <- names(w2)[action_idx]
  
  L_vec <- if (action_idx == 1L) c(0, 1) else c(1, 0)
  
  e0 <- sum(w2 * L_vec)
  
  # robust misclassification for this action
  e_rob <- sapply(C_grid, function(C) {
    robust_kl_risk_entropic(L_vec, w2, C)
  })
  
  # normalised robustness sensitivity curve
  if (e0 <= 0 || e0 >= 1 || !is.finite(e0)) {
    norm_curve <- rep(NA_real_, length(C_grid))
  } else {
    norm_curve <- (e_rob - e0) /
      (sqrt(2 * e0 * (1 - e0)) * sqrt(C_grid))
  }
  
  # robust risk for each possible action
  L_ER  <- c(0, 1)
  L_SBM <- c(1, 0)
  risk_ER <- sapply(C_grid, function(C) {
    robust_kl_risk_entropic(L_ER,  w2, C)
  })
  risk_SBM <- sapply(C_grid, function(C) {
    robust_kl_risk_entropic(L_SBM, w2, C)
  })
  action_rob_idx <- ifelse(risk_ER <= risk_SBM, 1L, 2L)
  action_rob     <- c("ER", "SBM")[action_rob_idx]
  
  idx_change <- which(action_rob_idx != action_idx)
  C_star <- if (length(idx_change) == 0L) Inf else C_grid[min(idx_change)]
  
  list(
    w_models        = w_models,
    p_ER            = w2["ER"],
    p_SBM           = w2["SBM"],
    K_best          = K_best,
    action_baseline = action_name,
    e0              = e0,
    e_rob           = e_rob,
    norm_curve      = norm_curve,
    C_grid          = C_grid,
    risk_ER         = risk_ER,
    risk_SBM        = risk_SBM,
    action_rob      = action_rob,
    C_star          = C_star
  )
}

## 6. Run Experiment 2 on all COBRE scans -----------------------

C_grid_model <- 10^seq(-4, -1, length.out = 12)
K_range      <- 2:4

n_scans <- length(A_list)

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
  if (i %% 20 == 0) message("Scan ", i, " / ", n_scans)
  
  A_i <- A_list[[i]]
  
  # ER vs SBM robustness
  model_res_list[[i]] <- analyze_scan_model_selection(
    A       = A_i,
    K_range = K_range,
    C_grid  = C_grid_model
  )
  
  # observed functionals
  f_i <- compute_network_functionals(A_i)
  fun_obs$C[i]       <- f_i$C
  fun_obs$L[i]       <- f_i$L
  fun_obs$S[i]       <- f_i$S
  fun_obs$lambda1[i] <- f_i$lambda1
}

## 7. Combine results and subject-level summaries ---------------

model_summary <- data.frame(
  scan      = seq_len(n_scans),
  subject   = subject_id,
  scan_id   = scan_id,
  p_ER      = sapply(model_res_list, `[[`, "p_ER"),
  p_SBM     = sapply(model_res_list, `[[`, "p_SBM"),
  K_best    = sapply(model_res_list, `[[`, "K_best"),
  action0   = sapply(model_res_list, `[[`, "action_baseline"),
  e0        = sapply(model_res_list, `[[`, "e0"),
  C_star    = sapply(model_res_list, `[[`, "C_star")
)

brain_results <- merge(
  fun_obs, model_summary,
  by = c("scan", "subject", "scan_id")
)
brain_results$group <- group

## Subject-level averages
subject_summary <- aggregate(
  cbind(e0, p_ER, p_SBM) ~ subject + group,
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

## 8. Illustrative figures -------------------------------------

## 8(a) Robustness curve for one scan + S vs p(SBM | G)

e0_all <- sapply(model_res_list, `[[`, "e0")
scan_example <- which(is.finite(e0_all))[1]
res_ex <- model_res_list[[scan_example]]

pdf("brain_exp2_small_radius_and_S_vs_pSBM.pdf",
    width = 7, height = 3.5)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1) + 0.1)

y <- res_ex$norm_curve
good <- is.finite(y)

if (any(good)) {
  plot(
    sqrt(res_ex$C_grid[good]),
    y[good],
    type = "b", pch = 19,
    xlab = expression(sqrt(C)),
    ylab = expression(
      (e[rob] - e[0]) /
        (sqrt(2 * e[0] * (1 - e[0])) * sqrt(C))
    ),
    main = paste0("COBRE scan ", scan_example, ": ER vs SBM")
  )
  abline(h = 1, lty = 2)
} else {
  ## fallback (should not happen with the code above)
  plot(
    sqrt(res_ex$C_grid),
    res_ex$e_rob,
    type = "b", pch = 19,
    xlab = expression(sqrt(C)),
    ylab = expression(e[rob](C)),
    main = paste0("COBRE scan ", scan_example, ": ER vs SBM")
  )
}

plot(
  brain_results$S,
  brain_results$p_SBM,
  pch  = 19,
  xlab = "Observed small-world index S",
  ylab = "Posterior probability of SBM",
  main = "Small-world structure vs ER/SBM choice"
)

dev.off()

## 8(b) C vs L coloured by ER/SBM preference (cf. Mantziou Fig. 25)

cluster_label <- ifelse(brain_results$p_SBM >= 0.5, 2, 1)

pdf("brain_exp2_C_vs_L_clusters.pdf", width = 5, height = 5)
plot(
  brain_results$C,
  brain_results$L,
  pch  = 19,
  col  = ifelse(cluster_label == 1, "darkgrey", "black"),
  xlab = "Clustering coefficient C",
  ylab = "Average shortest path length L",
  main = "Brain networks: small-world properties and ER/SBM choice"
)
legend(
  "bottomleft",
  legend = c("ER-like", "SBM-like"),
  col    = c("darkgrey", "black"),
  pch    = 19,
  bty    = "n"
)
dev.off()
