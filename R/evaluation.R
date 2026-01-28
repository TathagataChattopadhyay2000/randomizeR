## ------------------------------------------------------------
## Metric Calculation Functions
## ------------------------------------------------------------

# Average Standard Deviation (ASD)
calculate_asd <- function(assignments, Nsim, N, K, ratio, groups) {
  if (Nsim == 1L) return(rep(NA_real_, N))

  asd_values      <- numeric(N)
  all_proportions <- array(0, dim = c(Nsim, N, K))

  for (sim in seq_len(Nsim)) {
    sim_data <- assignments[assignments$Simulation == sim, ]
    counts   <- integer(K)

    for (m in seq_len(N)) {
      idx <- which(groups == as.character(sim_data$treatment[m]))
      counts[idx] <- counts[idx] + 1L
      all_proportions[sim, m, ] <- counts / m
    }
  }

  for (m in seq_len(N)) {
    group_props <- matrix(all_proportions[, m, ], nrow = Nsim, ncol = K)
    sd_values   <- apply(group_props, 2L, sd)
    asd_values[m] <- sqrt(m * sum(sd_values^2))
  }

  asd_values
}

# Imbalance at final step (not used directly inside evaluate_assignments,
# but kept as a utility)
calculate_imbalance <- function(assignments, ratio) {
  n      <- nrow(assignments)
  Ni     <- table(assignments$treatment)
  target <- n * ratio / sum(ratio)
  sqrt(sum((Ni - target)^2))
}

# FI_n for a single simulation given conditional prob matrix P and rho
# (not directly used now, but kept as separate utility)
calculate_FI_n <- function(P, target_proportions) {
  if (nrow(P) != length(target_proportions)) {
    stop("P row count mismatch with target_proportions.")
  }
  FIs <- apply(P, 2L, function(P_j) sqrt(sum((P_j - target_proportions)^2)))
  mean(FIs)
}

# Empirical allocation proportions π_i(m) = N_i(m)/m
calculate_pi_matrix <- function(assignments, N, K, groups) {
  pi_matrix <- matrix(0, nrow = K, ncol = N)
  for (j in seq_len(N)) {
    for (i in seq_len(K)) {
      count_ij <- sum(assignments$treatment[1:j] == groups[i])
      pi_matrix[i, j] <- count_ij / j
    }
  }
  pi_matrix
}

# Unified performance measures (UI, UR, G)
calculate_performance_measures <- function(
    MPM_n_xi, MPM_n_MaxEnt, MPM_n_CRD,
    AFI_n_xi, AFI_n_MaxEnt, AFI_n_CRD,
    wI = 1, wR = 1
) {
  UI <- (MPM_n_xi - MPM_n_MaxEnt) / (MPM_n_CRD - MPM_n_MaxEnt)
  UR <- (AFI_n_xi) / (AFI_n_MaxEnt)

  G <- sqrt((wI * UI)^2 + (wR * UR)^2) / sqrt(wI^2 + wR^2)

  list(UI = UI, UR = UR, G = G)
}

# Helper to get fields from S4 or list-like parameters
get_param <- function(params, name) {
  if (isS4(params)) {
    slot(params, name)
  } else {
    params[[name]]
  }
}

## ------------------------------------------------------------
## Reference Simulation Wrapper (CRD / MaxEnt)
## ------------------------------------------------------------

simulate_reference_metrics <- function(
    Nsim, ratio, groups, N, K,
    ref_method, allocation_func,
    eta_maxent = 1
) {
  cat("simulate_reference_metrics called for", ref_method, "\n")

  # minimal param object for CRD / MaxEnt
  params <- list(N = N, ratio = ratio, K = K, groups = groups)
  class(params) <- paste0(ref_method, "Par")

  # Wrap into function(params) → seq_obj, as expected by evaluate_assignments
  if (ref_method == "maxent") {
    allocation_wrapper <- function(p) {
      allocation_func(p$N, p$ratio, p$K, p$groups, eta = eta_maxent)
    }
  } else {
    allocation_wrapper <- function(p) {
      allocation_func(p$N, p$ratio, p$K, p$groups)
    }
  }

  result <- evaluate_assignments(
    allocation_function = allocation_wrapper,
    params              = params,
    Nsim                = Nsim,
    ratio               = ratio,
    groups              = groups,
    method              = ref_method,
    eta_maxent_ref      = eta_maxent,
    plot                = FALSE
  )

  list(MPM = result$MPM, AFI = result$AFI)
}

## ------------------------------------------------------------
## Internal helpers for sequence handling
## ------------------------------------------------------------

# Extract an assignments data.frame (patient, treatment) from a sequence object
# Works for:
#   - S4 with 'assignments' slot (your custom classes),
#   - S4 randSeq-like classes with 'M' and 'groups',
#   - list with $assignments as data.frame.
.extract_assignments <- function(seq_obj, N, K, groups) {
  if (isS4(seq_obj) && "assignments" %in% slotNames(seq_obj)) {
    # Your custom S4 sequence classes
    assignments_data <- slot(seq_obj, "assignments")
    assignments <- data.frame(
      patient   = seq_len(N),
      treatment = factor(assignments_data$treatment, levels = groups)
    )
    return(assignments)
  }

  if (isS4(seq_obj) && "M" %in% slotNames(seq_obj)) {
    # Generic randSeq / rRandSeq-style object with matrix M
    M <- slot(seq_obj, "M")
    if (nrow(M) < 1L) stop("Sequence matrix M has no rows.")
    seq_vec <- as.integer(M[1, ])

    # Try to infer coding: 0..K-1 or 1..K
    if (all(seq_vec %in% 0:(K - 1L))) {
      idx <- seq_vec + 1L
    } else if (all(seq_vec %in% 1:K)) {
      idx <- seq_vec
    } else {
      stop("Unable to infer treatment coding from sequence matrix M.")
    }

    trt_labels <- groups[idx]
    assignments <- data.frame(
      patient   = seq_len(N),
      treatment = factor(trt_labels, levels = groups)
    )
    return(assignments)
  }

  # List-based interface with $assignments data.frame
  if (!isS4(seq_obj) && !is.null(seq_obj$assignments)) {
    assignments_data <- seq_obj$assignments
    assignments <- data.frame(
      patient   = seq_len(N),
      treatment = factor(assignments_data$treatment, levels = groups)
    )
    return(assignments)
  }

  stop("Cannot extract assignments: unsupported sequence object format.")
}

# Extract the conditional probability matrix P (K x N) from a sequence object.
# Priority:
#   1) S4 + getProbMatrix(seq_obj) if available (your new designs).
#   2) List with $probabilities matrix (N x K, rows = patients).
#   3) Fallback: empirical allocation proportions via calculate_pi_matrix().
.extract_prob_matrix <- function(seq_obj, assignments, N, K, groups) {
  prb_mat <- NULL

  if (isS4(seq_obj)) {
    # Try getProbMatrix if defined
    prb_mat <- tryCatch(
      getProbMatrix(seq_obj),
      error = function(e) NULL
    )
  }

  if (is.null(prb_mat) && !isS4(seq_obj) && !is.null(seq_obj$probabilities)) {
    # List-based designs: probabilities is N x K, we need K x N
    prob_mat <- seq_obj$probabilities
    if (!is.matrix(prob_mat)) {
      stop("'probabilities' must be a matrix (N x K).")
    }
    prb_mat <- t(prob_mat)
  }

  if (is.null(prb_mat)) {
    # Fallback: use empirical allocation proportions (π matrix).
    # This does NOT change behaviour for designs where a probability matrix
    # already exists; it only enables approximate FI for other designs.
    prb_mat <- calculate_pi_matrix(assignments, N, K, groups)
  }

  prb_mat
}

## ------------------------------------------------------------
## Main Evaluation Function
## ------------------------------------------------------------

evaluate_assignments <- function(allocation_function = NULL,
                                 params,
                                 Nsim,
                                 ratio,
                                 groups,
                                 method,
                                 eta_maxent_ref = 1,
                                 plot = TRUE) {
  # If no allocation_function is provided:
  # - For S4 randPar objects, default to genSeq(par, r = 1)
  # - For others (CRD, MaxEnt, custom list-based), the user must supply one.
  if (is.null(allocation_function)) {
    if (isS4(params) && methods::is(params, "randPar")) {
      allocation_function <- function(par) genSeq(par, r = 1)
    } else {
      stop(
        "allocation_function is NULL and 'params' is not a 'randPar' object.\n",
        "Please provide an allocation_function for non-randPar designs ",
        "(e.g. CRD, MaxEnt)."
      )
    }
  }
  ## -------------------------
  ## Basic setup
  ## -------------------------
  N   <- get_param(params, "N")     # trial size n
  K   <- get_param(params, "K")     # number of treatments
  rho <- ratio / sum(ratio)         # target allocation proportions

  ## -------------------------
  ## Reference designs (CRD and MaxEnt(η = 1))
  ## Only simulate them if the current design is NOT itself CRD or MaxEnt
  ## -------------------------
  if (!(method %in% c("crd", "maxent"))) {
    ref_crd <- simulate_reference_metrics(
      Nsim           = Nsim,
      ratio          = ratio,
      groups         = groups,
      N              = N,
      K              = K,
      ref_method     = "crd",
      allocation_func = crd_allocation_function
    )

    ref_maxent <- simulate_reference_metrics(
      Nsim           = Nsim,
      ratio          = ratio,
      groups         = groups,
      N              = N,
      K              = K,
      ref_method     = "maxent",
      allocation_func = maxent_allocation_function,
      eta_maxent     = eta_maxent_ref
    )
  } else {
    ref_crd   <- NULL
    ref_maxent <- NULL
  }

  ## -------------------------
  ## Containers for all simulations
  ## -------------------------
  Imb_all    <- matrix(NA_real_, nrow = Nsim, ncol = N)  # Imb_s(m)
  FI_all     <- matrix(NA_real_, nrow = Nsim, ncol = N)  # FI_s(m)
  Imb_final  <- numeric(Nsim)                            # Imb_s(n)
  FI_final   <- numeric(Nsim)                            # FI_s(n)
  assignments_list <- vector("list", Nsim)

  last_prb_mat <- NULL

  ## -------------------------
  ## Main simulation loop
  ## -------------------------
  for (sim in seq_len(Nsim)) {

    ## 1) Generate one randomization sequence for this design
    seq_obj <- allocation_function(params)

    ## 2) Extract assignments (patient, treatment)
    assignments <- .extract_assignments(seq_obj, N, K, groups)
    assignments_list[[sim]] <- assignments

    ## 3) Imbalance by step: Imb_sim(m) for m = 1..N
    Imb_m <- numeric(N)
    for (m in seq_len(N)) {
      partial <- assignments[1:m, ]
      Ni_m <- table(factor(partial$treatment, levels = groups))

      if (K == 2L && ratio[1] == ratio[2]) {
        # 2-arm, 1:1 case: signed imbalance N1(m) - N2(m)
        Imb_m[m] <- as.numeric(Ni_m[1] - Ni_m[2])
      } else {
        # General K-arm / unequal ratio: Euclidean distance
        Imb_m[m] <- sqrt(sum((Ni_m - m * rho)^2))
      }
    }
    Imb_all[sim, ] <- Imb_m
    Imb_final[sim] <- Imb_m[N]

    ## 4) Conditional probabilities P(j) and FI_sim(m)
    prb_mat <- .extract_prob_matrix(seq_obj, assignments, N, K, groups)
    last_prb_mat <- prb_mat

    rho_mat <- matrix(rho,
                      nrow = nrow(prb_mat),
                      ncol = ncol(prb_mat),
                      byrow = FALSE)

    diff_mat <- prb_mat - rho_mat           # K x N
    FI_j     <- sqrt(colSums(diff_mat^2))   # length N

    FI_m <- cumsum(FI_j) / seq_len(N)       # FI_sim(m)
    FI_all[sim, ] <- FI_m
    FI_final[sim] <- FI_m[N]
  }

  ## -------------------------
  ## Step-wise expectations: MPM(m) and AFI(m)
  ## -------------------------
  MPM_by_step <- colMeans(Imb_all)  # MPM(m)
  AFI_by_step <- colMeans(FI_all)   # AFI(m)

  ## -------------------------
  ## Design-level imbalance and randomness for trial size n
  ## -------------------------
  # Imbalance measure: MPM(n)_ξ = (1/n) Σ_{m=1}^n MPM(m)
  MPM_design <- mean(MPM_by_step)

  # Randomness measure: AFI(n)_ξ = AFI(m) at m = n
  AFI_design <- AFI_by_step[N]

  # For reporting: mean Imb(n) and mean FI(n) across simulations
  Imb_mean  <- mean(Imb_final)
  FI_n_mean <- mean(FI_final)

  ## -------------------------
  ## Unified imbalance UI, unified randomness UR, overall G
  ## -------------------------
  if (!(method %in% c("crd", "maxent"))) {
    perf <- calculate_performance_measures(
      MPM_n_xi   = MPM_design,
      MPM_n_MaxEnt = ref_maxent$MPM,
      MPM_n_CRD    = ref_crd$MPM,
      AFI_n_xi   = AFI_design,
      AFI_n_MaxEnt = ref_maxent$AFI,
      AFI_n_CRD    = ref_crd$AFI
    )
    UI <- perf$UI
    UR <- perf$UR
    G  <- perf$G
  } else {
    UI <- UR <- G <- NA_real_
  }

  ## -------------------------
  ## ASD and plotting
  ## -------------------------
  assignments_df <- do.call(rbind, assignments_list)
  assignments_df$Simulation <- rep(seq_len(Nsim), each = N)

  asd_values <- calculate_asd(assignments_df, Nsim, N, K, ratio, groups)

  if (plot) {
    # 1) Conditional allocation probabilities over time
    matplot(t(last_prb_mat), type = "l", lty = 1, col = rainbow(K),
            main = paste("Allocation probabilities for", method),
            xlab = "Sample size (m)",
            ylab = "P(Treatment)",
            lwd  = 2)
    legend("topright", legend = groups, col = rainbow(K), lty = 1,
           title = "Groups")

    # 2) Final treatment distribution
    barplot(table(assignments_df$treatment) / (N * Nsim),
            main = paste("Final treatment distribution for", method),
            col  = rainbow(K),
            names.arg = groups,
            ylab = "Proportion")

    # 3) ASD over m
    plot(seq_len(N), asd_values, type = "l", col = "blue", lwd = 2,
         main = paste("ASD for", method),
         xlab = "Sample size (m)",
         ylab = "ASD")

    # 4) MPM(m) over m
    plot(seq_len(N), MPM_by_step, type = "l", lwd = 2,
         main = paste("MPM(m) for", method),
         xlab = "Sample size (m)",
         ylab = "MPM(m)")
    abline(h = 0, lty = 3)

    # 5) AFI(m) over m
    plot(seq_len(N), AFI_by_step, type = "l", lwd = 2,
         main = paste("AFI(m) for", method),
         xlab = "Sample size (m)",
         ylab = "AFI(m)")
    abline(h = 0, lty = 3)

    # 6) UI–UR–G point (3D) using scatterplot3d
    s3d <- scatterplot3d::scatterplot3d(
      x    = UR,
      y    = UI,
      z    = G,
      xlim = c(0, 1),
      ylim = c(0, 1),
      zlim = c(0, 1),
      xlab = "Unified Randomness (UR)",
      ylab = "Unified Imbalance (UI)",
      zlab = "Overall performance (G)",
      pch  = 19,
      main = sprintf("Balance/Randomness trade-off for %s", method)
    )

    pt <- s3d$xyz.convert(UR, UI, G)

    text(
      x      = pt$x,
      y      = pt$y,
      labels = sprintf("UR = %.3f\nUI = %.3f\nG = %.3f", UR, UI, G),
      pos    = 4,
      cex    = 0.9
    )
  }

  ## -------------------------
  ## Return all relevant metrics (unchanged structure)
  ## -------------------------
  list(
    Assignments   = assignments_df,
    ASD           = asd_values,
    Imbalance     = Imb_mean,       # mean Imb(n) over simulations
    FI_n          = FI_n_mean,      # mean FI(n) over simulations
    MPM           = MPM_design,     # MPM(n)_ξ
    AFI           = AFI_design,     # AFI(n)_ξ
    MPM_by_step   = MPM_by_step,    # MPM(m), m = 1..N
    AFI_by_step   = AFI_by_step,    # AFI(m), m = 1..N
    # ARP         = colMeans(ARP_values),  # left commented as in original
    UI            = UI,
    UR            = UR,
    G             = G,
    ref_crd       = ref_crd,
    ref_maxent    = ref_maxent
  )
}
