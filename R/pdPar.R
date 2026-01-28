#' @include randPar.R
NULL

###############################################
# --------------------------------------------#
# Class pdPar (Proportional Design)           #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Validity check
# --------------------------------------------

# Validity check function for objects of the pdPar class
#
# @inheritParams overview
#
# @return Returns TRUE if the settings of the object are valid.
validatepdpar <- function(object) {
  errors <- character()

  N      <- object@N
  K      <- object@K
  ratio  <- object@ratio
  pHigh  <- object@pHigh
  groups <- object@groups

  ## N: total sample size
  if (length(N) != 1 || !is.finite(N) || N <= 0 || round(N) != N) {
    errors <- c(errors, "'N' must be a single positive integer.")
  }

  ## K: number of arms
  if (length(K) != 1 || !is.finite(K) || K <= 0 || round(K) != K) {
    errors <- c(errors, "'K' must be a single positive integer.")
  }

  ## ratio vs K
  if (length(ratio) != K) {
    errors <- c(errors, "Length of 'ratio' must equal 'K'.")
  }
  if (any(!is.finite(ratio)) || any(ratio < 0)) {
    errors <- c(errors, "All elements of 'ratio' must be finite and non-negative.")
  }
  if (sum(ratio) == 0) {
    errors <- c(errors, "Sum of 'ratio' must be greater than zero.")
  }

  ## groups
  if (length(groups) != K) {
    errors <- c(errors, "Length of 'groups' must equal 'K'.")
  }
  if (any(is.na(groups)) || any(groups == "")) {
    errors <- c(errors, "'groups' must be non-empty character labels.")
  }
  if (any(duplicated(groups))) {
    errors <- c(errors, "'groups' must be unique.")
  }

  ## pHigh: vector of high probabilities
  if (!is.numeric(pHigh) || any(!is.finite(pHigh))) {
    errors <- c(errors, "'pHigh' must be a numeric vector of finite values.")
  } else {
    if (length(pHigh) != K) {
      errors <- c(errors, "Length of 'pHigh' must equal 'K'.")
    } else {
      q <- ratio / sum(ratio)
      if (any(pHigh <= q)) {
        errors <- c(errors,
                    "'pHigh' must be elementwise strictly greater than target proportions 'ratio/sum(ratio)'.")
      }
      if (any(pHigh >= 1)) {
        errors <- c(errors, "All elements of 'pHigh' must be strictly less than 1.")
      }
    }
  }

  if (length(errors) == 0) TRUE else errors
}

# --------------------------------------------
# Class definition for pdPar
# --------------------------------------------

# Randomization parameters for Proportional Design
setClass("pdPar",
         slots = c(pHigh = "numeric"),
         contains = "randPar",
         validity = validatepdpar)

# --------------------------------------------
# Constructor for pdPar
# --------------------------------------------

#' Parameters for Proportional Design
#'
#' Represents the parameter settings of the proportional design
#' (generalized Efron's biased coin in proportional form).
#'
#' @details
#' The design is defined for \code{K} treatment groups with target
#' allocation elements \code{ratio = (r_1, ..., r_K)} and a vector
#' of "high" probabilities \code{pHigh} used when a uniquely most
#' underrepresented arm is present.
#'
#' If \code{pHigh} is not supplied, the default
#' \deqn{pHigh(k) = (q_k + 1)/2, \quad q_k = r_k / \sum_j r_j}
#' is used.
#'
#' @family randomization procedures
#'
#' @inheritParams overview
#' @param ratio numeric vector of target allocation elements \eqn{r_j \ge 0}.
#' @param pHigh optional numeric vector of length \code{K} with high
#'   probabilities \eqn{pHigh_j}, each strictly between \eqn{q_j} and 1.
#'
#' @return
#' \code{S4} object of class \code{pdPar}.
#'
#' @export
pdPar <- function(N, K, ratio, pHigh = NULL, groups = LETTERS[1:K]) {
  if (length(ratio) != K) {
    stop("Length of 'ratio' must equal 'K'.")
  }
  if (length(groups) != K) {
    stop("Length of 'groups' must equal 'K'.")
  }

  if (is.null(pHigh)) {
    q <- ratio / sum(ratio)
    pHigh <- (q + 1) / 2
  }

  new("pdPar",
      N      = N,
      K      = K,
      ratio  = ratio,
      pHigh  = pHigh,
      groups = groups)
}

# --------------------------------------------
# Core sampling algorithm (Proportional Design)
# --------------------------------------------

# internal helper: generates a single numeric sequence 0..K-1
pdRand <- function(N, K, ratio, pHigh) {

  ## target proportions q
  q <- ratio / sum(ratio)

  alloc  <- integer(K)  # counts per arm
  result <- integer(N)  # sequence 0..K-1

  for (i in seq_len(N)) {
    n_current <- sum(alloc)

    if (n_current == 0L) {
      ## First patient: use target ratio
      probabilities <- q

    } else {
      ## Imbalance vector d = N_n - n * q
      d       <- alloc - n_current * q
      d_min   <- min(d)
      idx_min <- which(d == d_min)
      Ties    <- length(idx_min)

      if (Ties == K) {
        ## all d_k equal -> use target ratio
        probabilities <- q

      } else if (Ties == 1L) {
        ## unique most underrepresented arm k0
        k0 <- idx_min

        probabilities <- numeric(K)
        ## high probability for the most favorable arm
        probabilities[k0] <- pHigh[k0]

        ## remaining probability distributed among others in proportion to q
        others <- setdiff(seq_len(K), k0)
        if (length(others) > 0L) {
          q_others <- q[others]
          q_others <- q_others / sum(q_others)
          probabilities[others] <- (1 - pHigh[k0]) * q_others
        }

      } else {
        ## 1 < Ties < K : equal probability among tied arms
        probabilities <- numeric(K)
        probabilities[idx_min] <- 1 / Ties
      }
    }

    trt <- sample.int(K, size = 1L, prob = probabilities)
    result[i] <- trt - 1L
    alloc[trt] <- alloc[trt] + 1L
  }

  result
}

# --------------------------------------------
# Methods for pdPar
# --------------------------------------------

#' @rdname generateAllSequences
setMethod("getAllSeq", signature(obj = "pdPar"),
          function(obj) {
            M <- compltSet(obj)
            new("pdSeq",
                M      = M,
                N      = N(obj),
                pHigh  = obj@pHigh,
                ratio  = obj@ratio,
                groups = obj@groups,
                K      = K(obj))
          })

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "pdPar", r = "numeric", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            new("rPdSeq",
                M = t(sapply(1:r, function(x) {
                  pdRand(N(obj), K(obj), obj@ratio, obj@pHigh)
                })),
                N      = N(obj),
                pHigh  = obj@pHigh,
                K      = K(obj),
                ratio  = obj@ratio,
                groups = obj@groups,
                seed   = seed)
          })

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "pdPar", r = "missing", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            new("rPdSeq",
                M = t(pdRand(N(obj), K(obj), obj@ratio, obj@pHigh)),
                N      = N(obj),
                pHigh  = obj@pHigh,
                K      = K(obj),
                ratio  = obj@ratio,
                groups = obj@groups,
                seed   = seed)
          })

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "pdPar", r = "numeric", seed = "missing"),
          function(obj, r, seed) {
            seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            new("rPdSeq",
                M = t(sapply(1:r, function(x) {
                  pdRand(N(obj), K(obj), obj@ratio, obj@pHigh)
                })),
                N      = N(obj),
                pHigh  = obj@pHigh,
                K      = K(obj),
                ratio  = obj@ratio,
                groups = obj@groups,
                seed   = seed)
          })

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "pdPar", r = "missing", seed = "missing"),
          function(obj, r, seed) {
            seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            new("rPdSeq",
                M = t(pdRand(N(obj), K(obj), obj@ratio, obj@pHigh)),
                N      = N(obj),
                pHigh  = obj@pHigh,
                K      = K(obj),
                ratio  = obj@ratio,
                groups = obj@groups,
                seed   = seed)
          })

#' @rdname getDesign
setMethod("getDesign",
          signature(obj = "pdPar"),
          function(obj) {
            paste("PD(",
                  "ratio=", paste(obj@ratio, collapse = ":"),
                  ")",
                  sep = "")
          })
