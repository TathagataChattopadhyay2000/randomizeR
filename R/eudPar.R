#' @include randPar.R
NULL

###############################################
# --------------------------------------------#
# Class eudPar (Extended Urn Design)          #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the eudPar class
#
# @inheritParams overview
#
# @return Returns TRUE if the settings of the object are valid.
validateeudpar <- function(object) {
  errors <- character()

  a     <- object@a
  b     <- object@b
  ratio <- object@ratio
  K     <- object@K

  ## a: any finite numeric (length one)
  if (length(a) != 1 || !is.numeric(a) || !is.finite(a[1])) {
    errors <- c(errors, "a must be a single finite numeric value.")
  }

  ## b: > 0, length one
  if (length(b) != 1 || !is.numeric(b) || !is.finite(b[1]) || b[1] <= 0) {
    errors <- c(errors, "b must be a single finite numeric value > 0.")
  }

  ## ratio: allocation elements r_j, sum(ratio) > 1 for EUD
  if (!is.numeric(ratio) || any(!is.finite(ratio))) {
    errors <- c(errors, "ratio must be a numeric vector of finite values.")
  } else {
    if (length(ratio) != K) {
      msg <- paste("length(ratio) is ", length(ratio),
                   ". Should be equal to K = ", K, ".", sep = "")
      errors <- c(errors, msg)
    }
    if (any(ratio <= 0)) {
      errors <- c(errors, "All entries of ratio must be > 0 (allocation elements r_j).")
    }
    if (sum(ratio) <= 1) {
      errors <- c(errors, "For EUD, sum(ratio) = M must be > 1.")
    }
  }

  if (length(errors) == 0) TRUE else errors
}

# --------------------------------------------
# Class definition for eudPar
# --------------------------------------------

# Randomization parameters for extended Urn Design
setClass("eudPar",
         slots = c(a = "numeric", b = "numeric"),
         contains = "randPar",
         validity = validateeudpar)

# --------------------------------------------
# Constructor function for eudPar
# --------------------------------------------

#' Parameters for Extended Urn Design
#'
#' Represents the parameter settings of the extended Urn Design with
#' parameters \code{a} and \code{b > 0} and allocation elements \code{ratio}.
#'
#' @details
#' The design is defined for \code{K} treatment groups with allocation
#' elements \code{ratio = (r_1, ..., r_K)} and parameters \code{a} and
#' \code{b > 0}. The total sample size is given by \code{N}.
#'
#' @family randomization procedures
#'
#' @inheritParams overview
#' @param ratio numeric vector of allocation elements \eqn{r_j}.
#' @param a numeric parameter \eqn{a}.
#' @param b numeric parameter \eqn{b > 0}.
#'
#' @return
#' \code{S4} object of class \code{eudPar}.
#'
#' @export
eudPar <- function(N, ratio, a, b, groups = LETTERS[1:length(ratio)]) {
  if (missing(ratio)) stop("'ratio' must be provided, e.g. c(3, 1, 1).")
  if (missing(a))     stop("'a' must be provided.")
  if (missing(b))     stop("'b' must be provided.")

  K <- length(ratio)
  if (length(groups) != K) {
    stop("Length of 'groups' must equal length(ratio).")
  }

  new("eudPar",
      N      = N,
      K      = K,
      ratio  = ratio,
      a      = a,
      b      = b,
      groups = groups)
}

# --------------------------------------------
# Sampling algorithm for EUD (core sampler)
# --------------------------------------------

# internal helper (numeric codes 0..K-1)
eudRand <- function(N, K, ratio, a, b) {
  M <- sum(ratio)
  if (M <= 1) stop("For EUD, sum(ratio) must be > 1.")
  if (b <= 0) stop("For EUD, 'b' must be > 0.")

  R <- integer(N)
  n_counts <- integer(K)

  for (i in 1:N) {
    ## num_j = r_j + a*n_{i-1,j} + b*r_j*(i-1 - n_{i-1,j}) + b*n_{i-1,j}*(r_j - 1)
    num <- ratio +
      a * n_counts +
      b * ratio * ((i - 1L) - n_counts) +
      b * n_counts * (ratio - 1)

    ## den = M + a*(i-1) + b*(M - 1)*(i-1)
    den <- M + a * (i - 1L) + b * (M - 1L) * (i - 1L)

    num[num < 0] <- 0
    if (den <= 0) {
      p <- rep(1 / K, K)
    } else {
      p <- num / den
    }
    p[p < 0] <- 0

    s <- sum(p)
    if (s <= 0) {
      p <- rep(1 / K, K)
    } else {
      p <- p / s
    }

    trt_index <- sample.int(K, size = 1L, prob = p)
    R[i] <- trt_index - 1L
    n_counts[trt_index] <- n_counts[trt_index] + 1L
  }

  R
}

# --------------------------------------------
# Methods for eudPar
# --------------------------------------------

#' @rdname generateAllSequences
setMethod("getAllSeq", signature(obj = "eudPar"),
          function(obj) {
            M <- compltSet(obj)
            new("eudSeq",
                M      = M,
                N      = N(obj),
                a      = obj@a,
                b      = obj@b,
                ratio  = obj@ratio,
                groups = obj@groups,
                K      = K(obj))
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "eudPar", r = "numeric", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            new("rEudSeq",
                M = t(sapply(1:r, function(x) {
                  eudRand(N(obj), K(obj), obj@ratio, obj@a, obj@b)
                })),
                N      = N(obj),
                a      = obj@a,
                b      = obj@b,
                K      = K(obj),
                ratio  = obj@ratio,
                groups = obj@groups,
                seed   = seed)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "eudPar", r = "missing", seed = "numeric"),
          function(obj, r, seed) {
            set.seed(seed)
            new("rEudSeq",
                M = t(eudRand(N(obj), K(obj), obj@ratio, obj@a, obj@b)),
                N      = N(obj),
                a      = obj@a,
                b      = obj@b,
                K      = K(obj),
                ratio  = obj@ratio,
                groups = obj@groups,
                seed   = seed)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "eudPar", r = "numeric", seed = "missing"),
          function(obj, r, seed) {
            seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            new("rEudSeq",
                M = t(sapply(1:r, function(x) {
                  eudRand(N(obj), K(obj), obj@ratio, obj@a, obj@b)
                })),
                N      = N(obj),
                a      = obj@a,
                b      = obj@b,
                K      = K(obj),
                ratio  = obj@ratio,
                groups = obj@groups,
                seed   = seed)
          }
)

#' @rdname generateRandomSequences
setMethod("genSeq", signature(obj = "eudPar", r = "missing", seed = "missing"),
          function(obj, r, seed) {
            seed <- sample(.Machine$integer.max, 1)
            set.seed(seed)
            new("rEudSeq",
                M = t(eudRand(N(obj), K(obj), obj@ratio, obj@a, obj@b)),
                N      = N(obj),
                a      = obj@a,
                b      = obj@b,
                K      = K(obj),
                ratio  = obj@ratio,
                groups = obj@groups,
                seed   = seed)
          }
)

#' @rdname getDesign
setMethod("getDesign",
          signature(obj = "eudPar"),
          function(obj) {
            paste("EUD(",
                  "a=", obj@a, ",",
                  "b=", obj@b, ",",
                  "ratio=", paste(obj@ratio, collapse = ":"),
                  ")",
                  sep = "")
          }
)
