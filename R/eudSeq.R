#' @include randSeq.R eudPar.R
NULL

###############################################
# --------------------------------------------#
# Class eudSeq (Extended Urn Design)          #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Class definition for eudSeq
# --------------------------------------------

# Representation of sequences for Extended Urn Design (EUD)
#
# @description This set of classes provides functionality for storing
# randomization sequences of the EUD design along with the
# parameters representing the design.
#
# @slot N total number of patients included in the trial.
# @slot a numeric parameter a.
# @slot b numeric parameter b (> 0).
# @slot M matrix containing randomization sequences of length N in its rows.
setClass("eudSeq",
         slots = c(a = "numeric", b = "numeric"),
         contains = "randSeq")

# --------------------------------------------
# Class definition for rEudSeq
# --------------------------------------------

# Representation of random sequences for EUD
#
# @description Stores random randomization sequences of EUD along
# with the parameters representing the design.
setClass("rEudSeq", contains = c("rRandSeq", "eudSeq"))

# --------------------------------------------
# Methods for eudSeq
# --------------------------------------------

#' @rdname getProbabilities
setMethod("getProb", signature = c(obj = "eudSeq"),
          function(obj) {
            K     <- obj@K
            ratio <- obj@ratio
            a     <- obj@a
            b     <- obj@b

            Msum <- sum(ratio)

            apply(obj@M, 1, function(randSeq, K, ratio, a, b, Msum) {
              conditionalProb <- numeric(length(randSeq))
              n_counts <- integer(K)

              for (i in 1:length(randSeq)) {
                num <- ratio +
                  a * n_counts +
                  b * ratio * ((i - 1L) - n_counts) +
                  b * n_counts * (ratio - 1)

                den <- Msum + a * (i - 1L) + b * (Msum - 1L) * (i - 1L)

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

                idx <- randSeq[i] + 1L
                conditionalProb[i] <- p[idx]
                n_counts[idx] <- n_counts[idx] + 1L
              }
              prod(conditionalProb)
            }, K = K, ratio = ratio, a = a, b = b, Msum = Msum)
          }
)

#' @rdname getDesign
setMethod("getDesign",
          signature(obj = "eudSeq"),
          function(obj) {
            paste("EUD(",
                  "a=", obj@a, ",",
                  "b=", obj@b, ",",
                  "ratio=", paste(obj@ratio, collapse = ":"),
                  ")",
                  sep = "")
          }
)
