#' @include randSeq.R pdPar.R
NULL

###############################################
# --------------------------------------------#
# Class pdSeq (Proportional Design)           #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Class definition for pdSeq
# --------------------------------------------

# Representation of sequences for Proportional Design
setClass("pdSeq",
         slots = c(pHigh = "numeric"),
         contains = "randSeq")

# --------------------------------------------
# Class definition for rPdSeq
# --------------------------------------------

# Representation of random sequences for Proportional Design
setClass("rPdSeq", contains = c("rRandSeq", "pdSeq"))

# --------------------------------------------
# Methods for pdSeq
# --------------------------------------------

#' @rdname getProbabilities
setMethod("getProb", signature = c(obj = "pdSeq"),
          function(obj) {
            K      <- obj@K
            ratio  <- obj@ratio
            pHigh  <- obj@pHigh

            q <- ratio / sum(ratio)

            apply(obj@M, 1, function(randSeq, K, q, pHigh) {
              conditionalProb <- numeric(length(randSeq))
              alloc <- integer(K)

              for (i in seq_along(randSeq)) {
                n_current <- sum(alloc)

                if (n_current == 0L) {
                  probabilities <- q

                } else {
                  d       <- alloc - n_current * q
                  d_min   <- min(d)
                  idx_min <- which(d == d_min)
                  Ties    <- length(idx_min)

                  if (Ties == K) {
                    probabilities <- q

                  } else if (Ties == 1L) {
                    k0 <- idx_min

                    probabilities <- numeric(K)
                    probabilities[k0] <- pHigh[k0]

                    others <- setdiff(seq_len(K), k0)
                    if (length(others) > 0L) {
                      q_others <- q[others]
                      q_others <- q_others / sum(q_others)
                      probabilities[others] <- (1 - pHigh[k0]) * q_others
                    }

                  } else {
                    probabilities <- numeric(K)
                    probabilities[idx_min] <- 1 / Ties
                  }
                }

                idx <- randSeq[i] + 1L
                conditionalProb[i] <- probabilities[idx]
                alloc[idx] <- alloc[idx] + 1L
              }

              prod(conditionalProb)
            }, K = K, q = q, pHigh = pHigh)
          })

#' @rdname getDesign
setMethod("getDesign",
          signature(obj = "pdSeq"),
          function(obj) {
            paste("PD(",
                  "ratio=", paste(obj@ratio, collapse = ":"),
                  ")",
                  sep = "")
          })
