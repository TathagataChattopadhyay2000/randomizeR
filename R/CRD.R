# ------------------------------------------------------------------
# Completely Randomized Design (CRD)
# ------------------------------------------------------------------

# Internal helper: allocation probabilities for CRD
# rnd: list(target = ratio)
.allocation_prb_crd <- function(rnd) {
  w   <- rnd$target
  w / sum(w)
}

# Allocation function for CRD
#
# Arguments:
#   N      : total number of subjects
#   ratio  : target allocation counts, length K
#   K      : number of arms
#   groups : group labels (length K)
#
# Returns:
#   A list with components:
#     - assignments   : data.frame(treatment = ...)
#     - probabilities : N x K matrix with constant row = target probabilities
#   and class "crdPar" (just a tag, not S4).
# ------------------------------------------------------------------
crd_allocation_function <- function(N, ratio, K, groups) {

  probs <- .allocation_prb_crd(list(target = ratio))

  # Assignments
  assignments <- sample(groups, N, replace = TRUE, prob = probs)

  # Constant probability matrix for FI_n:
  # each row m has the same vector of probabilities
  prob_matrix <- matrix(
    rep(probs, each = N),
    nrow = N,
    ncol = K,
    byrow = FALSE,
    dimnames = list(NULL, groups)
  )

  structure(
    list(
      assignments   = data.frame(treatment = assignments),
      probabilities = prob_matrix
    ),
    class = "crdPar"
  )
}
