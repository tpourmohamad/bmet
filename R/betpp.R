#' Posterior predictive based Bayesian multigroup equivalence testing
#'
#' Function provides the necessary tools to carry out Bayesian multigroup equivalence testing based on sampling of the posterior predictive distribution.
#' The function returns posterior predictive samples of future differences amongst groups.
#'
#' @param values A vector of measurements sorted in the same order as the \code{groups} variable.
#' @param groups A vector of groups labels corresponding to the individual measurements in the \code{groups} variable.
#' @param em A c x 2 matrix of lower and upper equivalence margins. Here, c is the number of pairwise comparisons of interest.
#' @param A A c x k matrix of pairwise contrasts. Here, k is the number of groups, i.e., \code{length(unique(groups))}.
#' @param B A positive integer specifying the number of posterior predictive samples to draw. By default \code{B} is set to 10000.
#' @return The function returns a list object containing the following:
#' \itemize{
#'   \item prob: The probability that future differences fall within the equivalence margins.
#'   \item delta: A B x c matrix of posterior predictive samples of future differences for each pairwise comparison of interest.
#' }
#' @references Pourmohamad, T. and Lee, H.K.H. (2023). Equivalence Testing for Multiple Groups. Stat, e645.
#' @examples
#' ### Multigroup equivalence test for A vs. B and A vs. C
#' values <- rnorm(75)
#' groups <- rep(LETTERS[1:3], each = 25)
#'
#' mad1 <- 0.65  # The equivalence margin for A vs. B
#' mad2 <- 0.65  # The equivalence margin for A vs. C
#' mads <- c(mad1, mad2)
#' mads <- cbind(-mads, mads)
#'
#' A <- apc(3)
#' A <- A[1:2, ]
#'
#' out <- betpp(values, groups, mads, A, B = 10000)
#'
#' out$prob   # The probability that future differences
#'            # fall within the equivalence margins
#'
#' @export
betpp <- function(values, groups, em, A, B = 10000){

  if(((floor(B) - B) != 0) | B < 1){
    stop("The number of posterior samples, B, must be a positive integer.")
  }else if(length(groups) != length(values)){
    stop("The length of values must be the same as the length of groups.")
  }else{
    mads <- em
    n.contrasts <- nrow(A)
    group.labels <- unique(groups)
    groups <- factor(groups, levels = group.labels)
    m <- length(group.labels)
    n <- tapply(values, groups, length)

    E.mu <- tapply(values, groups, mean)

    sigma2 <- matrix(NA, ncol = m, nrow = B)
    mu <- matrix(NA, ncol = m, nrow = B)
    pp <- matrix(NA, ncol = n.contrasts, nrow = B)

    for(i in 1:B){
      for(j in 1:m){
        x <- values[groups == group.labels[j]]
        sigma2[i, j] <- MCMCpack::rinvgamma(1, n[j] / 2 - 1 / 2, sum((x - E.mu[j])^2) / 2)
        mu[i, j] <- stats::rnorm(1, E.mu[j], sigma2[i, j] / n[j])
      }

      Cov.mu <- diag(sigma2[i,])
      pp[i,] <- MASS::mvrnorm(1, mu = A %*% mu[i,], Sigma = A %*% Cov.mu %*% t(A))

    }

    prob <- sum(rowSums(t(t(pp) >= mads[,1]) + t(t(pp) <= mads[,2])) == 2 * n.contrasts) / B
    return(list(prob = prob, ppdraws = pp))

  }

}
