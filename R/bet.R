#' Posterior based Bayesian multigroup equivalence testing
#'
#' Function provides the necessary tools to carry out Bayesian multigroup equivalence testing based on sampling of the posterior distribution.
#' The function returns posterior samples of the average differences amongst groups, as well as posterior samples of group variances.
#'
#' @param values A vector of measurements sorted in the same order as the \code{groups} variable.
#' @param groups A vector of groups labels corresponding to the individual measurements in the \code{groups} variable.
#' @param em A c x 2 matrix of lower and upper equivalence margins. Here, c is the number of pairwise comparisons of interest.
#' @param A A c x k matrix of pairwise contrasts. Here, k is the number of groups, i.e., \code{length(unique(groups))}.
#' @param B A positive integer specifying the number of posterior samples to draw. By default \code{B} is set to 10000.
#' @param test Setting this to anything other than "mean" tells the function to not calculate the posterior probability that the average
#'             differences fall within the equivalence margins (applicable when testing equivalence based on something other than just
#'             average differences).
#' @return The function returns a list object containing the following:
#' \itemize{
#'   \item prob: The posterior probability that the average differences fall within the equivalence margins. Only returned if \code{test == "mean"}.
#'   \item delta: A B x c matrix of posterior samples of the average difference for each pairwise comparison of interest.
#'   \item sigma2: A B x k matrix of posterior samples of the variance for each group.
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
#' out <- bet(values, groups, mads, A, B = 10000)
#'
#' out$prob   # The posterior probability that the average
#'            # differences fall within the equivalence margins
#'
#' @export
bet <- function(values, groups, em, A, B = 10000, test = "mean"){

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
    delta <- matrix(NA, ncol = n.contrasts, nrow = B)

    for(i in 1:B){
      for(j in 1:m){
        x <- values[groups == group.labels[j]]
        sigma2[i, j] <- MCMCpack::rinvgamma(1, n[j] / 2 - 1 / 2, sum((x - E.mu[j])^2) / 2)
      }

      Cov.mu <- diag(sigma2[i,] / n)
      delta[i,] <- MASS::mvrnorm(1, mu = A %*% E.mu, Sigma = A %*% Cov.mu %*% t(A))

    }

    colnames(delta) <- rownames(A)
    colnames(sigma2) <- group.labels

    if(test == "mean"){
      prob <- sum(rowSums(t(t(delta) > mads[,1]) + t(t(delta) < mads[,2])) == 2 * n.contrasts) / B
    }else{
      prob <- NA
    }

    return(list(prob = prob, delta = delta, sigma2 = sigma2))

  }

}
