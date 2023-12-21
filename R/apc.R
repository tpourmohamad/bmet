#' All pairwise comparisons
#'
#' Function creates a contrast matrix for all pairwise comparisons
#'
#' @param ngroups A positive integer greater than 1 denoting the number of groups
#' @param labs A vector of groups labels with length equal to \code{ngroups}. The default is set to \code{NULL}, and if used, the labels will be set to \code{1:length(ngroups)}.
#' @return The function returns a matrix of all pairwise contrasts.
#' @examples
#' ### A contrast matrix based on all pairwise contrasts of 5 groups
#' apc(5)
#'
#' ### Adding group labels
#' apc(5, labs = paste("Group", 1:5, sep = " "))
#'
#' @export
apc <- function(ngroups, labs = NULL) {
  if(((floor(ngroups) - ngroups) != 0) | ngroups < 2){
    stop("The number of groups must be a positive integer greater than 1.")
  }else if(length(labs) != ngroups & !is.null(labs)){
    stop("The number of groups labels must be the same as the number of groups.")
  }else{
    lev <- labs
    lfm <- diag(ngroups)
    nlev <- nrow(lfm)
    rn <- rownames(lfm)
    a <- attr(lfm, "grid")
    if(is.null(lev)){
      if(!is.null(a)){
        lev <- apply(a, 1, paste, collapse = ":")
      }else if(!is.null(rn)){
        lev <- rn
      }else{
        lev <- as.character(1:nlev)
      }
    }
    cbn <- utils::combn(seq_along(lev), 2)
    M <- lfm[cbn[1, ], ] - lfm[cbn[2, ], ]
    if(is.vector(M)){
      dim(M) <- c(1, length(M))
    }
    rownames(M) <- paste(lev[cbn[1, ]], lev[cbn[2, ]], sep = "-")
    return(M)
  }
}


