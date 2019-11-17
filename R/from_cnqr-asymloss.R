#' Asymmetric absolute deviation function
#'
#' The loss function that defines the tau^th-quantile;
#' \code{asymloss(tau, u) = (tau - I(u<0))*u}. Vectorised over \code{u}.
#' Function \code{rho()} is deprecated; use \code{asymloss()} instead
#' (they're identical).
#'
#' @param tau Number in [0,1] representing the probability/order
#' of the quantile.
#' @param u Value to evaluate the loss function at. Vectorised.
#' @rdname asymloss
#' @export
asymloss <- function(tau, u){
    (tau - (u<0)) * u
    # ## Account for the possibility that tau = 0 or 1.
    # coeff <- (tau - (u<0))
    # nonzero <- coeff != 0
    # coeff[nonzero] <- (coeff * u)[nonzero]
    # return(coeff)
}

