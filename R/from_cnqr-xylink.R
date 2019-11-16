#' Extract information from a vine linking predictors and response
#'
#' Given a vine, this function will assume the first variable in the last
#' column is a response variable, and the variables below it are predictors.
#' This function will return the variable numbers of the predictors and
#' response, as well as the copulas linking them.
#'
#' @param obj Object of type \code{'rvine'}.
#' @return A list of the following items:
#'
#' \itemize{
#'      \item \code{$yvar}: Integer of the variable number of the response.
#'      \item \code{$xord}: Integer vector of the predictors listed below
#'      the response in the vine array (in order from upper to lower).
#'      \item \code{$p}: Integer; the number of predictors listed below the
#'      response.
#'      \item \code{$cops}: Character vector of the copula names linking the
#'      predictors and the response, in the order of \code{xord}.
#'      \item \code{$cpars}: List of numeric vectors, each vector being the
#'      copula parameters corresponding to the copulas in \code{cops}.
#' }
#' @export
xylink <- function(obj) {
    d <- ncol(obj$G)
    if (d == 0) return(NA)
    ## Y variable
    yvar <- obj$G[1, d]
    ## x variables (in order)
    xord <- obj$G[-1, d]
    xord <- xord[xord != 0]
    p <- length(xord)
    ## copulas and parameters
    cops <- obj$copmat[seq_len(p), d]
    cpars <- obj$cparmat[seq_len(p), d]
    ## Output
    list(yvar=yvar,
         xord=xord,
         p=p,
         cops=cops,
         cpars=cpars)
}
