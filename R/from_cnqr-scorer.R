#' Define a Scoring Rule
#'
#' Indicate which scoring rule you'd like to use. Note that the scoring rule
#' is defined by the quantile levels, the transformation function, and the
#' across-quantile weight function.
#'
#' @param tau Vector of quantile levels (in (0,1)) to score one.
#' @param g Vectorized transformation function of the response.
#' @param wtau Vector of weights to apply corresponding to each \code{tau}.
#'
#' @return A list of the following entries:
#'
#' \itemize{
#'    \item \code{$tau} The inputted \code{tau} vector.
#'    \item \code{$g} The inputted \code{g} function.
#'    \item \code{$wtau} Vector of length matching \code{tau} of weights.
#' }
#'
#' @note This function isn't especially glamorous; it's just intended as a
#' simple way to "bundle" information together so that you don't have to work
#' with as many "moving parts" in an analysis. It does, however, provide
#' a convenient framework for dealing with different ways of specifying
#' weights on \code{tau}.
#' @export
scorer <- function(tau, g=identity, wtau=1) {
    ## If wtau is a function, evaluate it at tau:
    if (is.function(wtau)) wtau <- wtau(tau)
    ## Ensure the length of wtau matches that of tau:
    wtau <- wtau * rep(1, length(tau))
    ## Return.
    list(tau=tau, g=g, wtau=wtau)
}
