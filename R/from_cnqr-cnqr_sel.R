#' Select Vines using CNQR
#'
#' Select from vine models using CNQR. The vines can only differ in the
#' last columns of their copula and copula parameter matrices, and are
#' assumed to have the response in the upper-right corner of each vine array.
#' The vine whose forecasts score optimally is selected.
#'
#' @param rv_list list of \code{'rvine'} objects, each of which may only differ
#' in the last columns of the copula and copula parameter matrices.
#' @param sc Scoring rule to use for the regression, as in the output
#' of \code{scorer}.
#' @param y Vector of response data.
#' @param uind Matrix of independent uniform predictors,
#' as in the output of \code{pcondseq}, of the predictors
#' in the last column of the (common) vine array. If we call "\code{b}" that
#' column without the first (response) variable, the input matrix should
#' be the PIT scores of variables
#' \code{b[1]}; \code{b[2]|b[1]}; \code{b[3]|b[1:2]}; etc.
#' @param QY Quantile function of the response \code{y}, which accepts a
#' vector of values (quantile levels) in (0,1). It should return
#' quantiles, either in the form of
#' a vector corresponding to the input,
#' or in the form of a matrix with columns corresponding to the inputted
#' quantile levels and rows corresponding to the observations of \code{y}
#' (thus allowing for each value in \code{y} to come from different
#' distributions).
#' @return Outputs the "best" entry of \code{rv_list}.
#' @note This function assumes that the inputted vines may only differ in the
#' last column of the copula and copula parameter matrices, and does not
#' check to see that this is the case.
#' @seealso \code{cnqr_est} for CNQR when the model space is a
#' parameter space.
#' @export
cnqr_sel <- function(rv_list, sc, y, uind, QY=identity) {
    if (length(rv_list) == 0) return(NA)
    if (length(rv_list) == 1) return(rv_list[[1]])
    tau <- sc$tau
    d <- ncol(rv_list[[1]]$G)
    ## Get scores on each rvine.
    scores <- sapply(rv_list, function(rv_) {
        ## Extract copulas and parameters.
        cops <- rv_$copmat[, d]
        sel <- cops != ""
        cops <- cops[sel]
        cpars <- rv_$cparmat[sel, d]
        ## Do the prediction
        yhat <- QYgX(tau, uind, cops, cpars, QY=QY)
        ## Score the predictions
        score_eval(y, yhat, sc=sc)
    })
    ## Which model provides the smallest score?
    best <- which(scores == min(scores))
    if (length(best) > 1)
        warning("More than one 'rvine' yields an optimal score. Choosing the first one.")
    rv_list[[best[1]]]
}
