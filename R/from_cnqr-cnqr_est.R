#' Estimate Vine Parameters using CNQR
#'
#' Estimate copula parameters using CNQR. More specifically, for a given
#' object of type \code{'rvine'} that you would like to extend by adding
#' more rows to its last column, estimates the parameters for the
#' specified copula families of each edge.
#'
#' @param rv Object of type \code{'rvine'}, to be augmented.
#' @param a Vector of variable numbers (of the predictors) to augment to the
#' end of the array of \code{rv} (i.e., the last column).
#' @param cop Character vector of copulas to augment \code{rv} with
#' (with reflections already chosen). Each entry corresponds to a new
#' edge/row in the vine (corresponding to \code{a}).
#' @param cpar_init Starting values for the copula families in \code{cop}.
#' Should be a list with entries being parameter vectors.
#' @param sc Scoring rule to use for the regression, as in the output
#' of \code{\link{scorer}}.
#' @param y Vector of response data.
#' @param uind Matrix of independent uniform predictors,
#' as in the output of \code{\link{pcondseq}}, of the predictors
#' up to the last entry in \code{a}.
#' So, if "\code{b}" is the non-zero entries of the last column of the array
#' in \code{rv} after removing the first variable (the response) and
#' appending \code{a}, the matrix will be the PIT scores of variables
#' \code{b[1]}; \code{b[2]|b[1]}; \code{b[3]|b[1:2]}; etc.
#' @param QY Quantile function of the response \code{y}, which accepts a
#' vector of values (quantile levels) in (0,1). It should return
#' quantiles, either in the form of
#' a vector corresponding to the input,
#' or in the form of a matrix with columns corresponding to the inputted
#' quantile levels and rows corresponding to the observations of \code{y}
#' (thus allowing for each value in \code{y} to come from different
#' distributions).
#' @param verbose verbose?
#' @return Returns a list of copula parameters (the same structure
#' as \code{cpar_init}), estimated with CNQR.
#' @details
#' Here's how the estimation is done.
#'
#' \enumerate{
#'      \item You begin with a vine \code{rv}, where the top-right corner
#'      of the array is the response variable. There may or may not be
#'      variables (representing the predictors) underneath this response. Note
#'      that, to be able to condition the response on predictors, those
#'      predictors must appear below the response variable in the vine array.
#'      \item Your objective is to add more variables (predictors, already
#'      existing in the vine) to the last column of the vine array. Package
#'      these additional variables in the argument \code{a}, with
#'      corresponding copula families \code{cop}.
#'      \item The copula parameters of \code{cop} need estimation. This is done
#'      by optimizing quantile predictions of the response,
#'      conditional on the predictors in \code{a} AND the predictors
#'      already listed below the response in the vine array of \code{rv}.
#'      \item The "optimization" in the previous step refers to optimizing a
#'      scoring rule, which is indicated in the \code{sc} argument.
#' }
#' @seealso \code{\link{cnqr_sel}} for CNQR when the model space is a
#' finite selection of vines.
#' @export
cnqr_est <- function(rv, a, cop, cpar_init, sc, y, uind, QY=identity, verbose=FALSE) {
    if (length(a) == 0) return(list())
    ## Initial parameters:
    len <- sapply(cpar_init, length)
    cparvec_init <- c(cpar_init, recursive=TRUE)
    if (length(cparvec_init) == 0) return(cpar_init)
    ## Get the already-fitted copulas.
    d <- ncol(rv$G)
    rows <- rv$copmat[, d] != ""
    cops_prev <- rv$copmat[rows, d]
    if (length(cops_prev) == 0) {
        cparvec_prev <- numeric(0)
        len_prev <- integer(0)
    } else {
        cpars_prev <- rv$cparmat[rows, d]
        cparvec_prev <- c(cpars_prev, recursive=TRUE)
        len_prev <- sapply(cpars_prev, length)
    }
    ## Some combined quantities
    len_all <- c(len_prev, len)
    cops_all <- c(cops_prev, cop)
    ## Get forecasts as a function of the model space parameters. Input a vector.
    tau <- sc$tau
    yhat <- function(cparvec) {
        cparvec_all <- c(cparvec_prev, cparvec)
        cpars_all <- cparvec2cpar(cparvec_all, len_all)
        QYgX(tau, uind, cops=cops_all, cpars=cpars_all, QY=QY)
    }
    ## Get parameter space, and adjust starting values is necessary.
    ## (Don't want starting value to be less than 0.01 units away from a
    ##  boundary. For example, with 'joe' copula, nlm() has problems
    ##  starting with a parameter value of 1.000001.)
    cparspace_ <- copsupp::cparspace(cop)
    bnds <- copsupp::cparspace(cop, fn = FALSE)
    cparvec_init <- pmax(bnds$lower + 0.01, cparvec_init)
    cparvec_init <- pmin(bnds$upper - 0.01, cparvec_init)
    ## Minimize the score on the data
    if (verbose) {
        obj <- function(cparvec) {
            cat("|")
            score_eval(y, yhat(cparvec), sc=sc)
        }
    } else {
        obj <- function(cparvec)
            score_eval(y, yhat(cparvec), sc=sc)
    }
    cparvec_hat <- try(copsupp::rnlm(obj, cparvec_init, cparspace_)$estimate)
    if (inherits(cparvec_hat, "try-error")) {
        warning(paste("Ignoring the error that 'nlm' threw, and just using",
                      "the starting value."))
        cparvec_hat <- cparvec_init
    }
    ## Put the resulting parameter estimate back into list form:
    copsupp::cparvec2cpar(cparvec_hat, len)
}
