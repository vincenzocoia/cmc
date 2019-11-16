#' Quantile Function for Y|X
#'
#' Computes the conditional quantile function (\code{QYgX})
#' or distribution function (\code{FYgX}) of a response
#' at specified quantile levels, given some observed predictors.
#' This is like the quantile version of the \code{pcondrvine} function
#' in the \code{copsupp} packages, except that it (1) only
#' takes the response to be the last entry in the vine, and (2) saves
#' computation time by accepting the predictors converted to an independent
#' set (thus bypassing the need to specify the distribution of the predictors).
#'
#' @param tau Vector of quantile levels to evaluate at. Or,
#' different vectors (of the same length) can be specified
#' for each observation by making this a matrix (with rows corresponding
#' to observations).
#' @param ucond Matrix of uniform independent predictors. The i'th
#' column corresponds to the PIT score of the i'th predictor linked to the
#' response, conditional on previously linked predictors. See
#' \code{pcondseq} in the \code{copsupp} package for this.
#' @param cops Vector of copula model names, corresponding to the columns of
#' \code{ucond}, that link (predictor, response) conditional on previous
#' predictors in the pairing order.
#' @param cpars List of parameter vectors corresponding to \code{cops}.
#' @param QY,FY Vectorized (marginal) quantile function/cdf of the response.
#' @return A matrix, with columns being the quantile levels, and rows
#' corresponding to the observations in \code{ucond}.
#' @note If copulas are not permutation-symmetric, ensure that the copula
#' specified has its fist argument corresponding to the predictor, and second
#' argument the response.
#' @import CopulaModel
#' @examples
#' (dat <- matrix(runif(3*6), ncol = 3))
#' tau <- c(0.1, 0.3, 0.7, 0.9)
#' cops <- c("gum", "bvtcop", "frk")
#' cpars <- list(3.5, c(0.6, 3), 2.5)
#'
#' ## Get tau quantiles for the predictors in `dat`
#' QYgX(tau, dat, cops = cops, cpars = cpars, QY = qexp)
#'
#' ## Get sample for these predictors:
#' (v <- runif(nrow(dat)))
#' (y <- QYgX(matrix(v), dat, cops = cops, cpars = cpars, QY = qexp))
#' FYgX(y, dat, cops, cpars, FY = pexp)  # Same as v.
#' @rdname dist_YgX
#' @export
QYgX <- function(tau, ucond, cops, cpars, QY = identity) {
	## We'll evaluate the quantile function recursively. Starting
	##  at the "tail" of ucond, we'll gradually "modify" tau
	##  according to the copula models.
	n <- nrow(ucond)
	d <- ncol(ucond)
	if (is.vector(tau)) tau <- matrix(tau, nrow = n, ncol = length(tau), byrow = TRUE)
	if (d == 0) {
		## There are no more predictors, so there's no need to modify tau.
		## Evaluate at the marginal.
		res <- apply(tau, 2, QY)
		if (nrow(ucond) == 1) res <- matrix(res, nrow=1)
		return(res)
	}
	## There are still some predictors. Use the last one to modify tau.
	u <- ucond[, d]
	## Get copula conditional quantile function. Assumes the second variable
	##  linked is the response.
	qcond <- utils::getFromNamespace(paste0("qcond", cops[d]), "CopulaModel")
	qcondfit <- function(p) qcond(p, u, cpars[[d]])
	taunew <- apply(tau, 2, qcondfit)
	## Remove a data column and evaluate the quantile function at taunew
	ucondnew <- ucond[, -d, drop=FALSE]
	return(QYgX(taunew, ucondnew, cops, cpars, QY))
}

#' @param y Vector of values to evaluate the cdf at.
#' @rdname dist_YgX
#' @export
FYgX <- function(y, ucond, cops, cpars, FY = identity) {
	d <- ncol(ucond)
	if (d == 0) {
		return(FY(y))
	}
	pcondcop <- utils::getFromNamespace(paste0("pcond", cops[d]), "CopulaModel")
	pcondcop(FYgX(y, ucond[, -d, drop=FALSE], cops, cpars, FY),
			 ucond[, d],
			 cpars[[d]])
}
