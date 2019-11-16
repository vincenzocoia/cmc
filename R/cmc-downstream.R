#' Evaluate cmc model
#'
#' Computes things useful for downstream computations.
#' @param object A cmc-fitted model.
#' @param newdata Data frame to operate on
eval.cmc <- function(object, newdata = NULL) {
	ycol  <- object$ycol
	x1col <- object$x1col
	x2col <- object$x2col
	xvine <- object$xvine
	if (is.null(newdata)) {
		newdata <- object$data
		xdat <- as.matrix(newdata[, c(ycol, x1col, x2col)])
		udat <- apply(xdat, 2, object$marginal$pdist)
		u2cond <- object$u2cond
	} else {
		xdat <- as.matrix(newdata[, c(ycol, x1col, x2col)])
		udat <- apply(xdat, 2, object$marginal$pdist)
		u2cond <- copsupp::pcondrvine(udat, xvine, var = 3, condset = 2)
	}
	u1u2cond <- cbind(udat[, 2], u2cond)
	list(
		newdata  = newdata,
		u1u2cond = u1u2cond
	)
}

eval <- function(object, newdata = NULL) UseMethod("eval")

#' Add structure to a null model
#'
#' Converts marginal Unif(0,1) distributions to have modelled distributions.
#' CURRENTLY ONLY HANDLES PIT SCORES AS VALUES.
#'
#' @param object A cmc-fitted model.
#' @param newdata Data frame to operate on
#' @param from_col Name of column containing the null values. For example,
#' probabilities (PIT scores) to evaluate the quantile functions at.
#' @param to_col Name of the column to append the output to. Leave blank if
#' you want a vector output.
#' @details If from_col values are numbers, these are interpreted as PIT
#' scores by applying the quantile function corresponding
#' to each predictive distribution.
#' If values are distributions, these are converted in
#' such a way that the predictive distributions are obtained if the values
#' are all Unif(0,1) distributions.
#' @rdname cmc_downstream
#' @export
construct.cmc <- function(object, newdata = NULL, from_col, to_col) {
	ycol  <- object$ycol
	x1col <- object$x1col
	x2col <- object$x2col
	eval <- eval(object, newdata)
	newdata  <- eval$newdata
	u1u2cond <- eval$u1u2cond
	pits <- newdata[[from_col]]
	constructed <- cnqr::QYgX(
		pits, u1u2cond,
		cops = c(object$vine1$copmat[1, 2],
				 object$vine2$copmat[1, 2]),
		cpars = list(object$vine1$cparmat[1, 2][[1]],
					 object$vine2$cparmat[1, 2][[1]]),
		QY = identity
	)
	if (missing(to_col)) {
		return(constructed)
	} else {
		newdata[[to_col]] <- constructed
		return(newdata)
	}
}

#' @rdname cmc_downstream
#' @export
construct <- function(object, newdata = NULL, from_col, to_col) UseMethod("construct")

#' @rdname cmc_downstream
#' @export
decompose.cmc <- function(object, newdata = NULL, from_col, to_col) {
	ycol  <- object$ycol
	x1col <- object$x1col
	x2col <- object$x2col
	eval <- eval(object, newdata)
	newdata  <- eval$newdata
	u1u2cond <- eval$u1u2cond
	v <- newdata[[from_col]]
	decomposed <- cnqr::FYgX(
		v, u1u2cond,
		cops = c(object$vine1$copmat[1, 2],
				 object$vine2$copmat[1, 2]),
		cpars = list(object$vine1$cparmat[1, 2][[1]],
					 object$vine2$cparmat[1, 2][[1]]),
		FY = identity
	)
	if (missing(to_col)) {
		return(decomposed)
	} else {
		newdata[[to_col]] <- decomposed
		return(newdata)
	}
}

#' @rdname cmc_downstream
#' @export
decompose <- function(object, newdata = NULL, from_col, to_col) UseMethod("decompose")


#' @param what What to predict. Could be "pdist", "qdist", or "rdist".
#' @param at Vector of values to evaluate the \code{what} function at.
#' @param ... Not used
#' @rdname cmc_downstream
#' @export
predict.cmc <- function(object, ..., newdata = NULL, what = "pdist", at, to_col) {
	ycol  <- object$ycol
	x1col <- object$x1col
	x2col <- object$x2col
	eval  <- eval(object, newdata)
	u1u2cond <- eval$u1u2cond
	if (what == "pdist") {
		funs <- apply(u1u2cond, 1, function(row) {
			this_u1u2cond <- matrix(row, nrow = 1)
			function(v) cnqr::FYgX(
				v, this_u1u2cond,
				cops = c(object$vine1$copmat[1, 2],
						 object$vine2$copmat[1, 2]),
				cpars = list(object$vine1$cparmat[1, 2][[1]],
							 object$vine2$cparmat[1, 2][[1]]),
				FY = identity
			)
		})
	} else if (what == "qdist" | what == "rdist") {
		funs <- apply(u1u2cond, 1, function(row) {
			this_u1u2cond <- matrix(row, nrow=1)
			function(p) cnqr::QYgX(
				p, this_u1u2cond,
				cops = c(object$vine1$copmat[1, 2],
						 object$vine2$copmat[1, 2]),
				cpars = list(object$vine1$cparmat[1, 2][[1]],
							 object$vine2$cparmat[1, 2][[1]]),
				QY = identity
			)[1, ]
		})
		if (what == "rdist") {
			qdists <- funs
			funs <- lapply(qdists, function(qdist)
				function(n) qdist(stats::runif(n))
			)
		}
	} else {
		stop("Don't know what to do with '", what,
			 "' as entry of 'what' argument.")
	}
	if (missing(at)) {
		list_output <- funs
	} else {
		list_output <- lapply(funs, function(fun) {
			res <- data.frame(.at  = at,
							  .fun = fun(at))
			names(res)[2] <- paste0(".", what)
			res
		})
	}
	if (missing(to_col)) {
		return(list_output)
	} else {
		newdata <- eval$newdata
		newdata[[to_col]] <- list_output
		return(newdata)
	}
}
