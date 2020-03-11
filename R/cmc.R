#' Fit CMC model via CNQR or MLE
#'
#' Provides a raw workhorse for CNQR (compute_cmc_cnqr) or MLE (compute_cmc_mle).
#' Need to provide it with highly processed input. Intended to be used
#' as internal functions. The "_ig" version forces (U1, V) to have an IG copula.
#'
#' @param v Vector of PIT scores of the response
#' @param u1 Vector of PIT scores of the first predictor to link to the response.
#' @param u2cond Vector of PIT scores of the second predictor, conditional on
#' the first.
#' @rdname compute_cmc_rd
#' @return Two vines: vine1 is the copula linking (X1, Y), and
#' vine2 is the copula linking (X2, Y)|X1.
compute_cmc_mle <- function(v, u1, u2cond) {
	u1v <- cbind(u1, v)
	vine1 <- copsupp::fitrvine_basic(u1v)
    vcond <- copsupp::pcondrvine(u1v, vine1, var = 2, condset = 1)
    u2vcond <- cbind(u2cond, vcond)
    vine2 <- copsupp::fitrvine_basic(u2vcond)
    list(
        vine1 = vine1,
        vine2 = vine2
    )
}

#' @rdname compute_cmc_rd
#' @param force_ig TRUE if you want
#' @param xvine Copula joining the 1st and 2nd columns of dmat.
#' @param method Passed to \code{fit_igcop_mle}.
compute_cmc_mle_igcop <- function(v, u1, u2cond, method) {
	cparigcop <- fit_igcop_mle(u1, v, method = method)
	cparmat <- copsupp::makevinemat(list(cparigcop), zerocol = TRUE)
	vine1 <- copsupp::rvine(
		G = matrix(c(1, 2,
					 0, 1), byrow = TRUE, ncol = 2),
		copmat = matrix(c("", "igcop"), nrow = 1),
		cparmat = cparmat
	)
	vcond <- copsupp::pcondigcop(v, u1, cparigcop)
	u2vcond <- cbind(u2cond, vcond)
	vine2 <- copsupp::fitrvine_basic(u2vcond)
	list(
		vine1 = vine1,
		vine2 = vine2
	)
}

#' @param verbose Passed to \code{cnqr_reduced}
#' @param copspace Passed to \code{cnqr_reduced}. \code{force_ig}
#' takes precedence.
#' @param sc Scorer, as in the output of \code{scorer}
#' @param families Vector of copula family names acting as a "pool"
#' to choose from when fitting.
#' @rdname compute_cmc_rd
compute_cmc_cnqr <- function(udat, force_ig, xvine, u2cond, sc,
						 verbose, copspace = NULL,
						 families = c("indepcop", "bvncop","bvtcop",
						 			 "mtcj","gum",
						 			 "frk","joe","bb1", "bskewncop",
						 			 "bskewncopp")) {
	if (force_ig) {
        if (is.null(copspace)) {
            copspace <- list("igcop", NULL)
        } else {
            copspace[[1]] <- "igcop"
        }
    }
    fit <- cnqr_reduced(1:3, dat = udat, sc = sc, basevine = xvine,
                        copspace = copspace, verbose = verbose,
    					families = families)
    vine1 <- copsupp::rvine(
        G       = matrix(ncol = 2, nrow = 2),
        copmat  = copsupp::makevinemat("", fit$copmat[1,3]),
        cparmat = copsupp::makevinemat(NULL, fit$cparmat[1,3])
    )
    vine2 <- copsupp::rvine(
        G       = matrix(ncol = 2, nrow = 2),
        copmat  = copsupp::makevinemat("", fit$copmat[2,3]),
        cparmat = copsupp::makevinemat(NULL, fit$cparmat[2,3])
    )
    list(
        vine1 = vine1,
        vine2 = vine2
    )
}


#' Fit a CMC model
#'
#' Intended to be used after a marginal model has been fit with
#' composite_dist(). So, the actual fitting that's done here is
#' strictly copula fitting, but the output is the entire model.
#' Limited capability here -- two predictors are assumed to be
#' lags, and therefore having the same marginals.
#'
#' @param ycol,x1col,x2col Character names of the columns.
#' @param data Data frame of data
#' @param force_ig Force Y and X1 to have dependence described by an IG copula?
#' Default is TRUE.
#' @param marginal Marginal model of class "dst" from the distplyr package.
#' Distribution is assumed to be the same for the response and
#' predictors. If NULL, assumes variables are already PIT scores.
#' @param method Method of fitting. For MLE (the default), one of
#' "optim" or "nlm", depending
#' on the numerical optimizer you want to use. If "cnqr",
#' fits by CNQR, and you need
#' to specify the sc argument.
#' @param sc Scorer obtained through scorer()
#' @param verbose Only works for CNQR. If TRUE, will output the fitting
#' status of CNQR.
#' @param copspace Only works for CNQR.
#' @param omit_na Omit rows of data where at least one of ycol, x1col,
#' and x2col are NA? Currently only TRUE works.
#' @param families Vector of copula family names acting as a "pool"
#' to choose from when fitting.
#' @export
cmc <- function(ycol, x1col, x2col, data, method = "optim",
                force_ig = TRUE, marginal = NULL, sc, verbose = FALSE,
                copspace = NULL, omit_na = TRUE,
				families = c("indepcop", "bvncop","bvtcop","mtcj","gum",
							 "frk","joe","bb1", "bskewncop", "bskewncopp")) {
    if (is.null(marginal)) marginal <- distplyr::dst_unif()
    pdist <- distplyr::get_cdf(marginal)
    if (omit_na) {
    	data <- na.omit(data[c(ycol, x1col, x2col)])
    } else {
    	stop("Sorry, must omit NA's at this time.")
    }
    udat <- as.matrix(data.frame(
        v  = pdist(data[[ycol]]),
        u1 = pdist(data[[x1col]]),
        u2 = pdist(data[[x2col]])
    ))
    xvine <- copsupp::fitrvine_basic(udat, vbls = 2:3)
    u2cond <- copsupp::pcondrvine(udat, xvine, var = 3, condset = 2)
    if (method == "cnqr") {
        fit <- compute_cmc_cnqr(udat = udat, force_ig = force_ig,
                            xvine = xvine, u2cond = u2cond, sc = sc,
                            verbose = verbose, copspace = copspace)
    } else {
        fit <- compute_cmc_mle(udat = udat, force_ig = force_ig, method = method,
                           xvine = xvine, u2cond = u2cond)
    }
    res <- list(xvine  = xvine,
                vine1  = fit$vine1,
                vine2  = fit$vine2,
                data   = data[c(ycol, x1col, x2col)],
                u2cond = u2cond,
                ycol   = ycol,
                x1col  = x1col,
                x2col  = x2col,
                marginal = marginal)
    if (method == "cnqr") res$scorer <- sc
    class(res) <- "cmc"
    res
}
