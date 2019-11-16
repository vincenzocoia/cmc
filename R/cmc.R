#' Fit CMC model via CNQR or MLE
#'
#' Provides a raw workhorse for CNQR (cmc_cnqr_raw) or MLE (cmc_mle_raw).
#' Need to provide it with highly processed input. Intended to be used
#' as internal functions.
#'
#' @param udat matrix of PIT score data, in the order of response (V),
#' first predictor (U1), second predictor (U2), linked by V-U1-U2.
#' @param force_ig TRUE if you want (U1, V) ~ IG copula.
#' @param method Passed to \code{fit_igcop_mle}.
#' @param xvine Copula joining X1 and X2, bundled into a vine object.
#' @param u2cond Vector of U2|U1.
#' @rdname cmc_raw
#' @return Two vines: vine1 is the copula linking (Y, X1), and
#' vine2 is the copula linking (Y, X2)|X1.
cmc_mle_raw <- function(udat, force_ig, method, xvine, u2cond) {
    if (force_ig) {
        cparigcop <- fit_igcop_mle(udat[, 2], udat[, 1], method = method)
        cparmat <- xvine$cparmat # as a template
        cparmat[1, 2] <- list(cparigcop)
        vine1 <- copsupp::rvine(
            G = matrix(c(1, 2,
                         0, 1), byrow = TRUE, ncol = 2),
            copmat = matrix(c("", "igcop"), nrow = 1),
            cparmat = cparmat
        )
    } else {
        vine1 <- copsupp::fitrvine_basic(udat, vbls = 1:2)
    }
    vcond <- copsupp::pcondrvine(udat, vine1, var = 1, condset = 2)
    vine2 <- copsupp::fitrvine_basic(cbind(u2cond, vcond))
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
#' @rdname cmc_raw
cmc_cnqr_raw <- function(udat, force_ig, xvine, u2cond, sc,
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
#' @param marginal Marginal model, assumed to be the same for the response and
#' predictors. If NULL, assumes variables are already PIT scores; otherwise,
#' really is expecting the output from the composite_dist() function.
#' @param method Method of fitting. For MLE (the default), one of
#' "optim" or "nlm", depending
#' on the numerical optimizer you want to use. If "cnqr",
#' fits by CNQR, and you need
#' to specify the sc argument.
#' @param sc Scorer obtained through scorer()
#' @param verbose Only works for CNQR. If TRUE, will output the fitting
#' status of CNQR.
#' @param copspace Only works for CNQR.
#' @param families Vector of copula family names acting as a "pool"
#' to choose from when fitting.
#' @export
cmc <- function(ycol, x1col, x2col, data, method = "optim",
                force_ig = TRUE, marginal = NULL, sc, verbose = FALSE,
                copspace = NULL,
				families = c("indepcop", "bvncop","bvtcop","mtcj","gum",
							 "frk","joe","bb1", "bskewncop", "bskewncopp")) {
    if (is.null(marginal)) marginal <- list(
        pdist = identity,
        qdist = identity
    )
    pdist <- marginal$pdist
    udat <- as.matrix(data.frame(
        v  = pdist(data[[ycol]]),
        u1 = pdist(data[[x1col]]),
        u2 = pdist(data[[x2col]])
    ))
    xvine <- copsupp::fitrvine_basic(udat, vbls = 2:3)
    u2cond <- copsupp::pcondrvine(udat, xvine, var = 3, condset = 2)
    if (method == "cnqr") {
        fit <- cmc_cnqr_raw(udat = udat, force_ig = force_ig,
                            xvine = xvine, u2cond = u2cond, sc = sc,
                            verbose = verbose, copspace = copspace)
    } else {
        fit <- cmc_mle_raw(udat = udat, force_ig = force_ig, method = method,
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
