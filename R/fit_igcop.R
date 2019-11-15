#' Fit IG copula via MLE or CNQR
#' @param method "nlm" and "optim" for MLE (representing the name of the
#' numerical optimizer); or "cnqr" for CNQR.
fit_igcop_mle <- function(u, v, method = "nlm", init=c(2,2), ...) {
	nllh <- function(cpar) {
		if (cpar[1] <= 0) return(Inf)
		if (cpar[2] <= 1.2) return(Inf)
		-sum(logdigcop(u, v, cpar))
	}
	if (method == "nlm") {
		.nlm <- nlm(nllh, init, ...)
		return(.nlm$estimate)
	}
	if (method == "optim") {
		.optim <- optim(init, nllh, ...)
		return(.optim$par)
	}
	stop("method must be either 'nlm' or 'optim'.")
}
