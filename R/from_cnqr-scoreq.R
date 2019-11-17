#' Score Individual Quantiles
#'
#' Score quantile forecasts against realizations.
#'
#' \code{scoreq} is deprecated. It doesn't allow for the computation of standard
#' error, and if asked to return the matrix of scores, it would incorporate
#' the across-observation weights into the scores, whereas \code{score_eval}
#' does not.
#'
#' @param y Vector of realizations
#' @param yhat Forecast matrix. Each row should be a forecast (of quantiles)
#' with corresponding realization in \code{y}; columns represent the level/index
#' of the quantile (corresponding to \code{tau}).
#' @param sc The scoring rule to use, as in the output of the function
#' \code{\link{scorer}}.
#' @param tau Vector of quantile indices of the quantiles.
#' @param w Vector of weights corresponding to the observations in
#' \code{yhat}. No need to normalize them. See details to see exactly how these
#' weights are used.
#' @param g Increasing univariate vectorized function to transform \code{y}
#' and \code{yhat} with (see computation details in the "details" section).
#' Default is the identity.
#' @param se Logical; should an estimate of the standard error of the mean
#' estimate be returned as well? \code{TRUE} if so (which also overrides
#' the \code{cmb} argument, which is taken to be \code{TRUE}).
#' @param cmb Logical; should the scores be combined (via average)? \code{TRUE}
#' if so, \code{FALSE} to output a score matrix for each observation (rows)
#' and each quantile level (column).
#' @param wtau Function that accepts a vector of quantile indices and returns
#' an equally lengthed vector of weights to multiply the corresponding
#' individual quantile scores by. \code{NULL} for equal weights.
#' @param na_omit Logical; should observations leading to an \code{NA} score
#' (for any \code{tau}) be removed? \code{TRUE} by default. Warning message
#' appears when observations are removed.
#'
#' @details Here's how the score for the \eqn{i}'th observation
#' and the \eqn{k}'th quantile forecast for that observation is computed:
#' \deqn{wtau(\tau_k) (\tau_k - I(y<yhat_{ik}))(g(y) - g(yhat_{ik})),}
#' which is a proper scoring rule as shown in Gneiting and Raftery (2007).
#'
#' To get a score for a particular observation, the average (not the sum)
#' is taken for each row. The scores aren't summed, so that the score doesn't
#' tend to infinity as we include more and more quantiles. Also, the
#' across-quantile weights, determined by the function \deqn{wtau}, are not
#' normalized, so that the individual scores don't tend to 0 as more
#' quantiles are included.
#' @return
#' If \code{se} is \code{TRUE}, returns a named vector of length two of the
#' average score (weighted by argument \code{w}) and standard error of the
#' average. The standard error is estimated by assuming iid scores, and is the
#' standard deviation of the scores times the root sum of squares of the
#' normalized weights \code{w}.
#'
#' Here's what is output if \code{se} is \code{FALSE} (always the case with
#' the deprecated \code{scoreq} function).
#' If \code{cmb} is \code{FALSE}, returns the score matrix (see details for
#' how each score is computed) (rows correspond
#' to observations, and columns correspond to quantile indices \code{tau}).
#' Otherwise, a single numeric score is combined that is the average of the
#' score matrix.
#' @note
#' You could consider having the transformation function \code{g} transform
#' each observation differentially, by forcing it to accept a vector of
#' length equal to your data. This is useful to add seasonal trends, for
#' example.
#' @references
#' \itemize{
#'      \item{Gneiting, T. and Raftery, A. E. (2007). Strictly proper scoring
#'      rules, prediction, and estimation. Journal of the American Statistical
#'      Association, 102(477):359–378.}
#' }
#' @rdname scoring
#' @export
scoreq <- function(y, yhat, tau, w = 1, g = identity, cmb = TRUE,
                   wtau = NULL, na_omit=TRUE) {
    n <- length(y)
    ntau <- length(tau)
    if (is.vector(yhat)) yhat <- matrix(yhat, ncol = 1)
    if (ntau != ncol(yhat))
        stop("quantile indices 'tau' and columns of yhat are not equal in number.")
    if (n != nrow(yhat))
        stop("number of observations in y and yhat do not match.")
    ## Normalize weights
    wtot <- sum(w * rep(1, n))
    wn <- w / wtot * n
    ## Get tau-weights
    if (is.null(wtau)) {
        wtau_comp <- rep(1, ntau)
    } else {
        wtau_comp <- wtau(tau)
    }
    if (length(wtau_comp) == 1) wtau_comp <- rep(wtau_comp, ntau)
    ## Get matrix of scores
    score <- sapply(1:length(tau), function(k)
        wn * wtau_comp[k] * asymloss(tau[k], g(y) - g(yhat[, k])))
    if (na_omit) {
        score2 <- na.omit(score)
        n_rmvd <- nrow(score) - nrow(score2)
        if (n_rmvd > 0) warning(paste("Removed", n_rmvd, "NA observations from score."))
        score <- score2
    }
    if (cmb) return(mean(score)) else return(score)
}


#' @rdname scoring
#' @export
score_eval <- function(y, yhat, sc, w=1, cmb=TRUE, na_omit=TRUE, se=FALSE) {
    n <- length(y)
    tau <- sc$tau
    wtau <- sc$wtau
    g <- sc$g
    ntau <- length(tau)
    if (is.vector(yhat)) yhat <- matrix(yhat, ncol = 1)
    if (ntau != ncol(yhat))
        stop("quantile indices 'tau' and columns of yhat are not equal in number.")
    if (n != nrow(yhat))
        stop("number of observations in y and yhat do not match.")
    ## Get matrix of scores
    scoremat <- sapply(1:length(tau), function(k)
        wtau[k] * asymloss(tau[k], g(y) - g(yhat[, k])))
    if (length(tau) == 1) scoremat <- matrix(scoremat)
    if (na_omit) {
        scoremat2 <- na.omit(scoremat)
        n_rmvd <- nrow(scoremat) - nrow(scoremat2)
        if (n_rmvd > 0) warning(paste("Removed", n_rmvd, "NA observations from score."))
        scoremat <- scoremat2
    }
    ## Compute standard error, if asked:
    scorevec <- apply(scoremat, 1, mean)
    ## Normalized weights (across-observation)
    if (length(w) == 1) {
        wn <- rep(1/n, n)
    } else {
        wtot <- sum(w * rep(1, n))
        wn <- w / wtot
    }
    if (se) {
        stderr <- sd(scorevec) * sqrt(sum(wn^2))
        res <- c(mean(scorevec), stderr)
        names(res) <- c("average", "se")
        return(res)
    }
    avg <- sum(wn*scorevec)
    if (cmb) return(avg) else return(scoremat)
}


