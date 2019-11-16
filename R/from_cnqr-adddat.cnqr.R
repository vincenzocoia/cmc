#' Add data to cnqr object
#'
#' Add more data to a \code{'cnqr'} object.
#'
#' @param obj Object of type \code{'cnqr'} that you wish to add
#' data to. See \code{\link{cnqr}}.
#' @param dat A list of data frames/matrices (of the raw data) to add
#' to \code{obj}. Or, a single data frame/matrix if there's only one.
#' @param QY List of quantile functions of length matching the number
#' of data matrices in \code{dat} (or a single quantile function if there's
#' only one data matrix). \code{NULL} if you'd like to use the quantile
#' function in \code{obj} -- in which case the response data are stationary.
#' @return The original \code{'cnqr'} object \code{obj}, with the
#' \code{$y}, \code{$uind}, \code{$yhat}, \code{$score}, and possibly
#' \code{$QY} (if non-stationary) entries extended to include the new data.
#' @note If \code{dat} is a named list (and not a matrix), then the data added to the
#' \code{'cnqr'} object will also have those names. Otherwise,
#' In the outputted \code{'cnqr'} object, the added data will have names
#' \code{'dat'} with the data number beside it. So for example,
#' if \code{obj} contains two data sets, and you want to add two more unnamed
#' data sets, they'll be named \code{'dat3'} and \code{'dat4'}.
#'
#' This function is useful for cutting down on
#' computation time. Particularly, it helps so that the independent
#' predictors are only needed to be computed once.
#' @rdname adddat
#' @export
adddat.cnqr <- function(obj, dat, QY=NULL) {
    info <- xylink(obj)
    tau <- obj$scorer$tau
    cdf <- obj$cdf
    ycol <- tail(obj$G[1, ], 1)
    ## Make dat a list.
    if (is.vector(dat) & !is.list(dat)) dat <- list(matrix(dat))
    if (is.data.frame(dat) | is.matrix(dat)) dat <- list(as.matrix(dat))
    ## Give the data names if they aren't already named.
    if (is.null(names(dat))) {
        nexist <- length(obj$y)
        names(dat) <- paste0("dat", nexist + 1:length(dat))
    }
    ndat <- length(dat)
    ## Extract response; uniformize data
    y <- lapply(dat, function(dat_) dat_[, ycol])
    dat <- lapply(dat, dat2udat, cdf=cdf)
    ## Get quantile function in list form, and figure out whether or not
    ##  responses are stationary.
    if (is.null(QY)) {
        stationary <- TRUE
        QY <- rep(list(obj$QY), ndat)
    } else {
        stationary <- FALSE
        if (is.function(QY)) {
            if (ndat > 1)
                warning(paste("Only one quantile function was input, and",
                              ndat, "data matrices. Using that quantile function",
                              "for each data matrix."))
            QY <- rep(list(QY), ndat)
        }
    }
    names(QY) <- names(dat)
    ## Get data, independent predictors, and predictions.
    uind <- lapply(dat, pcondseq, ord=info$xord, rv=obj)
    yhat <- mapply(function(uind_, QY_) {
        list(QYgX(tau, uind_, cops=info$cops, cpars=info$cpars, QY=QY_))
    }, uind, QY)
    ## Score
    scores <- mapply(function(y_, yhat_) score_eval(y_, yhat_, sc=obj$scorer),
                     y, yhat)
    ## Append new quantities to the cnqr object:
    obj$y <- c(obj$y, y)
    obj$uind <- c(obj$uind, uind)
    if (!stationary) obj$QY <- c(obj$QY, QY)
    obj$yhat <- c(obj$yhat, yhat)
    obj$score <- c(obj$score, scores)
    obj
}

#' @rdname adddat
#' @export
adddat <- function(obj, dat, QY=NULL) UseMethod("adddat")
