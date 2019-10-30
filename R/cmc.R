#' Fit IG copula via MLE or CNQR
#' @method "nlm" and "optim" for MLE (representing the name of the 
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

#' Fit CMC model via CNQR or MLE
#' 
#' Provides a raw workhorse for CNQR (cmc_cnqr_raw) or MLE (cmc_mle_raw).
#' Need to provide it with highly processed input. Intended to be used
#' as internal functions.
#' 
#' @param udat matrix of PIT score data, in the order of response (V),
#' first predictor (U1), second predictor (U2), linked by V-U1-U2.
#' @param force_ig TRUE if you want (U1, V) ~ IG copula.
cmc_mle_raw <- function(udat, force_ig, method, xvine, u2cond) {
    if (force_ig) {
        cparigcop <- fit_igcop_mle(udat[, 2], udat[, 1], method = method)
        cparmat <- xvine$cparmat # as a template
        cparmat[1, 2] <- list(cparigcop)
        vine1 <- rvine(
            G = matrix(c(1, 2,
                         0, 1), byrow = TRUE, ncol = 2),
            copmat = matrix(c("", "igcop"), nrow = 1),
            cparmat = cparmat
        )
    } else {
        vine1 <- fitrvine_basic(udat, vbls = 1:2)
    }
    vcond <- pcondrvine(udat, vine1, var = 1, condset = 2)
    vine2 <- fitrvine_basic(cbind(u2cond, vcond))
    list(
        vine1 = vine1,
        vine2 = vine2
    )
}

cmc_cnqr_raw <- function(udat, force_ig, method, xvine, u2cond, sc, 
                         verbose, copspace = NULL) {
    if (force_ig) {
        if (is.null(copspace)) {
            copspace <- list("igcop", NULL)
        } else {
            copspace[[1]] <- "igcop"
        }
    }
    fit <- cnqr_reduced(1:3, dat = udat, sc = sc, basevine = xvine, 
                        copspace = copspace, verbose = verbose)
    vine1 <- rvine(
        G       = diag(2),
        copmat  = makevinemat("", fit$copmat[1,3]),
        cparmat = makevinemat(NULL, fit$cparmat[1,3])
    )
    vine2 <- rvine(
        G       = diag(2),
        copmat  = makevinemat("", fit$copmat[2,3]),
        cparmat = makevinemat(NULL, fit$cparmat[2,3])
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
#' @param sc Scorer obtained through cnqr::scorer()
#' @param verbose Only works for CNQR. If TRUE, will output the fitting
#' status of CNQR.
#' @param copspace Only works for CNQR. 
cmc <- function(ycol, x1col, x2col, data, method = "optim",
                force_ig = TRUE, marginal = NULL, sc, verbose = FALSE, 
                copspace = NULL) {
    if (is.null(marginal)) marginal <- list(
        pdist = as_pdist(identity), 
        qdist = as_qdist(identity)
    )
    pdist <- marginal$pdist
    udat <- tibble(
        v  = pdist(data[[ycol]]),
        u1 = pdist(data[[x1col]]),
        u2 = pdist(data[[x2col]])
    ) %>% 
        as.matrix()
    xvine <- fitrvine_basic(udat, vbls = 2:3)
    u2cond <- pcondrvine(udat, xvine, var = 3, condset = 2)
    if (method == "cnqr") {
        fit <- cmc_cnqr_raw(udat = udat, force_ig = force_ig, method = method, 
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
    res
}


#' Predict a quantile from the CMC method
#' 
#' @param tau Single number for quantile level
#' @param type "response" to predict quantiles; "pit" to predict pdist values;
#' "qdist" to predict an entire quantile function; "pdist" to predict an entire
#' pdist (both functions are vectorized).
#' @param resp A response-related value to evaluate the pdist or quantile function
#' at. If type = "response", expecting either a number representing a quantile
#' level, or a column name of the data indicating a column of PIT scores. If
#' type = "pit", expecting either a response value to evaluate the pdist at, or
#' a column name of the data indicating a column of response values. If type = 
#' "pdist" or "qdist", this argument is ignored (since entire functions are returned).
#' @return A vector of function evaluations for each row if type = "response"
#' or "pit"; or, a list of functions for each row if type = "pdist" or "qdist".
# predict.cmc <- function(object, newdata = NULL, resp = 0.9, type = "response") {
#     if (is.null(newdata)) newdata <- object$data
#     
#     udat <- with(object, as.matrix(newdata[, c(ycol, x1col, x2col)]))
#     u1u2cond <- pcondseq(udat, 2:3, object$xvine)
#     QYgX(tau, u1u2cond, 
#          cops = c(object$vine1$copmat[1, 2], 
#                   object$vine2$copmat[1, 2]), 
#          cpars = list(object$vine1$cparmat[1, 2][[1]], 
#                       object$vine2$cparmat[1, 2][[1]]),
#          QY = identity
#     )
# }

#' 'augment' only predicts on the training data.
#' 
#' @param tau Vector of quantile levels to make predictions at
# augment.cmc <- function(object, newdata = NULL, tau = 0.9) {
#     if (is.null(newdata)) newdata <- object$data
#     yhat <- lapply(tau, function(.tau) predict.cmc(object, tau=.tau))
#     tidyr::unnest(dplyr::tibble(df = list(newdata), .tau = tau, yhat = yhat))
# }


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
        u2cond <- pcondrvine(udat, xvine, var = 3, condset = 2)
    }
    u1u2cond <- cbind(udat[, 2], u2cond)
    list(
        newdata  = newdata,
        u1u2cond = u1u2cond
    )
}

#' Add structure to a null model
#' 
#' Converts marginal Unif(0,1) distributions to have modelled distributions.
#' CURRENTLY ONLY HANDLES PIT SCORES AS VALUES.
#' 
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
construct.cmc <- function(object, newdata = NULL, from_col, to_col) {
    ycol  <- object$ycol
    x1col <- object$x1col
    x2col <- object$x2col
    eval <- eval.cmc(object, newdata)
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

#' @param from_col Name of column containing the quantiles
#' to evaluate the pdists at.
decompose.cmc <- function(object, newdata = NULL, from_col, to_col) {
    ycol  <- object$ycol
    x1col <- object$x1col
    x2col <- object$x2col
    eval <- eval.cmc(object, newdata)
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


predict.cmc <- function(object, newdata = NULL, what = "pdist", at, to_col) {
    ycol  <- object$ycol
    x1col <- object$x1col
    x2col <- object$x2col
    eval  <- eval.cmc(object, newdata)
    u1u2cond <- eval$u1u2cond
    if (what == "pdist") {
        funs <- apply(u1u2cond, 1, function(row) {
            this_u1u2cond <- matrix(row, nrow = 1)
            as_pdist(function(v) cnqr::FYgX(
                v, this_u1u2cond, 
                cops = c(object$vine1$copmat[1, 2], 
                         object$vine2$copmat[1, 2]), 
                cpars = list(object$vine1$cparmat[1, 2][[1]], 
                             object$vine2$cparmat[1, 2][[1]]),
                FY = identity
            ))
        })
    } else if (what == "qdist" | what == "rdist") {
        funs <- apply(u1u2cond, 1, function(row) {
            this_u1u2cond <- matrix(row, nrow=1)
            as_qdist(function(p) cnqr::QYgX(
                p, this_u1u2cond, 
                cops = c(object$vine1$copmat[1, 2], 
                         object$vine2$copmat[1, 2]), 
                cpars = list(object$vine1$cparmat[1, 2][[1]], 
                             object$vine2$cparmat[1, 2][[1]]),
                QY = identity
            )[1, ])
        })
        if (what == "rdist") {
            qdists <- funs
            funs <- lapply(qdists, function(qdist) 
                as_rdist(function(n) qdist(runif(n)))
            )
        }
    } else {
        stop("Don't know what to do with '", what, 
             "' as entry of 'what' argument.")
    }
    if (missing(at)) {
        list_output <- funs
    } else {
        list_output <- map(funs, function(fun) {
            res <- tibble(.at  = at,
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
