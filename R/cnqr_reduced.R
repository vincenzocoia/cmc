#' Reduced version of cnqr::cnqr().
#'
#' Only difference is that it doesn't try to compute a final score on the
#' fitted model.
#' @export
cnqr_reduced <- function(edges, dat, sc, basevine, pdist=identity,
                 QY=identity, copspace=NULL, refit=FALSE, verbose=FALSE,
                 families = c("indepcop", "bvncop","bvtcop","mtcj","gum",
                              "frk","joe","bb1", "bskewncop", "bskewncopp")) {
    ## --- Get appropriate basevine ---
    if (missing(basevine)) {
        ## basevine was not specified, which means that the predictors (if present)
        ##  are assumed to be independent.
        basevine <- copsupp::rvine(matrix(edges[-1], nrow=1))
    }
    if (is.vector(dat) & !is.list(dat)) dat <- matrix(dat)
    if (is.data.frame(dat) | is.matrix(dat)) {
        # Only one data set has been input.
        novaldat <- TRUE  # There's no validation data.
        dat <- list(as.matrix(dat))
    } else {
        novaldat <- FALSE
        dat <- lapply(dat, as.matrix)
    }
    ## --- Manipulate data ---
    y <- lapply(dat, function(dat) dat[, edges[1]])
    ytr <- y[[1]]
    if (novaldat) yval <- ytr else yval <- y[[2]]
    if (length(pdist) == 1) pdist <- rep(list(pdist), ncol(dat[[1]]))
    dat <- lapply(dat, cnqr::dat2udat, cdf = pdist)
    ## --- Learn about the Quantile Function ---
    ## Is the quantile function the same for each observation? Store the answer
    ##  in `stationary`. Either way, store the training and validation qdist's
    ##  separately.
    if (is.list(QY)) {  # If QY is a list, then no.
        stationary <- FALSE
        QYtr <- QY[[1]]
        QYval <- QY[[2]]
    } else {
        QYtr <- QY
        QYval <- QY
        if (is.matrix(QY(0.5))) { # If QY outputs a matrix, then no.
            stationary <- FALSE
        } else {  # Otherwise, quantile function is the same.
            stationary <- TRUE
        }
    }
    ## --- Extract useful data quantities ---
    ## (a) Which are predictors (and what order)? What's the response?
    xlab <- edges[-1]
    ylab <- edges[1]
    ## (b) Extract the response; put marginals back in.
    ##  NOTE: Uniform responses are needed for getting starting values with MLE.
    if (novaldat) {  # Only one data set has been input.
        ## (c) Define the two training and validation matrices.
        dattr <- dat[[1]]
        datval <- dattr
        ## (d) Map predictors to independent uniform set.
        if (verbose) cat("Computing independent predictors for training data.\n")
        uindtr <- copsupp::pcondseq(dattr, ord=xlab, rv=basevine)
        uindval <- uindtr
    } else {  # Both a training and validation set were input.
        ## (c) Define the two training and validation matrices.
        dattr <- dat[[1]]
        datval <- dat[[2]]
        ## (d) Map predictors to independent uniform set.
        if (verbose) cat("Computing independent predictors for training data.\n")
        uindtr <- copsupp::pcondseq(dattr, ord=xlab, rv=basevine)
        if (verbose) cat("Computing independent predictors for validation data.\n")
        uindval <- copsupp::pcondseq(datval, ord=xlab, rv=basevine)
    }
    vtr <- dattr[, ylab]
    vval <- datval[, ylab]
    ## --- Extract full model space ---
    res <- copsupp::augment(basevine, a=ylab)
    d <- ncol(res$G)
    p <- length(xlab) # This is not necessarily the number of predictors in the vine (=d-1).
    ## Fill-in copspace if not done already.
    nfam <- length(families)
    if (is.null(copspace)) copspace <- rep(list(NULL), p)
    copspace <- as.list(lapply(copspace, function(fams){
        if (is.null(fams)) families else fams
    }))
    ## --- Fitting and selection procedure ---
    for (i in seq_len(p)) {
        if (verbose) cat(paste("--- Fitting edge", i, "of", p, "---\n"))
        ## --- Fit a copula to edge i ---
        ## Extract the required independent predictors
        this_uindtr <- uindtr[, 1:i, drop=FALSE]
        ## Get conditional PIT score of response and most recent predictor,
        ##  given the predictors that have already been fit.
        ucondtr <- uindtr[, i]
        fittedcops <- res$copmat[seq_len(i-1), d]
        if (length(fittedcops) == 0) {
            ## There are no copula families fit yet.
            vcondtr <- vtr
        } else {
            vcondtr <- FYgX(y=vtr,
                            ucond=this_uindtr[, -i, drop=F],  # rmv last ('active') col.
                            cops=fittedcops,
                            cpars=res$cparmat[seq_len(i-1), d],
                            FY=identity)
        }
        ## Loop through candidate copula families, fitting each one to the running vine.
        res_cand <- lapply(copspace[[i]], function(cop){
            if (verbose) cat(paste0("Fitting copula '", cop, "'.\n"))
            ## Get initial parameter estimates, and select copula reflection/permutation.
            init <- cnqr::cpar_init(ucondtr, vcondtr, cop)
            ## Let's just use the copula family that comes out of cpar_init(),
            ##  which may be different than the requested family (due to a
            ##  restriction of VineCopula's BiCopSelect()).
            cop <- init$cop
            cpar <- list(init$cpar)
            ## Get parameter estimates using CNQR on this copula, with the
            ##  training data.
            cparhat <- cnqr::cnqr_est(res, a=xlab[i], cop=cop, cpar_init=cpar, sc=sc,
                                y=ytr, uind=this_uindtr, QY=QYtr, verbose=verbose)
            if (verbose) cat(paste0("\nParameter: (", paste(cparhat[[1]], collapse=", "), ")\n"))
            ## Augment running vine with this fit. The result is a candidate model.
            copsupp::augment(res, a=xlab[i], cop=cop, cpar=cparhat, col=d)
        })
        ## Select the best candidate model on the validation set.
        if (verbose) cat("Selecting best copula family for this edge.\n")
        this_uindval <- uindval[, 1:i, drop=FALSE]
        res <- cnqr::cnqr_sel(res_cand, sc=sc, y=yval, uind=this_uindval, QY=QYval)
        chosen_cop <- tail(cnqr::xylink(res)$cops, 1)
        if (verbose) cat(paste0("Selected '", chosen_cop, "' copula.\n"))
    }
    ## --- Refit entire column --- (if asked)
    if (refit & p>0) {
        if (verbose) print("Refitting the selected copula families altogether.\n")
        ## Get copula families for each edge, and use their parameters as
        ##  starting values:
        cops <- res$copmat[seq_len(p), d]
        cpars <- res$cparmat[seq_len(p), d]
        ## Re-start `res` as having no links with the predictors.
        res <- copsupp::augment(basevine, a=ylab)
        ## Get parameter estimates using *all* the data:
        if (novaldat) {
            ## There's no validation data. Just use training data.
            cparhat <- cnqr::cnqr_est(res, a=xlab, cop=cops, cpar_init=cpars,
                                sc=sc, y=ytr, uind=uindtr, QY=QYtr, verbose=verbose)
        } else {
            ## There's separate validation data. Combine training and validation
            ##  data to use in the estimation.
            y <- c(ytr, yval)
            uind <- rbind(uindtr, uindval)
            if (stationary) {
                QYall <- QY
            } else {
                QYall <- function(tau) rbind(QYtr(tau), QYval(tau))
            }
            cparhat <- cnqr::cnqr_est(res, a=xlab, cop=cops, cpar_init=cpars,
                                sc=sc, y=y, uind=uind, QY=QYall, verbose=verbose)
        }
        ## Bind parameter estimates to the vine.
        res <- copsupp::augment(res, a=xlab, cop=cops, cpar=cparhat, col=d)
    }
    ## --- Convert vine to `cnqr` object ---
    ## (c) Append items to the vine.
    if (novaldat) {
        res$y    <- list(tr = ytr)
        res$uind <- list(tr = uindtr)
    } else {
        res$y    <- list(tr = ytr, val = yval)
        res$uind <- list(tr = uindtr, val = uindval)
    }
    if (!stationary) {
        if (novaldat) {
            QY <- list(tr = QYtr)
        } else {
            QY <- list(tr = QYtr, val = QYval)
        }
    }
    res$QY     <- QY
    res$pdist  <- pdist
    res$scorer <- sc
    class(res) <- c("cnqr", "rvine")
    if (length(dat) > 2) {
        if (stationary) {
            res <- cnqr::adddat(res, dat[-(1:2)])
        } else {
            res <- cnqr::adddat(res, dat[-(1:2)], QY = QY[-(1:2)])
        }
    }
    return(res)
}
