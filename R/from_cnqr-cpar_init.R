#' Get Initial Values for Copula
#'
#' Obtains initial parameter estimates for a bivariate copula, using MLE.
#' This includes the copula reflection. Uses \code{\link{fitbicop_lh}}.
#' If the IG, IGL, or skew normal copula families are requested, then MLE
#' is not used (since the copulas are not recognized by the VineCopula pkg).
#' Pre-chosen values are output instead, and no reflections are considered.
#'
#' @param u,v Vectors of uniform values.
#' @param cop Character of a (single) copula family to consider fitting.
#' @export
cpar_init <- function(u, v, cop) {
    ## Use the VineCopula package's `BiCopSelect()` to get initial parameter
    ##  values, and the copula reflection. Do this through copsupp's
    ##  `fitbicop_lh()` function, which handles the odd output of
    ##  BiCopSelect (which, for example, makes the parameter negative
    ##  when there's a reflection). Since it can't handle the
    ##  bivariate skew normal copula, and the IG and IGL copulas, just
    ##  deal with those separately.
    if (cop == "bskewncop" | cop == "bskewncopp")
        return(list(cop=cop, cpar=c(0.1,0.1,0.1)))
    if (cop == "igcop")
        return(list(cop=cop, cpar=c(2.1, 3.1)))
    if (cop == "iglcop")
        return(list(cop=cop, cpar=3.1))
    if (cop %in% c("bb1rsk", "bb1sk", "bb1vsk", "bb1skp")) {
        res <- copsupp::fitbicop_lh(u, v, families="bb1")
        return(list(cop=cop, cpar=c(res$cpar, 0.8)))
    }
    res <- copsupp::fitbicop_lh(u, v, families=cop)
    list(cop=res$cop, cpar=res$cpar)
}
