#' Quick Normal Scores Plot
#'
#' Produces a quick base R normal scores plot.
#' @param .data Data frame containing the data.
#' @param x_col,y_col Names of the columns containing the
#' X and Y variables, respectively.
#' @export
nscoreplot <- function(.data, x_col, y_col) {
	x <- .data[[x_col]]
	y <- .data[[y_col]]
	u <- stats::ecdf(x)(x)
	v <- stats::ecdf(y)(y)
	graphics::plot(stats::qnorm(u), stats::qnorm(v))
}
