set.seed(100)
n <- 100
vine <- copsupp::rvine(
	copsupp::AtoG(diag(1:2)),
	copmat  = copsupp::makevinemat("bvtcop", zerocol = TRUE),
	cparmat = copsupp::makevinemat(list(c(0.5, 4)), zerocol = TRUE)
)
udat <- copsupp::rrvine(n, vine)

test_that("cnqr_reduced provides sensible output", {
	fit_gum <- cnqr_reduced(1:2, udat, sc = scorer(1:9/10),
							 families = "gum")
	fit_nor <- cnqr_reduced(1:2, udat, sc = scorer(1:9/10),
							 families = "bvncop")
	expect_true(fit_gum$copmat[1, 2] %in% c("gum", "gumr"))
	expect_true(fit_nor$copmat[1, 2] %in% "bvncop")
	expect_gt(fit_nor$cparmat[1, 2][[1]], 0)
	expect_identical(dim(fit_gum$G), c(2L, 2L))
})
