set.seed(1)
copmat <- copsupp::makevinemat(
	"gum",
	c("bvtcop", "frk"),
	zerocol = TRUE
)
cparmat <- copsupp::makevinemat(
	3.1,
	list(c(0.5, 4), 2.3),
	zerocol = TRUE
)
G <- copsupp::AtoG(CopulaModel::Dvinearray(3))
vine_full <- copsupp::rvine(G, copmat, cparmat)
vine_x <- copsupp::subset(vine_full, 2:3)
dat <- copsupp::rrvine(100, vine_full)
u2cond <- copsupp::pcondrvine(dat, vine_full, 3, 2)


# n <- 100
# u1 <- runif(n)
# x1 <- qnorm(u1)
# tau2 <- runif(n)
# u2 <- CopulaModel::qcondbvtcop(tau2, u1, c(0.5, 4))
# x2 <- qnorm(u2)
# tauy <- runif(n)
# v <- diag(QYgX(tauy, data.frame(u1, tau2),
# 		  cops = c("bvncop", "bvncop"),
# 		  cpars = list(0.5, 0.5)))
# y <- qnorm(v)
# dat <- data.frame(x1 = x1, x2 = x2, y = y,
# 				  u1 = u1, u2 = u2, v = v)

# foodat <- dat
# foodat[, 2] <- 1 - foodat[, 2]
test_that("cmc_mle_raw gives something sensible", {
	families <- c("bvncop", "bvtcop", "mtcj", "gum", "frk", "joe", "bb1", "bb7",
				  "bb8")
	fit <- cmc_mle_raw(dat, force_ig = FALSE, xvine = vine_x, u2cond = u2cond)
	expect_identical(fit$vine1$G, matrix(c(1, 0, 2, 1), ncol = 2))
	expect_true(fit$vine1$copmat[1,2] %in% c(families, paste0(families, "r")))
	expect_true(fit$vine2$copmat[1,2] %in% c(families, paste0(families, "r")))
})

test_that("cmc_cnqr_raw gives something sensible", {
	families <- c("bvncop", "bvtcop", "gum")
	fit <- cmc_cnqr_raw(dat, force_ig = FALSE, xvine = vine_x, u2cond = u2cond,
						sc = cnqr::scorer(1:9/10), verbose = TRUE,
						families = families)
	# expect_identical(fit$vine1$G, matrix(c(1, 0, 2, 1), ncol = 2))
	expect_true(fit$vine1$copmat[1,2] %in% c(families, paste0(families, "r")))
	expect_true(fit$vine2$copmat[1,2] %in% c(families, paste0(families, "r")))
	fit_norm <- cmc_cnqr_raw(dat, force_ig = FALSE, xvine = vine_x, u2cond = u2cond,
							 sc = cnqr::scorer(1:9/10), verbose = TRUE,
							 families = "bvncop")
	expect_gt(fit_norm$vine1$cparmat[1,2][[1]], 0)
	expect_gt(fit_norm$vine2$cparmat[1,2][[1]], 0)
})
