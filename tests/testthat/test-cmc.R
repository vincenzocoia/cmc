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
						sc = scorer(1:9/10), verbose = FALSE,
						families = families)
	# expect_identical(fit$vine1$G, matrix(c(1, 0, 2, 1), ncol = 2))
	expect_true(fit$vine1$copmat[1,2] %in% c(families, paste0(families, "r")))
	expect_true(fit$vine2$copmat[1,2] %in% c(families, paste0(families, "r")))
	fit_norm <- cmc_cnqr_raw(dat, force_ig = FALSE, xvine = vine_x, u2cond = u2cond,
							 sc = scorer(1:9/10), verbose = FALSE,
							 families = "bvncop")
	expect_gt(fit_norm$vine1$cparmat[1,2][[1]], 0)
	expect_gt(fit_norm$vine2$cparmat[1,2][[1]], 0)
})

test_that("cmc() matches the raw versions", {
	families <- c("bvncop", "bvtcop", "mtcj", "gum", "frk", "joe", "bb1", "bb7",
				  "bb8")
	dat <- as.data.frame(dat)
	dat <- stats::setNames(dat, c("v", "u1", "u2"))
	dat_shuff <- dat[, c("u1", "u2", "v")] # checking robustness to order
	fit <- cmc("v", "u1", "u2", data = dat_shuff, force_ig = FALSE)
	vine_x_est <- copsupp::fitrvine_basic(dat_shuff, vbls = 1:2)
	pcondxcop <- getFromNamespace(paste0("pcond", vine_x_est$copmat[1, 2]),
								  ns = "CopulaModel")
	xcop_cpar <- vine_x_est$cparmat[1, 2][[1]]
	u2cond <- pcondxcop(dat_shuff$u2, dat_shuff$u1, xcop_cpar)
	fit_mle_raw <- cmc_mle_raw(dat, xvine = vine_x_est, force_ig = FALSE, u2cond = u2cond)
	expect_identical(fit_mle_raw$vine1$copmat, fit$vine1$copmat)
	expect_identical(fit_mle_raw$vine1$cparmat, fit$vine1$cparmat)
	expect_identical(fit_mle_raw$vine2$copmat, fit$vine2$copmat)
	expect_identical(fit_mle_raw$vine2$cparmat, fit$vine2$cparmat)
})

test_that("predictions are sensible", {
	set.seed(100)
	n <- 100
	x <- runif(n)
	dat <- data.frame(v  = x,
					  u1 = x,
					  u2 = x)
	fit <- cmc("v", "u1", "u2", data = dat, method = "cnqr", sc = scorer(1:9/10))
})
