test_that("dugamma works", {

  x <- 0.5
  mu <- stats::runif(1, 0.4, 1)
  sigma <- stats::runif(1, 0.4, 1)

  # Inconsistent arguments
  expect_equal(dugamma(-x, mu, sigma), 0)
  expect_equal(dugamma(3, mu, sigma), 0)
  expect_equal(dugamma(0, mu, sigma), 0)
  expect_equal(dugamma(1, mu, sigma), 0)

  expect_warning(expect_true(is.nan(dugamma(x, -mu, sigma))))
  expect_warning(expect_true(is.nan(dugamma(x, mu, -sigma))))

  # Vectorization
  x <- matrix(stats::runif(3 * 4, 0.4, 1), ncol = 4)
  expect_equal(sum(dugamma(x, mu, sigma) -
                     matrix(c(
                       dugamma(x[1, 1], mu, sigma),
                       dugamma(x[2, 1], mu, sigma),
                       dugamma(x[3, 1], mu, sigma),
                       dugamma(x[1, 2], mu, sigma),
                       dugamma(x[2, 2], mu, sigma),
                       dugamma(x[3, 2], mu, sigma),
                       dugamma(x[1, 3], mu, sigma),
                       dugamma(x[2, 3], mu, sigma),
                       dugamma(x[3, 3], mu, sigma),
                       dugamma(x[1, 4], mu, sigma),
                       dugamma(x[2, 4], mu, sigma),
                       dugamma(x[3, 4], mu, sigma)
                     ), ncol = 4)), 0)

  mu <- matrix(stats::runif(3 * 4, 0.4, 1), ncol = 4)
  expect_equal(sum(dugamma(x, mu, sigma) -
                     matrix(c(
                       dugamma(x[1, 1], mu[1, 1], sigma),
                       dugamma(x[2, 1], mu[2, 1], sigma),
                       dugamma(x[3, 1], mu[3, 1], sigma),
                       dugamma(x[1, 2], mu[1, 2], sigma),
                       dugamma(x[2, 2], mu[2, 2], sigma),
                       dugamma(x[3, 2], mu[3, 2], sigma),
                       dugamma(x[1, 3], mu[1, 3], sigma),
                       dugamma(x[2, 3], mu[2, 3], sigma),
                       dugamma(x[3, 3], mu[3, 3], sigma),
                       dugamma(x[1, 4], mu[1, 4], sigma),
                       dugamma(x[2, 4], mu[2, 4], sigma),
                       dugamma(x[3, 4], mu[3, 4], sigma)
                     ), ncol = 4)), 0)


  sigma <- matrix(stats::runif(3 * 4, 0.4, 1), ncol = 4)
  expect_equal(sum(dugamma(x, mu, sigma) -
                     matrix(c(
                       dugamma(x[1, 1], mu[1, 1], sigma[1, 1]),
                       dugamma(x[2, 1], mu[2, 1], sigma[2, 1]),
                       dugamma(x[3, 1], mu[3, 1], sigma[3, 1]),
                       dugamma(x[1, 2], mu[1, 2], sigma[1, 2]),
                       dugamma(x[2, 2], mu[2, 2], sigma[2, 2]),
                       dugamma(x[3, 2], mu[3, 2], sigma[3, 2]),
                       dugamma(x[1, 3], mu[1, 3], sigma[1, 3]),
                       dugamma(x[2, 3], mu[2, 3], sigma[2, 3]),
                       dugamma(x[3, 3], mu[3, 3], sigma[3, 3]),
                       dugamma(x[1, 4], mu[1, 4], sigma[1, 4]),
                       dugamma(x[2, 4], mu[2, 4], sigma[2, 4]),
                       dugamma(x[3, 4], mu[3, 4], sigma[3, 4])
                     ), ncol = 4)), 0)


  # Inconsistent arguments and vectorization
  x[1, 1] <- -2
  mu[1, 2] <- -2
  mu[2, 3] <- -1
  sigma[2, 1] <- -2
  sigma[3, 4] <- -2

  expect_warning(expect_warning(expect_equal(sum(is.nan(dugamma(x, mu, sigma))), 4)))
  expect_warning(expect_warning(expect_equal(sum(dugamma(x, mu, sigma) == 0, na.rm = TRUE), 1)))

})
