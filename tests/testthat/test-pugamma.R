test_that("pugamma works", {

  x <- 0.5
  mu <- stats::runif(1, 0.4, 1)
  sigma <- stats::runif(1, 0.4, 1)

  # Inconsistent and special arguments
  expect_equal(pugamma(-x, mu, sigma), 0)
  expect_equal(pugamma(3, mu, sigma), 1)
  expect_equal(pugamma(0, mu, sigma), 0)
  expect_equal(pugamma(1, mu, sigma), 1)

  expect_warning(expect_true(is.nan(dugamma(x, -mu, sigma))))
  expect_warning(expect_true(is.nan(dugamma(x, mu, -sigma))))

  # Vectorization
  x <- matrix(stats::runif(3 * 4, 0.4, 1), ncol = 4)
  expect_equal(sum(pugamma(x, mu, sigma) -
                     matrix(c(
                       pugamma(x[1, 1], mu, sigma),
                       pugamma(x[2, 1], mu, sigma),
                       pugamma(x[3, 1], mu, sigma),
                       pugamma(x[1, 2], mu, sigma),
                       pugamma(x[2, 2], mu, sigma),
                       pugamma(x[3, 2], mu, sigma),
                       pugamma(x[1, 3], mu, sigma),
                       pugamma(x[2, 3], mu, sigma),
                       pugamma(x[3, 3], mu, sigma),
                       pugamma(x[1, 4], mu, sigma),
                       pugamma(x[2, 4], mu, sigma),
                       pugamma(x[3, 4], mu, sigma)
                     ), ncol = 4)), 0)

  mu <- matrix(stats::runif(3 * 4, 0.4, 1), ncol = 4)
  expect_equal(sum(pugamma(x, mu, sigma) -
                     matrix(c(
                       pugamma(x[1, 1], mu[1, 1], sigma),
                       pugamma(x[2, 1], mu[2, 1], sigma),
                       pugamma(x[3, 1], mu[3, 1], sigma),
                       pugamma(x[1, 2], mu[1, 2], sigma),
                       pugamma(x[2, 2], mu[2, 2], sigma),
                       pugamma(x[3, 2], mu[3, 2], sigma),
                       pugamma(x[1, 3], mu[1, 3], sigma),
                       pugamma(x[2, 3], mu[2, 3], sigma),
                       pugamma(x[3, 3], mu[3, 3], sigma),
                       pugamma(x[1, 4], mu[1, 4], sigma),
                       pugamma(x[2, 4], mu[2, 4], sigma),
                       pugamma(x[3, 4], mu[3, 4], sigma)
                     ), ncol = 4)), 0)


  sigma <- matrix(stats::runif(3 * 4, 0.4, 1), ncol = 4)
  expect_equal(sum(pugamma(x, mu, sigma) -
                     matrix(c(
                       pugamma(x[1, 1], mu[1, 1], sigma[1, 1]),
                       pugamma(x[2, 1], mu[2, 1], sigma[2, 1]),
                       pugamma(x[3, 1], mu[3, 1], sigma[3, 1]),
                       pugamma(x[1, 2], mu[1, 2], sigma[1, 2]),
                       pugamma(x[2, 2], mu[2, 2], sigma[2, 2]),
                       pugamma(x[3, 2], mu[3, 2], sigma[3, 2]),
                       pugamma(x[1, 3], mu[1, 3], sigma[1, 3]),
                       pugamma(x[2, 3], mu[2, 3], sigma[2, 3]),
                       pugamma(x[3, 3], mu[3, 3], sigma[3, 3]),
                       pugamma(x[1, 4], mu[1, 4], sigma[1, 4]),
                       pugamma(x[2, 4], mu[2, 4], sigma[2, 4]),
                       pugamma(x[3, 4], mu[3, 4], sigma[3, 4])
                     ), ncol = 4)), 0)


  # Inconsistent arguments and vectorization
  x[1, 1] <- -2
  mu[1, 2] <- -2
  mu[2, 3] <- -1
  sigma[2, 1] <- -2
  sigma[3, 4] <- -2

  expect_warning(expect_warning(expect_equal(sum(is.nan(pugamma(x, mu, sigma))), 4)))
  expect_warning(expect_warning(expect_equal(sum(pugamma(x, mu, sigma) == 0, na.rm = TRUE), 1)))
})
