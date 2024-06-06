# Quantile function --------------------------------------------------------------------------------
test_that("qugamma works", {

  x <- stats::runif(1, 0, 1)
  mu <- stats::runif(1, 0.4, 1)
  sigma <- stats::runif(1, 0.4, 1)

  # Inconsistent arguments
  expect_true(is.nan(qugamma(stats::runif(1, -100, 0), mu, sigma)))
  expect_true(is.nan(qugamma(stats::runif(1, 1, 100), mu, sigma)))
  expect_true(qugamma(0, mu, sigma) == 0)
  expect_true(qugamma(1, mu, sigma) == 1)
  expect_warning(expect_true(is.nan(qugamma(x, -mu, sigma))))
  expect_warning(expect_true(is.nan(qugamma(x, mu, -sigma))))

  # Vectorization
  x <- matrix(stats::runif(3 * 4, 0.4, 1), ncol = 4)
  expect_equal(sum(qugamma(x, mu, sigma) -
                     matrix(c(
                       qugamma(x[1, 1], mu, sigma),
                       qugamma(x[2, 1], mu, sigma),
                       qugamma(x[3, 1], mu, sigma),
                       qugamma(x[1, 2], mu, sigma),
                       qugamma(x[2, 2], mu, sigma),
                       qugamma(x[3, 2], mu, sigma),
                       qugamma(x[1, 3], mu, sigma),
                       qugamma(x[2, 3], mu, sigma),
                       qugamma(x[3, 3], mu, sigma),
                       qugamma(x[1, 4], mu, sigma),
                       qugamma(x[2, 4], mu, sigma),
                       qugamma(x[3, 4], mu, sigma)
                     ), ncol = 4)), 0)

  mu <- matrix(stats::runif(3 * 4, 0.4, 1), ncol = 4)
  expect_equal(sum(qugamma(x, mu, sigma) -
                     matrix(c(
                       qugamma(x[1, 1], mu[1, 1], sigma),
                       qugamma(x[2, 1], mu[2, 1], sigma),
                       qugamma(x[3, 1], mu[3, 1], sigma),
                       qugamma(x[1, 2], mu[1, 2], sigma),
                       qugamma(x[2, 2], mu[2, 2], sigma),
                       qugamma(x[3, 2], mu[3, 2], sigma),
                       qugamma(x[1, 3], mu[1, 3], sigma),
                       qugamma(x[2, 3], mu[2, 3], sigma),
                       qugamma(x[3, 3], mu[3, 3], sigma),
                       qugamma(x[1, 4], mu[1, 4], sigma),
                       qugamma(x[2, 4], mu[2, 4], sigma),
                       qugamma(x[3, 4], mu[3, 4], sigma)
                     ), ncol = 4)), 0)


  sigma <- matrix(stats::runif(3 * 4, 0.4, 1), ncol = 4)
  expect_equal(sum(qugamma(x, mu, sigma) -
                     matrix(c(
                       qugamma(x[1, 1], mu[1, 1], sigma[1, 1]),
                       qugamma(x[2, 1], mu[2, 1], sigma[2, 1]),
                       qugamma(x[3, 1], mu[3, 1], sigma[3, 1]),
                       qugamma(x[1, 2], mu[1, 2], sigma[1, 2]),
                       qugamma(x[2, 2], mu[2, 2], sigma[2, 2]),
                       qugamma(x[3, 2], mu[3, 2], sigma[3, 2]),
                       qugamma(x[1, 3], mu[1, 3], sigma[1, 3]),
                       qugamma(x[2, 3], mu[2, 3], sigma[2, 3]),
                       qugamma(x[3, 3], mu[3, 3], sigma[3, 3]),
                       qugamma(x[1, 4], mu[1, 4], sigma[1, 4]),
                       qugamma(x[2, 4], mu[2, 4], sigma[2, 4]),
                       qugamma(x[3, 4], mu[3, 4], sigma[3, 4])
                     ), ncol = 4)), 0)


  # Inconsistent arguments and vectorization
  x[1, 1] <- -2
  mu[1, 2] <- -2
  mu[2, 3] <- -1
  sigma[2, 1] <- -2
  sigma[3, 4] <- -2

  expect_warning(expect_warning(expect_equal(sum(is.nan(qugamma(x, mu, sigma))), 5)))
})
