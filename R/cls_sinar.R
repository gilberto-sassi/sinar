#' Simulating  SINAR(1,1) process with innovations from a poison distribution.
#'
#' \code{sinar_pois} returns a matrix representing a simulated regular lattice
#' from a SINAR(1,1) process with innovations from a poison distribution.
#'
#' This function simulates a regular lattice from the model
#' \deqn{X_{i,j}= a_{10} X_{i-1,j} + a_{01} X_{i,j-1} + a_{11} X_{i-1, j-1} +
#'  \epsilon_{i,j}}
#' where \eqn{\epsilon_{i,j}} is an iid process with poison distribution. Note
#' the \eqn{a_{10}, a_{01}, a_{11}} must belong to the interval \eqn{[0,1]}.
#'
#' @param n_row Number of rows in the simulated lattice.
#' @param n_col Number of columns in the simulated lattice.
#' @param a10 Coefficient from the element \eqn{X_{i-1, j}}.
#' @param a01 Coefficient from the element \eqn{X_{i, j-1}}.
#' @param a11 Coefficient from the element \eqn{X_{i-1, j-1}}.
#' @param l Mean of the poison distribution used as innovations.
#' @return A integer matrix.
#' @export
#'
#' @examples
#' n_row <- 20
#' n_col <- 50
#' a10 <- 0.2
#' a01 <- 0.2
#' a11 <-  0.5
#' l <- 1
#' sinar_pois(n_row, n_col, a10, a01, a11, l)
sinar_pois <- function(n_row, n_col, a10, a01, a11, l) {
  n_row1  <-  n_row + 1e+2
  n_col1 <- n_col + 1e+2
  e <- matrix(stats::rpois(n_row1 * n_col1, l),
              nrow = n_row1, ncol = n_col1)
  x <-  e

  for (j in 2:n_col1) {
    for (i in 2:n_row1) {
      x[i, j] <- stats::rbinom(1, x[i - 1, j], a10) +
        stats::rbinom(1, x[i, j - 1], a01) +
        stats::rbinom(1, x[i - 1, j - 1], a11) + e[i, j]
    }
  }

  x[(n_row1 - n_row + 1):n_row1, (n_col1 - n_col + 1):n_col1]
}

#' Conditional least square estimates for a SINAR(1,1)  process.
#'
#' \code{cls} computes the conditional least square for a process described
#' by
#' \deqn{X_{i,j}= a_{10} X_{i-1,j} + a_{01} X_{i,j-1} + a_{11} X_{i-1, j-1} +
#'  \epsilon_{i,j}}
#' where \eqn{\epsilon_{i,j}} is an iid process with poison distribution. Note
#' the \eqn{a_{10}, a_{01}, a_{11}} must belong to the interval \eqn{[0,1]}.
#' We obtain estimates for \eqn{a_{10}, a_{01}, a_{11}} and \eqn{\mu_\epsilon}.
#' We do not make any asumption about the distribution of the innovation in the
#' process.
#'
#' @param X A integer matrix where each cell is the observed value in the
#' regular lattice.
#' @return a vector with the estimates of \eqn{a_{10}, a_{01}, a_{11}, \mu}.
#' @export
#'
#' @examples
#'
#' data("nematodes")
#' cls(nematodes)
cls <- function(X) {
  n_row <- nrow(X)
  n_col <- ncol(X)

  mA <- matrix(0, nrow = 4, ncol = 4)
  for (i in 2:n_row) {
    for (j in 2:n_col) {
      mA <- mA +
        rbind(c(X[i - 1, j], X[i, j - 1], X[i - 1, j - 1], 1) * X[i - 1, j],
              c(X[i - 1, j], X[i, j - 1], X[i - 1, j - 1], 1) * X[i, j - 1],
              c(X[i - 1, j], X[i, j - 1], X[i - 1, j - 1], 1) * X[i - 1, j - 1],
              c(X[i - 1, j], X[i, j - 1], X[i - 1, j - 1], 1))
    }
  }


  mb <- vector("double", 4)
  for (i in 2:n_row) {
    for (j in 2:n_col) {
      mb <- mb +
        c(X[i - 1, j] * X[i, j],
          X[i, j - 1] * X[i, j],
          X[i - 1, j - 1] * X[i, j],
          X[i, j])
    }
  }

  est <- solve(mA, mb)
  names(est) <- c("a10", "a01", "a11", "mu")

  est
}

#' Empirical estimate for the matrix V in the Klimko-Nelson.
#'
#' \code{emp_V} is the matrix in the Klimko-Nelson seminal paper. Basically,
#' we know
#' \deqn{\sqrt{n}(\hat{a}_{10} - a_{10}, \hat{a}_{01} - a_{01}, \hat{a}_{11} -
#'  a_{11}, \hat{\mu}_\epsilon - \mu_\epsilon)^\top \sim MNV(0, \Sigma)}
#' where
#' \deqn{\Sigma = V^{-1}W V^{-1}.}
#' For more details, check Klimko and Nelson (1978).
#'
#' @param X A integer matrix where each cell is the observed value in the
#' regular lattice.
#' @return The matrix V estimated empirically.
#' @export
#'
#' @examples
#'
#' data("nematodes")
#' emp_V(nematodes)
emp_V <- function(X) {

  V <- matrix(0, 4, 4)
  for (i in 2:nrow(X)) {
    for (j in 2:ncol(X)) {
      V <- V +
        rbind(c(X[i - 1, j], X[i, j - 1], X[i - 1, j - 1], 1) * X[i - 1, j],
              c(X[i - 1, j], X[i, j - 1], X[i - 1, j - 1], 1) * X[i, j - 1],
              c(X[i - 1, j], X[i, j - 1], X[i - 1, j - 1], 1) * X[i - 1, j - 1],
              c(X[i - 1, j], X[i, j - 1], X[i - 1, j - 1], 1))
    }
  }

  V / ((nrow(X) - 1) * (ncol(X) - 1))
}

#' Compute the value of matrix V using the coefficients.
#'
#' \code{V} is the theoretical matrix from Klimko-Nelson for the SINAR(1,1)
#' model. Basically, we know
#' \deqn{\sqrt{n}(\hat{a}_{10} - a_{10}, \hat{a}_{01} - a_{01}, \hat{a}_{11} -
#'  a_{11}, \hat{\mu}_\epsilon - \mu_\epsilon)^\top \sim MNV(0, \Sigma)}
#' where
#' \deqn{\Sigma = V^{-1}W V^{-1}.}
#' For more details, check Klimko and Nelson (1978).
#'
#' @param a10 is the parameter in the equation \eqn{X[i, j]a_{10}X[i - 1, j] +
#' a_{01}X[i, j - 1] + a_{11}X[i - 1, j - 1] + \epsilon_{i,j}}
#' @param a01 is the parameter in the equation \eqn{X[i, j]a_{10}X[i - 1, j] +
#' a_{01}X[i, j - 1] + a_{11}X[i - 1, j - 1] + \epsilon_{i,j}}
#' @param a11 is the parameter in the equation \eqn{X[i, j]a_{10}X[i - 1, j] +
#' a_{01}X[i, j - 1] + a_{11}X[i - 1, j - 1] + \epsilon_{i,j}}
#' @param mu_e is the mean of the innovations \eqn{\epsilon_{i,j}}
#' @param s2_e is the standar deviation of the innovations \eqn{\epsilon_{i,j}}
#' @return The matrix V estimated empirically.
#' @export
#'
#' @examples
#'
#' n_row <- 20
#' n_col <- 50
#' a10 <- 0.2
#' a01 <- 0.2
#' a11 <-  0.5
#' l <- 1 # mean and variance for poison innovations
#'
#' teo_V(a10, a01, a11, l, sqrt(l))
teo_V <- function(a10, a01, a11, mu_e, s2_e) {

  # mean of X_{i,j}
  mu  <-  mu_e / (1 - a10 - a01 - a11)

  A <- rbind(c(-1, a10, a01, a11),
             c(a10, a11, -1, a01),
             c(a01, -1, a11, a10),
             c(a11, a10, a01, -1))

  b <- c(-mu_e^2 - s2_e - 2 * mu * mu_e * (a10 + a01 + a11),
         -mu * mu_e, -mu * mu_e, -mu * mu_e)

  K <- solve(A, b)
  names(K) <- c("k00", "k01", "k10", "k11")

  V <- rbind(c(K["k00"], K["k11"], K["k01"], mu),
             c(K["k11"], K["k00"], K["k10"], mu),
             c(K["k01"], K["k10"], K["k00"], mu),
             c(mu, mu, mu, 1))
  V
}

#' Empirical estimate for the matrix W in the Klimko-Nelson.
#'
#' \code{emp_W} is the matrix in the Klimko-Nelson seminal paper. Basically,
#' we know
#' \deqn{\sqrt{n}(\hat{a}_{10} - a_{10}, \hat{a}_{01} - a_{01}, \hat{a}_{11} -
#'  a_{11}, \hat{\mu}_\epsilon - \mu_\epsilon)^\top \sim MNV(0, \Sigma)}
#' where
#' \deqn{\Sigma = V^{-1}W V^{-1}.}
#' For more details, check Klimko and Nelson (1978).
#'
#' @param X A integer matrix where each cell is the observed value in the
#' regular lattice.
#' @return The matrix \code{W} estimated empirically.
#' @export
#'
#' @examples
#'
#' data("nematodes")
#' emp_V(nematodes)
emp_W <- function(X) {

  theta <- cls(X) # estimates for a10, a01, a11 and mu

  W <- matrix(0, 4, 4)
  for (i in 2:nrow(X)) {
    for (j in 2:ncol(X)) {
        U <-  X[i, j] - (theta["a10"] * X[i - 1, j] + theta["a01"] *
                           X[i, j - 1] +
                         theta["a11"] * X[i - 1, j - 1] + theta["mu"])
      W <- W +
        U^2 * rbind(c(X[i - 1, j], X[i, j - 1], X[i - 1, j - 1], 1) *
                      X[i - 1, j],
                    c(X[i - 1, j], X[i, j - 1], X[i - 1, j - 1], 1) *
                      X[i, j - 1],
                    c(X[i - 1, j], X[i, j - 1], X[i - 1, j - 1], 1) *
                      X[i - 1, j - 1],
                    c(X[i - 1, j], X[i, j - 1], X[i - 1, j - 1], 1))
    }
  }

  W / ((nrow(X) - 1) * (ncol(X) - 1))
}

#' Empirical estimate for the variance of innovations.
#'
#' \eqn{\sigma^2_\epsilon} is the variance the innovations for the
#' \eqn{SINAR(1,1)} model.
#'
#' @param X A integer matrix where each cell is the observed value in the
#' regular lattice.
#' @return The estimated standard deviation in the \eqn{SINAR(1,1)}.
#' @export
#'
#' @examples
#'
#' data("nematodes")
#' var_sinar(nematodes)
var_sinar <- function(X) {

  theta <- cls(X) # estimates for a10, a01, a11 and mu

  lambda <- (1 + theta["a10"]^2 - theta["a01"]^2 - theta["a11"]^2 -
               sqrt((1 + theta["a10"]^2 - theta["a01"]^2 - theta["a11"]^2)^2 -
                      4 * (theta["a10"] + theta["a01"] * theta["a11"])^2)) /
    (2 * (theta["a10"] + theta["a01"] * theta["a11"]))

  eta <- (theta["a01"] + theta["a11"] * lambda) /
    (1 - theta["a10"] * lambda)

  g00 <- mean((X - mean(X))^2)


  s2 <- g00 *
    (1 - (theta["a10"] + theta["a01"] * theta["a11"]) * lambda -
       (theta["a01"] + theta["a10"] * theta["a11"]) * eta -
       theta["a11"]^2) - mean(X) *
    (theta["a10"] * (1 - theta["a10"]) +
       theta["a01"] * (1 - theta["a01"]) +
       theta["a11"] * (1 - theta["a11"]))

  names(s2) <- NULL

  s2
}

#' Empirical estimate for the Covariance matrix in the Klimko-Nelson.
#'
#' \eqn{\Sigma} is the covariance matrix in the Klimko-Nelson seminal paper.
#' Basically, we know
#' \deqn{\sqrt{n}(\hat{a}_{10} - a_{10}, \hat{a}_{01} - a_{01}, \hat{a}_{11} -
#'  a_{11}, \hat{\mu}_\epsilon - \mu_\epsilon)^\top \sim MNV(0, \Sigma)}
#' where
#' \deqn{\Sigma = V^{-1}W V^{-1}.}
#' For more details, check Klimko and Nelson (1978).
#'
#' @param X A integer matrix where each cell is the observed value in the
#' regular lattice.
#' @return The covariance matrix estimated empirically.
#' @export
#'
#' @examples
#'
#' data("nematodes")
#' emp_cov(nematodes)
emp_cov <- function(X) {
  cov <- MASS::ginv(emp_V(X)) %*% emp_W(X) %*% MASS::ginv(emp_V(X)) /
    ((nrow(X) - 1) * (ncol(X) - 1))
  colnames(cov) <- c("a10", "a01", "a11", "mu")
  rownames(cov) <- c("a10", "a01", "a11", "mu")

  cov
}

#' Variance of standard deviation of epsilon.
#'
#' \eqn{\hat{\sigma}_\epsilon} is the standard deviation of \eqn{SINAR(1,1)}
#' model.
#'
#' @param X A integer matrix where each cell is the observed value in the
#' regular lattice.
#' @return The variance of standard deviation of the estimate of
#' \eqn{\sigma_\epsilon}.
#' @export
#'
#' @examples
#'
#' data("nematodes")
#' var_hat_sigma(nematodes)
var_hat_sigma <- function(X) {

  m <- mean(X)
  g00 <- mean((X - mean(X))^2)
  f <- function(theta) {
    names(theta) <- c("a10", "a01", "a11", "mu")

    lambda <- (1 + theta["a10"]^2 - theta["a01"]^2 - theta["a11"]^2 -
                 sqrt((1 + theta["a10"]^2 - theta["a01"]^2 - theta["a11"]^2)^2 -
                        4 * (theta["a10"] + theta["a01"] * theta["a11"])^2)) /
      (2 * (theta["a10"] + theta["a01"] * theta["a11"]))

    eta <- (theta["a01"] + theta["a11"] * lambda) /
      (1 - theta["a10"] * lambda)

    s2 <- g00 *
      (1 - (theta["a10"] + theta["a01"] * theta["a11"]) * lambda -
         (theta["a01"] + theta["a10"] * theta["a11"]) * eta -
         theta["a11"]^2) - m *
      (theta["a10"] * (1 - theta["a10"]) +
         theta["a01"] * (1 - theta["a01"]) +
         theta["a11"] * (1 - theta["a11"]))

    names(s2) <- NULL

    s2
  }

  g <- matrix(numDeriv::grad(f, cls(X)), ncol = 1)
  t(g) %*% emp_cov(X) %*% g
}
