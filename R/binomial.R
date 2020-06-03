#' Title
#'
#' @param Y
#' @param X
#' @param dimV
#' @param nu
#' @param p0
#' @param w0
#' @param optim_method
#'
#' @return
#' @export
#'
#' @examples
gce_mult <- function(Y, X, dimV, nu, p0, w0, optim_method = "BFGS") {

  dimX <- dim(X)
  N <- dimX[1]
  K <- dimX[2]
  if (missing(dimV)) {
    M <- 5
  } else {
    M <- dimV
  }
  v <- seq(from = -1, to = 1, length.out = M)
  if (missing(nu)) nu <- 0.5
  if (missing(p0)) {
    p0 <- rep(1 / N, N)
    p0_is_uniform <- TRUE
  } else {
    p0_is_uniform <- FALSE
  }
  if (missing(w0)) w0 <- array(1 / M, c(N, M))

  lambda0 <- rep(0, K)
  gce_optim <- optim(lambda0, gce_bin_obj, gce_bin_grad,
                     Y = Y, X = X, v = v, nu = nu, p0 = p0,
                     w0 = w0, N = N, K = K, M = M,
                     method = optim_method, hessian = TRUE)

  lambda <- as.matrix(gce_optim$par)

  temp  <- X %*% lambda

  p    <-  p0 * exp(temp / (1 - nu))
  Phi  <-  sum(p)
  p    <-  p / Phi

  w    <-  w0 * exp(matrix(kronecker(temp, v), nrow = N, byrow = T) / nu)
  Psi  <-  apply(w, 1, sum)
  w    <-  w / Psi

  e    <-  w %*% v

  # Marginal effects
  dp_dx <- array(0, c(N, J, K))
  for (i in 1:N) {
    for (j in 1:J) {
      for (k in 1:K) {
        dp_dx[i, j, k] <- (p[i, j] / (1 - nu)) * (lambda[k, j] - sum(lambda[k, ] * p[i, ]))
      }
    }
  }
  ave_dp_dx <- apply(dp_dx, c(2, 3), mean)

  # Information measures, Golan (1988)
  if (p0_is_uniform == TRUE) {
    Sp <- -sum(p * log(p)) / (N * log(J))
  } else {
    Sp <-  sum(p * log(p)) / sum(p0 * log(p0))
  }

  S_p_i <- rep(0, N)
  for (i in 1:N) {
    if (p0_is_uniform == TRUE) {
      S_p_i[i] <- -sum(p[i, ] * log(p[i, ])) / log(J)
    } else {
      S_p_i[i] <-  sum(p[i, ] * log(p[i, ])) / sum(p0[i, ] * log(p0[i, ]))
    }
  }

  # Error bounds based on the Fano's [weaker] inequality
  p_e_i <- rep(0, N)
  for (i in 1:N) {
    p_e_i[i] <- S_p_i[i] - log(J)
  }

  # Entropy ratio Statistic
  if (p0_is_uniform == TRUE) {
    ER <- 2 * N * log(J) * (1 - Sp)
  } else {
    ER <- 2 * (-sum(p0[i, ] * log(p0[i, ]))) * (1 - Sp)
  }

  # Pseudo R-squared
  R2 <- 1 - Sp

  info_estim_all <- list("lambda" = lambda, "hess" = gce_optim$hessian, "p" = p, "w" = w, "e" = e,
                         "marg_eff" = ave_dp_dx,"Sp" = Sp, "S_p_i" = S_p_i, "p_e_i" = p_e_i,
                         "H_p_w" = - gce_optim$value, "ER" = ER, "Pseudo_R2" = R2,
                         "CM" = CM, "conv" = gce_optim$convergence)

  return(info_estim_all)

}



gce_bin_obj <- function(lambda, y, X, v, nu, p0, w0, N, K, M) {

  temp <-  X %*% lambda

  p    <-  p0 * exp(temp / (1 - nu))
  Phi  <-  sum(p)

  w    <-  w0 * exp(matrix(kronecker(temp, v), nrow = N, byrow = T) / nu)
  Psi  <-  apply(w, 1, sum)

  like <- -lambda %*% (X %*% y) + log(Phi) + sum(log(Psi))

}

gce_bin_grad <- function(lambda, y, X, v, nu, p0, w0, N, K, M) {

  temp <-  X %*% lambda

  p    <-  p0 * exp(temp / (1 - nu))
  Phi  <-  sum(p)
  p    <-  p / Phi

  w    <-  w0 * exp(matrix(kronecker(temp, v), nrow = N, byrow = T) / nu)
  Psi  <-  apply(w, 1, sum)
  w    <-  w / Psi
  e    <-  w %*% v

  grad <- -X %*% y + X %*% p + X %*% e

}

