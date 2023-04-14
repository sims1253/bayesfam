rasym_laplace_custom <- function(n, mu = 0, sigma = 1, beta = 1) {
  if(!isNat_len(n)) {
    stop("n has to be a scalar natural number")
  }
  if(isTRUE(any(sigma <= 0))) {
    stop("sigma has to be bigger than 0")
  }
  if(isTRUE(any(beta <= 0))) {
    stop("beta has to be bigger than 0")
  }
  if(!lenEqual(list(mu, sigma, beta), scalars_allowed = TRUE,
               type_check = is.numeric, na_allowed = FALSE)) {
    stop("Type error, or unequal length vectors with n > 1")
  }
  return(brms::rexgaussian(n,  mu - beta, sigma, beta))
}
