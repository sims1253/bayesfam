#library(brms)


#pdf
#' Title
#'
#' @param x
#' @param mu
#' @param shape
#' @param shift
#' @param log
#'
#' @return
#' @export
#'
#' @examples
dshifted_inv_gaussian <- function(x, mu, shape, shift, log = FALSE) {
  if(!isLogic_len(log)) {
    stop("The log argument has to be a boolean")
  }
  if(isTRUE(any(x <= 0))) {
    stop("Argument has to be > 0")
  }
  if(isTRUE(any(mu <= 0))) {
    stop("Argument has to be > 0")
  }
  if(isTRUE(any(shape <= 0))) {
    stop("Argument has to be > 0")
  }
  if(!lenEqual(list(x, mu, shape, shift), type_check=is.numeric, scalars_allowed = TRUE)) {
    stop("Either type error, or length missmatch")
  }
  brms::dinv_gaussian(x-shift, mu, shape, log)
}

#rng
#' Title
#'
#' @param n
#' @param mu
#' @param shape
#' @param shift
#'
#' @return
#' @export
#'
#' @examples
rshifted_inv_gaussian <- function(n, mu = 1, shape = 1, shift = 0) {
  if(!isNat_len(n)) {
    stop("n has to be a natural scalar")
  }
  if(isTRUE(any(mu <= 0))) {
    stop("Argument has to be > 0")
  }
  if(isTRUE(any(shape <= 0))) {
    stop("Argument has to be > 0")
  }
  if(!lenEqual(list(mu, shape, shift), type_check=is.numeric, scalars_allowed = TRUE)) {
    stop("Either type error, or length missmatch")
  }
  brms::rinv_gaussian(n, mu, shape) + shift
}

posterior_epred_shifted_inv_gaussian <- function(prep) {
    with(prep$dpars, mu + ndt)
}

posterior_predict_shifted_inv_gaussian <- function(i, prep, ...) {
    n <- prep$ndraws
    mu <- brms::get_dpar(prep, "mu", i = i)
    shape <- brms::get_dpar(prep, "shape", i = i)
    ndt <- brms::get_dpar(prep, "ndt", i = i)
    brms::rshifted_inv_gaussian(n, mu, shape, ndt)
}

log_lik_shifted_inv_gaussian <- function(i, prep) {
    mu <- brms::get_dpar(prep, "mu", i = i)
    shape <- brms::get_dpar(prep, "shape", i = i)
    ndt <- brms::get_dpar(prep, "ndt", i = i)
    y <- prep$data$Y[i]
    dshifted_inv_gaussian(y, mu, shape, ndt, log = TRUE)
}

#' Title
#'
#' @param link
#' @param link_shape
#' @param link_ndt
#'
#' @return
#' @export
#'
#' @examples
shifted_inv_gaussian <- function(link = "1/mu^2", link_shape = "log", link_ndt = "identity"){
  family <- brms::custom_family(
      "shifted_inv_gaussian",
      dpars = c("mu", "shape", "ndt"),
      links = c(link, link_shape, link_ndt),
      lb = c(0, 0, -NA),
      ub = c(NA, NA, 0),
      type = "real",
      log_lik = log_lik_shifted_inv_gaussian,
      posterior_predict = posterior_predict_shifted_inv_gaussian,
      posterior_epred = posterior_epred_shifted_inv_gaussian
    )
  #return inv_gaussian_lpdf(y - ndt| mu, shape);
  family$stanvars <- brms::stanvar(
    scode = "
      real shifted_inv_gaussian_lpdf(real y, real mu, real shape, real ndt) {

          real x = y - ndt;
          return 0.5 * log(shape / (2 * pi())) - 1.5 * log(x)
               - 0.5 * shape * square((x - mu) / (mu * sqrt(x)));
      }

      real shifted_inv_gaussian_rng(real mu, real shape, real ndt) {
          real y = normal_rng(0, 1) ^ 2;
          real x = mu + (mu^2 * y) / (2 * shape) - mu / (2 * shape) *
                    sqrt(4 * mu * shape * y + mu^2 * y^2);
          real z = uniform_rng(0,1);

          if(z <= mu / (mu + x))
            return x + ndt;
          else
            return mu^2 / x + ndt;
      }",
    # code ported from BRMS
    block = "functions"
  )
  return(family)
}
