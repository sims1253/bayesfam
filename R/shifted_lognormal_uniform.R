#' Probability density and RNG function for the mixture of shifted lognormal and
#'  uniform distributions.
#'
#' @source Some background, discussion and examples at
#' <http://www.martinmodrak.cz/2021/04/01/using-brms-to-model-reaction-times-contaminated-with-errors/>
#'
#' @details The  mixture of shifted lognormal and uniform can be described as
#' \deqn{y_i =
#'   \begin{cases}
#'   u_i  & \mathrm{if} \quad z_i = 0 \\
#'   s_i + r_i  &  \mathrm{if} \quad z_i = 1
#'   \end{cases}
#'   \\
#'   u_i \sim Uniform(0, \alpha) \\
#'   \log(r_i) \sim Normal(\mu_i, \sigma) \\
#'   P(z_i = 0) = \theta}
#'
#'
#' Where θ corresponds to  `mix`, α to  `max_uniform`
#' and \eqn{s_i} to  `shift`.
#'
#' @param n the number of values to draw from the RNG
#' @param y the observed value
#' @param meanlog the mean of the lognormal component
#' @param sdlog the sd of the lognormal component
#' @param mix the probability of the value comming from the uniform component
#' @param shift the shift the lognormal distribution
#' @param max_uniform the maximum value of the uniform component
#' @name shifted_lognormal_uniform_distribution
NULL

#' @rdname shifted_lognormal_uniform_distribution
#' @export
rshifted_lognormal_uniform <- function(
  n,
  meanlog = 0,
  sdlog = 1,
  mix = 0.1,
  shift = 0,
  max_uniform = 100
) {
  stopifnot(is.numeric(n) & length(n) == 1 & n >= 0)
  n <- as.integer(n)
  stopifnot(all(sdlog > 0))
  stopifnot(all(mix >= 0 & mix <= 1))
  stopifnot(all(shift >= 0))
  stopifnot(all(max_uniform > 0))

  ifelse(
    runif(n) < mix,
    runif(n, 0, max_uniform),
    shift + rlnorm(n, meanlog = meanlog, sdlog = sdlog)
  )
}


#' @rdname shifted_lognormal_uniform_distribution
#' @export
dshifted_lognormal_uniform <- function(
  y,
  meanlog = 0,
  sdlog = 1,
  mix = 0.1,
  shift = 0,
  max_uniform = 100
) {
  stopifnot(all(y > 0))
  stopifnot(all(sdlog > 0))
  stopifnot(all(mix >= 0 & mix <= 1))
  stopifnot(all(shift >= 0))
  stopifnot(all(max_uniform > 0))

  unif_llh = dunif(y, min = 0, max = max_uniform, log = TRUE)
  lognormal_llh = dlnorm(
    y - shift,
    meanlog = meanlog,
    sdlog = sdlog,
    log = TRUE
  ) -
    plnorm(max_uniform - shift, meanlog = meanlog, sdlog = sdlog, log.p = TRUE)

  # Computing logsumexp(log(mix) + unif_llh, log1p(-mix) + lognormal_llh)
  # but vectorized
  llh_matrix <- array(
    NA_real_,
    dim = c(2, max(length(unif_llh), length(lognormal_llh)))
  )
  llh_matrix[1, ] <- log(mix) + unif_llh
  llh_matrix[2, ] <- log1p(-mix) + lognormal_llh
  return(apply(llh_matrix, MARGIN = 2, FUN = logsumexp))
}

posterior_predict_shifted_lognormal_uniform <- function(i, prep, ...) {
  if (
    (!is.null(prep$data$lb) && prep$data$lb[i] > 0) ||
      (!is.null(prep$data$ub) && prep$data$ub[i] < Inf)
  ) {
    stop("Predictions for truncated distributions not supported")
  }

  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  mix <- brms::get_dpar(prep, "mix", i = i)
  shiftprop <- brms::get_dpar(prep, "shiftprop", i = i)

  max_shift <- prep$data$vreal1[i]
  max_uniform <- prep$data$vreal2[i]
  shift = shiftprop * max_shift

  return(
    rshifted_lognormal_uniform(
      prep$ndraws,
      meanlog = mu,
      sdlog = sigma,
      mix = mix,
      shift = shift,
      max_uniform = max_uniform
    )
  )
}


posterior_epred_shifted_lognormal_uniform <- function(prep) {
  if (
    (!is.null(prep$data$lb) && any(prep$data$lb > 0)) ||
      (!is.null(prep$data$ub) && any(prep$data$ub < Inf))
  ) {
    stop("Predictions for truncated distributions not supported")
  }

  mu <- brms::get_dpar(prep, "mu")
  sigma <- brms::get_dpar(prep, "sigma")
  mix <- brms::get_dpar(prep, "mix")
  shiftprop <- brms::get_dpar(prep, "shiftprop")

  max_shift <- prep$data$vreal1
  max_uniform <- prep$data$vreal2
  shift = shiftprop * max_shift

  shifted_lognormal_mean <- shift + exp(mu + sigma^2 / 2)
  uniform_mean <- 0.5 * max_uniform

  return(
    mix * uniform_mean + (1 - mix) * shifted_lognormal_mean
  )
}


log_lik_shifted_lognormal_uniform <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  mix <- brms::get_dpar(prep, "mix", i = i)
  shiftprop <- brms::get_dpar(prep, "shiftprop", i = i)

  max_shift <- prep$data$vreal1[i]
  max_uniform <- prep$data$vreal2[i]
  shift = shiftprop * max_shift

  y <- prep$data$Y[i]
  dshifted_lognormal_uniform(
    y,
    meanlog = mu,
    sdlog = sigma,
    mix = mix,
    shift = shift,
    max_uniform = max_uniform
  )
}


#' A mixture of shifted lognormal and uniform distribution suitable for
#' modelling reaction times.
#'
#' A contaminated response time distribution. The mixture can be described as
#' \deqn{y_i =
#'   \begin{cases}
#'   u_i  & \mathrm{if} \quad z_i = 0 \\
#'   p_i s_i  + r_i  &  \mathrm{if} \quad z_i = 1
#'   \end{cases}
#'   \\
#'   u_i \sim Uniform(0, \alpha) \\
#'   \log(r_i) \sim Normal(\mu_i, \sigma) \\
#'   P(z_i = 0) = \theta \\
#'   0 < p_i < 1
#'   }
#' Here \eqn{\mu, \sigma, \theta} (`mix`) and \eqn{p} (`shiftprop`)
#' are estimated, whereas \eqn{s_i} (`max_shift`)
#' and \eqn{\alpha} (`max_uniform`) are given as data via `vreal()`.
#'
#' @details
#' Note that you cannot build this distribution with the built-in support
#' for mixtures in `brms`,
#' because the uniform component is effectively a zero-parameter distribution
#' which cannot be expressed in  `brms`.
#'

#' @source Idea by
#' Nathaniel Haines (https://twitter.com/Nate__Haines), code by Martin Modrák.
#' Some background, discussion and examples at
#' <http://www.martinmodrak.cz/2021/04/01/using-brms-to-model-reaction-times-contaminated-with-errors/>
#'
#' @param link Link function for the location parameter (default: "identity")
#' @param link_sigma Link function for the scale parameter (default: "log")
#' @param link_mix Link function for the mixture parameter (default: "logit")
#' @param link_shiftprop Link function for the shift proportion parameter (default: "logit")
#'
#' @examples library(brms)
#' set.seed(31546522)
#' # Bounds of the data
#' max_shift <- 0.3
#' shift <- runif(1) * max_shift
#' max_uniform <- 10
#' mix <- 0.1
#'
#' # Generate parameters
#' N <- 100
#' Intercept <- 0.3
#' beta <- 0.5
#' X <- rnorm(N)
#' mu <- rep(Intercept, N) + beta * X
#' sigma <- 0.5
#'
#' rt <- rshifted_lognormal_uniform(N, meanlog = mu, sdlog = sigma, mix = mix,
#'                                  shift = shift, max_uniform = max_uniform)
#'
#' dd <- data.frame(rt = rt, x = X,
#'                  max_shift = max_shift, max_uniform = max_uniform)
#'
#' fam <- shifted_lognormal_uniform()
#' fit_mix <- brm(rt | vreal(max_shift, max_uniform) ~ x, data = dd, family = fam,
#'                stanvars = fam$stanvars,
#'                prior = c(prior(beta(1, 5), class = "mix")))
#' plot(fit_mix)
#' @export
shifted_lognormal_uniform <- function(
  link = "identity",
  link_sigma = "log",
  link_mix = "logit",
  link_shiftprop = "logit"
) {
  fam <- brms::custom_family(
    "shifted_lognormal_uniform",
    dpars = c("mu", "sigma", "mix", "shiftprop"), # Those will be estimated
    links = c(link, link_sigma, link_mix, link_shiftprop),
    type = "real",
    lb = c(NA, 0, 0, 0), # bounds for the parameters
    ub = c(NA, NA, 1, 1),
    vars = c("vreal1[n]", "vreal2[n]"), # Data for max_shift and max_uniform (known)
    posterior_predict = posterior_predict_shifted_lognormal_uniform,
    posterior_epred = posterior_epred_shifted_lognormal_uniform,
    log_lik = log_lik_shifted_lognormal_uniform
  )

  fam$stanvars <- brms::stanvar(
    block = "functions",
    scode = "
  real shifted_lognormal_uniform_lpdf(real y, real mu, real sigma, real mix,
                      real shiftprop, real max_shift, real max_uniform) {
    real shift = shiftprop * max_shift;
    if(y <= shift) {
      // Could only be created by the contamination
      return log(mix) + uniform_lpdf(y | 0, max_uniform);
    } else if(y >= max_uniform) {
      // Could only come from the lognormal
      return log1m(mix) + lognormal_lpdf(y - shift | mu, sigma);
    } else {
      // Actually mixing
      real lognormal_llh = lognormal_lpdf(y - shift | mu, sigma);
      real uniform_llh = uniform_lpdf(y | 0, max_uniform);
      return log_mix(mix, uniform_llh, lognormal_llh);
    }
  }

  real shifted_lognormal_uniform_lcdf(real y, real mu, real sigma, real mix,
                      real shiftprop, real max_shift, real max_uniform) {
    real shift = shiftprop * max_shift;
    if(y <= shift) {
      return log(mix) + uniform_lcdf(y | 0, max_uniform);
    } else if(y >= max_uniform) {
      // The whole uniform part is below, so the mixture part is log(1) = 0
      return log_mix(mix, 0, lognormal_lcdf(y - shift | mu, sigma));
    } else {
      real lognormal_llh = lognormal_lcdf(y - shift | mu, sigma);
      real uniform_llh = uniform_lcdf(y | 0, max_uniform);
      return log_mix(mix, uniform_llh, lognormal_llh);
    }
  }

  real shifted_lognormal_uniform_lccdf(real y, real mu, real sigma, real mix,
                      real shiftprop, real max_shift, real max_uniform) {

    real shift = shiftprop * max_shift;
    if(y <= shift) {
      // The whole lognormal part is above, so the mixture part is log(1) = 0
      return log_mix(mix, uniform_lccdf(y | 0, max_uniform), 0);
    } else if(y >= max_uniform) {
      return log1m(mix) + lognormal_lccdf(y - shift | mu, sigma);
    } else {
      real lognormal_llh = lognormal_lccdf(y - shift | mu, sigma);
      real uniform_llh = uniform_lccdf(y | 0, max_uniform);
      return log_mix(mix, uniform_llh, lognormal_llh);
    }

  }
"
  )
  return(fam)
}
