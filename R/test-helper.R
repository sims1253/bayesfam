# used only for an internal function

#' Check that |a - b| < eps. Works with scalars and vectors on any input.
#' For vector input, one may defined, how many |a - b| >= eps are acceptable
#' with r argument. Given it is a "workhorse" used in almost all test-files,
#' this function is also very well tested.
#'
#' @param a numeric scalar or vector a to be compared
#' @param b numeric scalar or vector b to be compared
#' @param eps numeric scalar or vector, setting the max differences, eps > 0
#' @param r optional numeric scalar (r = 0 in default), relative number of values, that may have an difference > eps.
#' @param note optional parameter used for debugging.
#' @param relative bool argument, if set will take the normale difference with euler metric. Default = FALSE
#' @param debug bool argument, if set will printout the difference, in case of failure. Default = FALSE
#' Used internally, too calculate acceptable amount of deviances. Calculated absolute value will be floored.
#'
#' @md
#' @details For vector/scalar combinations, allowed are (r has to always be a scalar!):
#'   * all scalar
#'   * a vector, b scalar, eps scalar -> compare a values to be close enough to b
#'   * a scalar, b vector, eps scalar -> compare b values to be close enough to a
#'   * a scalar, b scalar, eps vector -> compare the difference a, b too all eps
#'   -> case not intended, but would work
#'   * a vector, b scalar, eps vector -> compare each vector-difference entry against each eps
#'   * all vectors -> each entry of |a - b| is compared to the same entry in eps
#'   -> different vector lengths != 1 dissallowed!expect_brms_family
#'
#' @return success or failure with message
#'
#' @examples print(bayesfam:::expect_eps(1, 1.1, 0.2)) # should pass
#' # print(expect_error(bayesfam:::expect_eps(c(0, 1, 3), c(1, 1, 2), 1e-4, 1 / 3)))
#' # should fail (2/3 were wrong, but only 1/3 was allowed)
expect_eps <- function(a, b, eps, r = 0, relative = FALSE, note = NULL, debug = FALSE) {
  # then check, that r is only a scalar. Also check, that r is in range [0, 1)
  if (isFALSE(isNum_len(r) && r >= 0 && r < 1)) {
    stop("The relative number of tolerated deviances r has to be a scalar in [0, 1).")
  }
  if (!isLogic_len(relative)) {
    stop("The relative argument has to be a single boolean value")
  }
  # then check, that all eps are >= 0
  # (changed to >= given the r might also prove interesting for i.e. integer comparisons)
  if (isTRUE(any(eps < 0))) {
    stop("Tried checking against negative differences.")
  }
  if (relative && isTRUE(any(eps >= 1.0))) {
    stop("In relative mode, the eps should be in [0, 1)")
  }
  # then check, if vectors are of same length or length 1
  if (!lenEqual(list(a, b, eps), scalars_allowed = TRUE, type_check = is.numeric)) {
    stop("Used different length of numeric vectors in test. (Or vectors containing NAs)")
  }
  if (!is.null(note) && !isSingleString(note)) {
    stop("If note is to be used, it has to be a single string argument")
  }

  # For the test, calculate the absolute value, of how many entries in |a - b|
  # may be > eps. This is especially important for RNG tests, which will not always be precise
  # (obviously). Opposed to this in PDF comparisons, this figure is usually 0, as they
  # should always produce comparable results within a few eps machine precision.
  vector_length <- max(length(a), length(b), length(eps))
  tolerated_deviances <- floor(vector_length * r)

  # last check, the actual value difference, either as an absolute or relative error
  if (relative) {
    difference <- normale_difference(a, b)
  } else {
    difference <- abs(a - b)
  }
  eps_comparison_wrong <- (difference > eps)

  # convert the logical vector in a sum of how many entries were wrong
  number_deviances <- sum(eps_comparison_wrong, na.rm = TRUE)

  # at the end, check only the allowed number of entries were wrong
  if (isTRUE(number_deviances <= tolerated_deviances)) {
    testthat::succeed()
  } else {
    # in case of a failure, print how many entries have been incorrect and how many were allowed.
    if (isTRUE(tolerated_deviances > 0)) {
      message <- paste0(
        "In expect_eps: ", toString(number_deviances), " of ",
        toString(vector_length), " were bigger, then eps. Only ",
        toString(tolerated_deviances), " allowed!"
      )
    } else {
      message <- paste0(
        "In expect_eps: ", toString(number_deviances), " of ",
        toString(vector_length), " were bigger, then eps. None were allowed!"
      )
    }
    if (debug) {
      print(paste0("relative: ", relative))
      print(paste0("a: ", a, " b: ", b, " dif: ", difference, " eps: ", eps))
    }
    testthat::fail(paste(message, "\nWith relative:", relative, "the max difference was:", max(difference)))
  }
}

#' Uses euler metric for denominator
#'
#' @param va Numeric scalar or vector of entries
#' @param vb Numeric scalar or vector of entries
#' If both va and vb are no scalars, their lengths have to be equal
#'
#' @return Vector of normalized differences
#'
#' @examples print(bayesfam:::normale_difference(c(1, 1, 1, 1, 1), c(-1, 0, 1, 2, 3)))
normale_difference <- function(va, vb) {
  if (!lenEqual(list(va, vb), scalars_allowed = TRUE, type_check = is.numeric)) {
    stop("In normale_difference function, both vector va and vb have to be numeric and of same len (or scalar)")
  }
  difference <- abs(va - vb)
  denominator <- (va^2 + vb^2)^0.5
  # I like euler as a compromise, between va or vb, given we usually do not know
  # which is the correct one.

  result <- difference / denominator
  result[denominator == 0.0] <- 0.0 # those would be NAs, but are clearly valid 0!
  # I think, this should be the only point, where this formula would fail

  return(result)
}

#' Tests if an RNG can achieve a good enough location often enough
#'
#' @param rng_fun RNG function under test
#' @param metric_mu Metric to be used on RNG data (usually mean or median)
#' @param n Sample size for the rng test.
#' @param mu_list Metric data used as RNG argument and to be compared to (usually mean or median)
#' @param aux_list Auxiliary parameter
#' @param mu_eps Acceptable difference of |mu - metric_mu(rng_fun)
#' @param p_acceptable_failures Acceptable rate of failure, relative value of difference bigger mu_eps
#' @param mu_link Default=identity, optional link-function argument, for example
#' useful in link-normal-distributions
#' @param relative True if the error should be relative to the mu_list, Default = FALSE
#' @param debug bool argument, if set will printout the difference, in case of failure. Default = FALSE
#'
#' @return Nothing actually, just wraps the test
#'
#' @examples eps <- 1e-6
#' mu_list <- seq(from = 1 + eps, to = 20, length.out = 10)
#' phis <- seq(from = 2 + eps, to = 20, length.out = 10)
#' result <- bayesfam:::test_rng(
#'   rng_fun = rbetaprime, metric_mu = mean, n = 10000, mu_list = mu_list,
#'   aux_list = phis, mu_eps = 0.2, p_acceptable_failures = 0.05
#' )
#' print(result)
test_rng <- function(rng_fun,
                     metric_mu,
                     n,
                     mu_list,
                     aux_list,
                     mu_eps,
                     p_acceptable_failures,
                     mu_link = identity,
                     relative = FALSE,
                     debug = TRUE) {
  # check, that all function arguments are actually functions
  # TODO: (Is it possible, to check, if they also take the correct arguments?)
  if (isFALSE(is.function(rng_fun) && is.function(metric_mu) && is.function(mu_link))) {
    stop("RNG-, Metric- or mu_link-function argument was not a function!")
  }
  # check the number of samples to be generated and cecked
  if (!(isInt_len(n) && n >= 1)) {
    stop("n must be an integer, positive scalar!")
  }
  # check, that compare eps is a scalar
  # all used moments should only deviate by eps in most tests)
  if (isFALSE(length(mu_eps) == 1)) {
    stop("mu_eps has to be a scalar in this test-function!")
  }
  # all other arguments are passed and then checked in the called functions!

  # prepare the data, use a vector for ease of use
  # allows re-using the expect_eps.
  # As opposed to using a matrix, which would just complicate implementation and comparison.
  len_mu <- length(mu_list)
  len_aux <- length(aux_list)
  expected_mus <- rep(mu_list, times = len_aux)
  rng_mu_list <- vector(mode = "numeric", length = len_aux * len_mu)

  # calculate rng data
  for (i in seq_along(aux_list)) {
    for (j in seq_along(mu_list)) {
      rng_mu_list[(i - 1) * len_mu + j] <-
        metric_mu(
          rng_fun(
            n,
            mu = mu_link(mu_list[j]), aux_list[i]
          )
        )
    }
  }
  # print(normale_difference(rng_mus, expected_mus))
  # now the data was written, compare it
  expect_eps(
    a = rng_mu_list,
    b = expected_mus,
    eps = mu_eps,
    r = p_acceptable_failures,
    relative = relative,
    debug = debug
  )
}

#' Test if an RNG asymptotically approaches the true location
#'
#' Is currently not in a reliable state useful for testing. Needs more thought.
#'
#' @param rng_fun RNG function under test
#' @param metric_mu Metric to be used on RNG data (usually mean or median)
#' @param n_samples Default=c(10, 10000), sample sizes for rng test, len >= 2. Gets sorted in routine.
#' @param mu_list Metric data used as RNG argument and to be compared to (usually mean or median)
#' @param aux_list Auxiliary parameter value list.
#' @param mu_link Default=identity, optional link-function argument, for example
#' @param allowed_failures Default=0.05 marks, that 5 percent of all test cases are allowed
#' to fail for the whole test to still succeed.
#'
#'
#' @return Success or failure with message
#'
#' @examples eps <- 0.001
#' mu_list <- seq(from = 1 + eps, to = 20, length.out = 10)
#' phis <- seq(from = 2 + eps, to = 20, length.out = 10)
#' result <- bayesfam:::test_rng_asym(
#'   rng_fun = rbetaprime,
#'   metric_mu = mean,
#'   mu_list = mu_list,
#'   aux_list = phis,
#' )
#' print(result)
test_rng_asym <- function(rng_fun,
                          metric_mu,
                          n_samples = c(10, 10000),
                          mu_list,
                          aux_list,
                          mu_link = identity,
                          allowed_failures = 0.05) {
  len_n <- length(n_samples)
  if (len_n < 2 || !isNat_len(n_samples, len = len_n)) {
    stop("n_samples to be a vector of at least two positive integer entries")
  }

  n <- sort(n_samples)

  num_failures <- 0

  # Generate a list of mus per mu, aux combination for growing sample sizes
  for (mu in mu_list) {
    for (aux in aux_list) {
      loop_mu_list <- vector(mode = "numeric", length = len_n)
      for (i in seq_along(n)) {
        loop_mu_list[i] <- mu_link(
          metric_mu(
            rng_fun(
              n[i],
              mu = mu, aux
            )
          )
        )
      }
      # Tests if the resulting distances to the true mu reduce with growing sample size
      if (!identical(
        abs(loop_mu_list - mu),
        sort(abs(loop_mu_list - mu), decreasing = TRUE)
      )) {
        num_failures <- num_failures + 1
      }
    }
  }

  num_tests <- length(mu_list) * length(aux_list)
  allowed_failures_abs <- ceiling(num_tests * allowed_failures)
  if (num_failures <= allowed_failures_abs) {
    testthat::succeed()
  } else {
    testthat::fail(paste(
      "Number of allowed failures in asymp test was violated\n",
      "Allowed were", allowed_failures_abs, "of", num_tests, "to fail\n",
      "but actually", num_failures, "number of tests did fail"
    ))
  }
}



#' Tests if an RNG can recover the true quantiles within a margin of error
#'
#' @param rng_fun RNG function under test
#' @param quantile_fun Quantile function related to the rng under test.
#' @param n Sample size for the rng test.
#' @param mu_list Metric data used as RNG argument and to be compared to
#'                (usually mean or median)
#' @param aux_list Auxiliary parameter value list.
#' @param eps Acceptable difference of |mu - metric_mu(rng_fun)
#' @param quantiles Quantiles to test for recovery.
#' @param p_acceptable_failures Acceptable rate of failure, relative value of
#'                              difference bigger mu_eps
#' @param mu_link Default=identity, optional link-function argument, for example
#' useful in link-normal-distributions
#' @param relative True if the error should be relative to the mu_list
#'
#' @return Nothing actually, just wraps the test
#' @examples eps <- 0.001
#' mu_list <- seq(from = 1 + eps, to = 20, length.out = 10)
#' phi_list <- seq(from = 2 + eps, to = 20, length.out = 10)
#' # if working as expected, this test should not print any errors
#' bayesfam:::test_rng_quantiles(
#'   rng_fun = rbetaprime,
#'   quantile_fun = qbetaprime,
#'   n = 10000,
#'   mu_list = mu_list,
#'   aux_list = phi_list,
#'   eps = 0.1,
#'   quantiles = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99),
#'   p_acceptable_failures = 0.1,
#'   relative = TRUE
#' )
test_rng_quantiles <- function(rng_fun,
                               quantile_fun,
                               n,
                               mu_list,
                               aux_list,
                               eps,
                               quantiles,
                               p_acceptable_failures,
                               mu_link = identity,
                               relative = FALSE) {
  for (mu in mu_list) {
    for (aux in aux_list) {
      sample <- rng_fun(
        n,
        mu = mu_link(mu),
        aux
      )
      true_quantiles <- do.call(quantile_fun, list(quantiles, mu_link(mu), aux))
      expect_eps(
        a = true_quantiles,
        b = quantile(sample, quantiles),
        eps = eps,
        r = p_acceptable_failures,
        relative = relative
      )
    }
  }
}

#' brms family expect recovery. Tries linear baysian model y ~ 1.
#' Checking of the arguments done in construct_brms.
#'
#' @param n_data_sampels How many samples per chain. Positive integer scalar. Default = 1000.
#' @param intercept Intercept for data generating RNG.
#' @param ref_intercept Reference intercept to compare model against. If NULL (default) uses the given intercept.
#' @param aux_par Auxiliary parameter of each distribution.
#' @param rng_link Link function pointer used for data generation. Mainly for transformed normal distributions.
#' @param parameter_link Link function pointer for the latent parameters. Used to transform for comparison with ref_intercept
#' @param family brms family under test.
#' @param rng function pointer of bespoke RNG for the family to be tested.
#' @param aux_name brms string of aux_par argument name. Single string.
#' @param seed Seed argument, so that input data is always the same in each test.
#' brms test does not test RNG and is not guaranteed to fit on all data. Positive Integer scalar, Default = 1235813.
#' Seed is stored before test and restored after it finished. If wants not to use a seed set to NA.
#' @param data_threshold Usually unused. But in rare cases, data too close at the boundary may cause trouble.
#' If so, set a two entry real vector c(lower, upper). If one of them is NA, the data will not be capped for that boundary.
#' Default = Null, will be in R terms "invisible" and will not cap any input data.
#' @param thresh Acceptable threshold for quantiles of recovered arguments.
#' Scalar or 2-entry real vector within (0, 1).
#' Vector is used as is, scalar will be interpreted as c(thresh, 1-thresh).
#' Default = 0.05
#' @param debug Scalar Boolean argument, whether debug info is printed or not. Default = False.
#'
#' @return None
#'
#' @examples result <- bayesfam:::expect_brms_family(
#'   intercept = 5,
#'   aux_par = 2,
#'   ref_intercept = 5,
#'   rng_link = identity,
#'   parameter_link = log,
#'   family = betaprime,
#'   rng = rbetaprime,
#'   aux_name = "phi"
#' )
#' print(result)
expect_brms_family <- function(n_data_sampels = 1000,
                               intercept,
                               ref_intercept = NULL,
                               aux_par,
                               rng_link,
                               parameter_link,
                               family,
                               rng,
                               aux_name,
                               seed = 1235813,
                               data_threshold = NULL,
                               thresh = 0.05,
                               debug = FALSE) {
  if (!isSingleString(aux_name)) {
    stop("The aux_par name argument has to be a single string")
  }
  if (is.null(ref_intercept)) {
    ref_intercept <- intercept
  }
  posterior_fit <- construct_brms(n_data_sampels,
    intercept,
    aux_par,
    rng_link = rng_link,
    family,
    rng,
    seed = seed,
    data_threshold = data_threshold
  )

  intercept_recovered <- test_brms_quantile(
    posterior_fit, "b_Intercept", parameter_link(ref_intercept), thresh, debug
  )
  aux_par_recovered <- test_brms_quantile(
    posterior_fit, aux_name, aux_par, thresh, debug
  )
  success <- intercept_recovered & aux_par_recovered

  if (debug & !success) {
    print("Data were not recovered correctly! Print plot.")
    fam_name <- posterior_fit$family$name
    debug <- paste0(
      fam_name,
      " expect_brms_family failed with inputs intercept = ",
      intercept,
      " and aux_par = ",
      aux_par
    )
    plot(posterior_fit, main = debug)
  }

  if (success) {
    testthat::succeed()
  } else {
    testthat::fail("One or more variables have not been recovered correctly! You may set debug true and check the plot.")
  }
}

#' Construct brms family for simple linear y ~ 1 model.
#'
#' @param n_data_sampels How many samples per chain. Positive integer scalar.
#' @param intercept Intercept data argument, real scalar.
#' @param aux_par aux_par argument of each distribution.
#' @param rng_link Link function pointer used data. For positive bounded uses exp as example.
#' @param family brms family under test.
#' @param rng function pointer of bespoke RNG for the family to be tested.
#' @param seed Seed argument, so that input data is always the same in each test.
#' brms test does not test RNG and is not guaranteed to fit on all data.
#' Positive Integer scalar, Default = NA will do nothing. Seed is stored before and restored after.
#' @param data_threshold Usually unused. But in rare cases, data too close at the boundary may cause trouble.
#' If so, set a two entry real vector c(lower, upper). If one of them is NA, the data will not be capped for that boundary.
#' Default = Null, will be in R terms "invisible" and will not cap any input data.
#'
#' @return brms model for the specified family.
#'
#' @examples posterior_fit <- bayesfam:::construct_brms(
#'   n_data_sampels = 1000,
#'   intercept = 5.0,
#'   aux_par = 2.0,
#'   rng_link = identity,
#'   family = betaprime,
#'   rng = rbetaprime
#' )
#' plot(posterior_fit)
#' # beta_prime uses log-link for Intercept
construct_brms <- function(n_data_sampels,
                           intercept,
                           aux_par,
                           rng_link,
                           family,
                           rng,
                           seed = NULL,
                           data_threshold = NULL) {
  if (!(is.function(family) && is.function(rng) && is.function(rng_link))) {
    stop("family, rng or rng_link argument were not a function!")
  }
  if (!(isNat_len(n_data_sampels))) {
    stop("n_data_sampels has to be a positive integer scalar")
  }
  if (!isNum_len(intercept)) {
    stop("intercept argument has to be a real scalar")
  }
  if (!isNum_len(aux_par)) {
    stop("aux_par argument has to be a real scalar")
  }
  if (!(isNum_len(seed) || is.null(seed))) {
    stop("seed argument if used has to be a real scalar. Else it is let default as NULL,
         which will not change the current RNG seed")
  }


  if (!is.null(seed)) {
    old_seed <- .Random.seed
    set.seed(seed)
  }

  y_data <- rng(n_data_sampels, rng_link(intercept), aux_par)
  if (!is.null(data_threshold)) {
    y_data <- limit_data(y_data, data_threshold)
  }

  if (!is.null(seed)) {
    set.seed(old_seed)
  }


  data <- list(y = y_data)

  posterior_fit <- brms::brm(
    y ~ 1,
    data = data,
    family = family(),
    stanvars = family()$stanvars,
    chains = 2,
    cores = 2,
    silent = 2,
    refresh = 0,
    init = 0.1
  )

  return(posterior_fit)
}

#' Check, that data of the posterior is close enough to the reference data.
#'
#' @param posterior_data Data fitted and drawn, a brms object, which is a list in R terms
#' @param arg_name Name of the argument variable to check, as single string
#' @param reference Reference value to check against, single real scalar
#' @param thresh real scalar or 2-length vector of quantile bounds.
#' For scalar constructs bound as [thresh, 1-thresh]
#' thresh has to be inside the Unit-Interval.
#' @param debug True for verbose output of test results.
#'
#' @return Single boolean success, fail or error
#'
#' @examples fit <- bayesfam:::construct_brms(
#'   n_data_sampels = 1000,
#'   intercept = 5.0,
#'   aux_par = 2.0,
#'   rng_link = identity,
#'   family = betaprime,
#'   rng = rbetaprime
#' )
#' result <- bayesfam:::test_brms_quantile(
#'   posterior_data = fit, arg_name = "phi", 2.0, 0.025
#' )
#' plot(fit)
#' # beta_prime uses log-link for Intercept
test_brms_quantile <- function(posterior_data, arg_name, reference, thresh, debug = FALSE) {
  if (!is.list(posterior_data)) {
    stop("The posterior_data frame has to be brms data, which itself is in R of type list")
  }
  if (!isSingleString(arg_name)) {
    stop("The variable arg_name argument has to be a single string")
  }
  if (!isNum_len(reference)) {
    stop("The reference data has to be a single real scalar")
  }
  if (isTRUE(any(thresh <= 0 | thresh >= 1))) {
    stop("The threshold has to be in the Unit-Interval")
  }
  if (!isLogic_len(debug)) {
    stop("The argument debug has to be a single boolean")
  }
  if (isNum_len(thresh)) {
    if (thresh > 0.5) {
      stop("Any threshold scalar bigger 0.5 would create a bigger lower bound than upper bound!")
    }
    bounds <- c(thresh, 1 - thresh)
  } else if (isNum_len(thresh, 2)) {
    if (isFALSE(thresh[1] <= thresh[2])) {
      stop("If a 2 entry vector is used for the bounds, the first entry is the lower bound")
    }
    bounds <- thresh
  } else {
    stop("The quantile-thresholds can only be a scalar, or a 2 length vector")
  }

  calculated <- tryCatch(
    {
      posterior::extract_variable_matrix(posterior_data, variable = arg_name)
    },
    error = function(e) {
      # if the extract fails, most probably cause is a wrong variable arg_name string
      # instead of giving an unreadable error, throw clear warning and return false
      warning(paste0("In test_brms_quantile, the variable extraction of ", arg_name, " failed.
                   Most probable cause, the variable-string was written wrong or did not exist somehow.
                   Return FALSE in this case."))
      # warning(paste0("Original error was: ", e))
      return(FALSE)
    }
  )

  quantiles <- unname(quantile(calculated, probs = bounds))
  if (debug) {
    print(paste("At: ", bounds, " the quantile is: ", quantiles))
    print(paste("The supplied reference value was: ", reference))
  }

  # if(quantiles[1] < reference && reference < quantiles[2])
  #   succeed()
  # else
  #   fail("The reference value was not within the quantiles of the given posterior data!")
  value <- quantiles[1] < reference && reference < quantiles[2]
  return(isTRUE(value))
}

#' Check if a vector is numeric, has no na entries and has length len
#'
#' @param x Numeric vector to be checked
#' @param len Length of vector, default argument is 1
#'
#' @return Boolean, whether x was numeric and of correct size
#'
#' @examples bayesfam:::isNum_len(c(1.1, 2.2), 2) # should be TRUE
#' bayesfam:::isNum_len(0.2) # should be TRUE
#'
isNum_len <- function(x, len = 1) {
  value <- all(!is.na(x)) && all(is.numeric(x)) && length(x) == len
  return(isTRUE(value))
}

#' Integer vector check
#'
#' @param int Integer vector to be checked
#' @param len Length of vector, default argument is 1
#'
#' @return Boolean, whether int was Integer and of correct size
#' @export
#'
#' @examples bayesfam:::isInt_len(c(1, 2), 2) # should be TRUE
#' bayesfam:::isInt_len(1, 2) # should be FALSE, wrong length
isInt_len <- function(int, len = 1) {
  all_numeric <- all(!is.na(int)) && all(is.numeric(int))
  if (isTRUE(all_numeric)) {
    # only check for integer, if the type is numeric!
    value <- all(int %% 1 == 0) && length(int) == len
    return(isTRUE(value))
  } else {
    return(FALSE)
  }
  # One might also do this all in a single AND beginning with is.numeric.
  # In a single test, this worked fine, given if !is.numeric, the other boolean,
  # checks were not done (because FALSE & x <=> FALSE)
  # this did prevent errors (from "r" %% 1 == 0) at the least

  # But I would prefer only checking numerics,
  # given logical AND may not necessarily follow this behavior!
}

#' Natural number vector check n >= 0.
#'
#' @param int Integer vector to be checked
#' @param len Length of vector, default argument is 1
#'
#' @return Boolean, whether int was Integer >= 0 and of correct size
#' @export
#'
#' @examples bayesfam:::isNat_len(c(1, 2), 2) # should be TRUE
#' bayesfam:::isNat_len(-1) # should be FALSE
isNat_len <- function(int, len = 1) {
  return(isInt_len(int, len) && all(int >= 0))
  # May require Unit-Tests (though as a wrapper to isInt_len, may makes it almost redundant)
}

#' Boolean vector check
#'
#' @param logic Logic vector to be checked
#' @param len Length of vector, default argument is 1
#'
#' @return Boolean, whether logic was Boolean and of correct size
#' @export
#'
#' @examples bayesfam:::isLogic_len(c(TRUE, FALSE), 2) # should be TRUE
#' bayesfam:::isLogic_len(0, len = 1) # should be FALSE, 0 and 1 are numeric
isLogic_len <- function(logic, len = 1) {
  if (any(is.function(logic))) {
    # other comparable functions threw warnings for function-ptr.
    # This did not, so I added it in manually.
    warning("Function type given instead of boolean in isLogic_len")
    return(FALSE)
  }
  value <- all(is.logical(logic)) && length(logic) == len
  return(isTRUE(value))
}

#' Check, if the input is a single string
#'
#' @param input String argument
#'
#' @return Is a string and only one string
#' @export
#'
#' @examples bayesfam:::isSingleString("abc") # should be TRUE
#' bayesfam:::isSingleString(c("abc", "def")) # should be FALSE, not a single string
isSingleString <- function(input) {
  value <- all(!is.na(input)) && is.character(input) && length(input) == 1
  return(isTRUE(value))
}

#' Length and data check function
#'
#' @param list_of_vectors List of vectors of any type you want
#' @param scalars_allowed In many applications, mixing scalars with vectors may be fine,
#' so this argument controls, if scalars (length == 1) would be fine to. Default = FALSE
#' @param type_check Function pointer argument of type to be checked. All vectors have to
#' be of this type, if defined. Default = NULL
#' @param na_allowed True if test shall pass with NAs in the data.
#'
#' @return boolean, given the arguments above
#'
#' @examples va <- c(1, 2, 3)
#' vb <- c(4, 5, 6)
#' vc <- c(7, 8)
#' bayesfam:::lenEqual(list(va, vb)) # both got 3 entries
#' bayesfam:::lenEqual(list(va, vb, vc)) # not all vectors have the same number of entries
lenEqual <- function(list_of_vectors, scalars_allowed = FALSE, type_check = NULL, na_allowed = FALSE) {
  if (!isLogic_len(scalars_allowed)) {
    stop("scalars_allowed has to be a single boolean value")
  }
  if (!isLogic_len(na_allowed)) {
    stop("na_allowed has to be a single boolean value")
  }
  if (!is.null(type_check) && !is.function(type_check)) {
    stop("If type_check is to be used, it has to be a function pointer, to a type checking function
         (for example is.numeric)")
  }

  maxLen <- 0
  for (vector in list_of_vectors) {
    currentLen <- length(vector)
    if (currentLen > maxLen) {
      maxLen <- length(vector)
    }
    # gets the vector of the greatest length

    if (is.function(type_check) && !type_check(vector)) {
      warning("At least one vector was not of the specified type! Return FALSE immediatly.")
      return(FALSE)
    }
    if (!na_allowed && any(is.na(vector))) {
      # if NAs are disallowed and the input contains any NAs, the function returns FALSE immediatly.
      warning("NAs disallowed, but at least one entry in the vectors was a NA! Return FALSE immediatly.")
      return(FALSE)
    }
  }

  for (vector in list_of_vectors) {
    currentLen <- length(vector)
    scalar_correct <- scalars_allowed && currentLen == 1
    if (!scalar_correct && currentLen != maxLen) {
      # if not scalar, nor the correct max len, then return false immediatly.
      return(FALSE)
    }
  }

  # if no issues occured so far, the result must be of the correct length.
  return(TRUE)
}

#' Data limit function
#'
#' @param data Data to be limited
#' @param limits Limits to be used. Vector with 2 real entries, limits[1] <= limits[2]
#' If the lower bound does not have to be restricted, set it to NA and vice versa.
#' Sets data outside those bounds to those bounds.
#'
#' @return data limited by the limits
#' @export
#'
#' @examples input <- c(1, 2, 3, 4)
#' print(bayesfam:::limit_data(input, c(2, 3))) # lower and upper bounds
#' print(bayesfam:::limit_data(input, c(2, NA))) # only lower bound
limit_data <- function(data, limits) {
  # check that the limit is usable
  if (length(limits) != 2) {
    stop("If the limits is to be used, it has to be of size 2.")
  }
  # isNum_len also checks, that data does not contain any NAs
  if (!isNum_len(data, len = length(data))) {
    stop("Some data was not numeric, or was NA")
  }

  # if both bounds are used, check the order of them
  if (isNum_len(limits, 2)) {
    if (limits[1] > limits[2]) {
      stop("In limit_data, the first limit is the lower limit, so it has to be
           smaller than the second limit.")
    }
  }

  if (isNum_len(limits[1])) {
    data[data < limits[1]] <- limits[1]
  }
  if (isNum_len(limits[2])) {
    data[data > limits[2]] <- limits[2]
  }

  # now return the data
  return(data)
}
