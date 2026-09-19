# Coverage and consistency of the cross-layer family inventory (issue #28).
#
# These tests fail loudly when a family is missing from the inventory, when a
# custom family does not register a callback, or when a registered callback has
# no declared layer-3 coverage ("ok", "unsupported" with a documented reason,
# or "skip" referencing a wave-2 issue).

test_that("inventory enumerates every custom and helper family", {
  inventory_names <- names(cross_layer_inventory)
  expect_setequal(
    inventory_names,
    c(cross_layer_custom_families, cross_layer_helper_families)
  )
})

test_that("custom families expose all three brms callbacks", {
  for (name in cross_layer_custom_families) {
    entry <- cross_layer_inventory[[name]]
    fam <- entry$family_fun()
    expect_true(inherits(fam, "brmsfamily"), label = name)
    expect_true(is.function(fam$log_lik), label = paste(name, "log_lik"))
    expect_true(
      is.function(fam$posterior_predict),
      label = paste(name, "posterior_predict")
    )
    expect_true(
      is.function(fam$posterior_epred),
      label = paste(name, "posterior_epred")
    )
  }
})

test_that("every callback of every custom family has declared coverage", {
  for (name in cross_layer_custom_families) {
    entry <- cross_layer_inventory[[name]]
    for (cb in c("ll3", "pp3", "ep3")) {
      status <- entry[[cb]]$status
      expect_false(is.null(status), label = paste(name, cb))
      expect_true(
        status %in% c("ok", "warning", "unsupported", "skip"),
        info = paste(name, cb, "has invalid status", status)
      )
      if (status %in% c("unsupported", "skip")) {
        expect_true(
          !is.null(entry[[cb]]$reason) &&
            nzchar(entry[[cb]]$reason),
          info = paste(name, cb, "needs a documented reason")
        )
      }
    }
  }
})

test_that("inventory declares densities, RNGs and Stan functions correctly", {
  for (entry in cross_layer_inventory) {
    if (!is.null(entry$d_fun)) {
      expect_true(is.function(entry$d_fun), label = entry$name)
    }
    expect_true(is.function(entry$r_fun), label = entry$name)
    if (entry$custom) {
      expect_false(
        is.na(entry$stan_lpdf),
        info = paste(entry$name, "must declare its injected Stan lpdf name")
      )
      scode <- cross_layer_stan_scode(entry)
      expect_match(
        scode,
        sprintf("%s\\(", entry$stan_lpdf),
        info = paste(entry$name, "lpdf missing from injected scode")
      )
      if (!is.na(entry$stan_rng)) {
        expect_match(
          scode,
          sprintf("%s\\(", entry$stan_rng),
          info = paste(entry$name, "rng missing from injected scode")
        )
      }
    } else {
      expect_true(
        is.null(entry$family_fun),
        info = paste(entry$name, "is not a custom brms family")
      )
    }
  }
})

test_that("helper families accompany built-in brms families without callbacks", {
  for (name in cross_layer_helper_families) {
    entry <- cross_layer_inventory[[name]]
    expect_false(entry$custom, label = name)
    expect_null(entry$family_fun, info = paste(name, "family_fun"))
    expect_true(is.na(entry$stan_lpdf), info = paste(name, "stan_lpdf"))
  }
})

test_that("parities are covered by at least the representative families", {
  parity_families <- names(Filter(
    function(e) isTRUE(e$parity),
    cross_layer_inventory
  ))
  expect_gte(length(parity_families), 3)
  expect_true("lomax" %in% parity_families)
  expect_true("generalized_normal" %in% parity_families)
  expect_true("generalized_gamma" %in% parity_families)
})
