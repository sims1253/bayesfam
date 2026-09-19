# Layer 2 of the cross-layer suite (issue #28): parity between the R log
# densities and the Stan lpdfs injected as strings in the family R files.
#
# The injected functions of all inventory families with parity = TRUE are
# compiled into one cmdstanr program and evaluated in generated quantities at
# the per-family parity grids. Without a Stan toolchain the test skips cleanly
# (the comparison logic stays written so CI with CmdStan runs it).
#
# Parity for shifted_lognormal_uniform compares Stan against the mathematically
# correct reference, not against the R density, which is wrong until #24
# (wave-2 PR) lands on this branch.

test_that("injected Stan lpdfs match the R log densities", {
  skip_if_no_stan_toolchain()

  entries <- Filter(
    function(e) isTRUE(e$parity) && !is.na(e$stan_lpdf),
    cross_layer_inventory
  )
  program <- build_parity_program()
  model <- compile_cross_layer_stan(program, "cross-layer-parity")

  data <- list()
  for (entry in entries) {
    fam_tag <- gsub(".", "_", entry$name, fixed = TRUE)
    grid <- entry$parity_grid
    data[[paste0("M_", fam_tag)]] <- length(grid$y)
    data[[paste0("y_", fam_tag)]] <- grid$y
    for (an in names(grid)[-1]) {
      data[[paste0(an, "_", fam_tag)]] <- grid[[an]]
    }
  }
  ll_vars <- vapply(
    entries,
    function(e) paste0("ll_", gsub(".", "_", e$name, fixed = TRUE)),
    character(1)
  )
  fit <- run_fixed_param(model, data = data)
  expect_s3_class(fit, "CmdStanMCMC")
  total_points <- 0

  refs <- parity_reference()
  for (k in seq_along(entries)) {
    entry <- entries[[k]]
    fam_tag <- gsub(".", "_", entry$name, fixed = TRUE)
    stan_ll <- stan_gq(fit, paste0("ll_", fam_tag))
    expect_length(stan_ll, length(refs[[k]]$ll))
    total_points <- total_points + length(stan_ll)
    diff <- abs(stan_ll - refs[[k]]$ll)
    expect_true(
      all(diff < 1e-6),
      info = sprintf(
        "%s: max |Stan - R| log density difference %.3g at grid point %s",
        entry$name,
        max(diff),
        which.max(diff)
      )
    )
  }
  expect_gte(total_points, 50)
  expect_identical(length(ll_vars), length(entries))
})

test_that("injected scode of every custom family is extractable", {
  for (name in cross_layer_custom_families) {
    entry <- cross_layer_inventory[[name]]
    scode <- cross_layer_stan_scode(entry)
    expect_match(scode, entry$stan_lpdf, info = name)
  }
})
