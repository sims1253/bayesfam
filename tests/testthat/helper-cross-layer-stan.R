# Stan toolchain helpers for the cross-layer test suite (issues #28 / #29).
#
# Layer 2 tests extract the Stan functions injected in the family R files,
# compile them with cmdstanr (when a CmdStan toolchain is available), and
# evaluate them in generated quantities. Without a toolchain the tests skip
# cleanly with the comparison logic fully written so they light up in CI.

skip_if_no_stan_toolchain <- function() {
  testthat::skip_if_not_installed("cmdstanr")
  has_cmdstan <- tryCatch(
    {
      cmdstanr::cmdstan_path()
      TRUE
    },
    error = function(e) FALSE
  )
  if (!has_cmdstan) {
    testthat::skip("Stan toolchain unavailable (no CmdStan installation found)")
  }
  invisible(TRUE)
}

# Extract the injected functions-block scode of a family as a single string.
cross_layer_stan_scode <- function(entry) {
  fam <- entry$family_fun()
  scodes <- vapply(fam$stanvars, function(v) v$scode, character(1))
  paste(scodes, collapse = "\n")
}

# Compile a Stan program, caching the executable in the session temp dir so
# repeated / parameterized runs do not recompile.
compile_cross_layer_stan <- function(program, name) {
  dir <- file.path(tempdir(), "bayesfam-cross-layer-stan")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  stan_file <- file.path(dir, paste0(name, ".stan"))
  writeLines(program, stan_file)
  cmdstanr::cmdstan_model(stan_file, dir = dir)
}

# Run a parameter-free Stan program in fixed_param mode and return the fit.
run_fixed_param <- function(
  model,
  data,
  seed = 20240928,
  iter_sampling = 1
) {
  model$sample(
    data = data,
    chains = 1,
    iter_warmup = 0,
    iter_sampling = iter_sampling,
    fixed_param = TRUE,
    seed = seed,
    refresh = 0,
    show_messages = FALSE,
    show_exceptions = FALSE
  )
}

# Extract one generated-quantity variable as a flat numeric vector (vector
# GQs come back indexed as variable[1], variable[2], ... and are flattened in
# index order).
stan_gq <- function(fit, variable) {
  as.vector(fit$draws(variables = variable))
}

# Build one combined Stan program evaluating the injected lpdf of every
# inventory family with parity = TRUE at its parity grid. Argument vectors are
# named y_<family> / <arg>_<family>; results come back as ll_<family>.
build_parity_program <- function() {
  entries <- Filter(
    function(e) isTRUE(e$parity) && !is.na(e$stan_lpdf),
    cross_layer_inventory
  )
  data_lines <- character(0)
  gq_lines <- character(0)
  for (entry in entries) {
    fam_tag <- gsub(".", "_", entry$name, fixed = TRUE)
    grid <- entry$parity_grid
    m <- length(grid$y)
    data_lines <- c(
      data_lines,
      sprintf("int<lower=1> M_%s;", fam_tag),
      sprintf("vector[M_%s] y_%s;", fam_tag, fam_tag),
      vapply(
        names(grid)[-1],
        function(an) {
          sprintf("vector[M_%s] %s_%s;", fam_tag, an, fam_tag)
        },
        character(1)
      )
    )
    arg_expr <- vapply(
      names(grid)[-1],
      function(an) sprintf("%s_%s[m]", an, fam_tag),
      character(1)
    )
    call <- sprintf(
      "%s(y_%s[m] | %s)",
      entry$stan_lpdf,
      fam_tag,
      paste(arg_expr, collapse = ", ")
    )
    gq_lines <- c(
      gq_lines,
      sprintf("vector[M_%s] ll_%s;", fam_tag, fam_tag),
      sprintf(
        "for (m in 1:M_%s) ll_%s[m] = %s;",
        fam_tag,
        fam_tag,
        call
      )
    )
  }
  paste0(
    "functions {\n",
    paste(
      vapply(entries, cross_layer_stan_scode, character(1)),
      collapse = "\n"
    ),
    "\n}\ndata {\n",
    paste(data_lines, collapse = "\n"),
    "\n}\nmodel {\n}\ngenerated quantities {\n",
    paste(gq_lines, collapse = "\n"),
    "\n}\n"
  )
}

# R-side reference values for the combined parity program, one list per family
# entry with `ll` = per-point R log density.
parity_reference <- function() {
  entries <- Filter(
    function(e) isTRUE(e$parity) && !is.na(e$stan_lpdf),
    cross_layer_inventory
  )
  lapply(entries, function(entry) {
    grid <- entry$parity_grid
    m <- length(grid$y)
    ll <- vapply(seq_len(m), function(j) {
      point <- lapply(grid[-1], function(v) v[j])
      if (!is.null(entry$parity_ref)) {
        return(entry$parity_ref(c(list(y = grid$y[j]), point)))
      }
      # grid parameters follow dpars order; densities may use different
      # argument names (e.g. ndt -> shift for shifted_inv_gaussian)
      if (length(point) == length(entry$d_argnames)) {
        names(point) <- entry$d_argnames
      }
      cross_layer_r_loglik(entry, grid$y[j], point)
    }, numeric(1))
    list(name = entry$name, ll = ll)
  })
}
