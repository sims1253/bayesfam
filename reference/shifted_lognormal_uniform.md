# A mixture of shifted lognormal and uniform distribution suitable for modelling reaction times.

A contaminated response time distribution. The mixture can be described
as \$\$y_i = \begin{cases} u_i & \mathrm{if} \quad z_i = 0 \\ p_i s_i +
r_i & \mathrm{if} \quad z_i = 1 \end{cases} \\ u_i \sim Uniform(0,
\alpha) \\ \log(r_i) \sim Normal(\mu_i, \sigma) \\ P(z_i = 0) = \theta
\\ 0 \< p_i \< 1 \$\$ Here \\\mu, \sigma, \theta\\ (`mix`) and \\p\\
(`shiftprop`) are estimated, whereas \\s_i\\ (`max_shift`) and
\\\alpha\\ (`max_uniform`) are given as data via `vreal()`.

## Usage

``` r
shifted_lognormal_uniform(
  link = "identity",
  link_sigma = "log",
  link_mix = "logit",
  link_shiftprop = "logit"
)
```

## Source

Idea by Nathaniel Haines (https://twitter.com/Nate\_\_Haines), code by
Martin Modrák. Some background, discussion and examples at
<http://www.martinmodrak.cz/2021/04/01/using-brms-to-model-reaction-times-contaminated-with-errors/>

## Arguments

- link:

  Link function for the location parameter (default: "identity")

- link_sigma:

  Link function for the scale parameter (default: "log")

- link_mix:

  Link function for the mixture parameter (default: "logit")

- link_shiftprop:

  Link function for the shift proportion parameter (default: "logit")

## Details

Note that you cannot build this distribution with the built-in support
for mixtures in `brms`, because the uniform component is effectively a
zero-parameter distribution which cannot be expressed in `brms`.

## Examples

``` r
library(brms)
#> Loading required package: Rcpp
#> Loading 'brms' package (version 2.23.1). Useful instructions
#> can be found by typing help('brms'). A more detailed introduction
#> to the package is available through vignette('brms_overview').
#> 
#> Attaching package: ‘brms’
#> The following object is masked from ‘package:stats’:
#> 
#>     ar
set.seed(31546522)
# Bounds of the data
max_shift <- 0.3
shift <- runif(1) * max_shift
max_uniform <- 10
mix <- 0.1

# Generate parameters
N <- 100
Intercept <- 0.3
beta <- 0.5
X <- rnorm(N)
mu <- rep(Intercept, N) + beta * X
sigma <- 0.5

rt <- rshifted_lognormal_uniform(N, meanlog = mu, sdlog = sigma, mix = mix,
                                 shift = shift, max_uniform = max_uniform)

dd <- data.frame(rt = rt, x = X,
                 max_shift = max_shift, max_uniform = max_uniform)

fam <- shifted_lognormal_uniform()
fit_mix <- brm(rt | vreal(max_shift, max_uniform) ~ x, data = dd, family = fam,
               stanvars = fam$stanvars,
               prior = c(prior(beta(1, 5), class = "mix")))
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')
plot(fit_mix)
#> Error: object 'fit_mix' not found
```
