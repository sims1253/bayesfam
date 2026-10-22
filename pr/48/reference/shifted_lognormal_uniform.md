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
#> Loading 'brms' package (version 2.23.2). Useful instructions
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
#> Start sampling
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 4.4e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.44 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.225 seconds (Warm-up)
#> Chain 1:                0.177 seconds (Sampling)
#> Chain 1:                0.402 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 2.9e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.29 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.246 seconds (Warm-up)
#> Chain 2:                0.201 seconds (Sampling)
#> Chain 2:                0.447 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 2.9e-05 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.29 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 0.249 seconds (Warm-up)
#> Chain 3:                0.157 seconds (Sampling)
#> Chain 3:                0.406 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 2.6e-05 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.26 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 0.229 seconds (Warm-up)
#> Chain 4:                0.179 seconds (Sampling)
#> Chain 4:                0.408 seconds (Total)
#> Chain 4: 
plot(fit_mix)
```
