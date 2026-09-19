# Probability density and RNG function for the mixture of shifted lognormal and uniform distributions.

Probability density and RNG function for the mixture of shifted
lognormal and uniform distributions.

## Usage

``` r
rshifted_lognormal_uniform(
  n,
  meanlog = 0,
  sdlog = 1,
  mix = 0.1,
  shift = 0,
  max_uniform = 100
)

dshifted_lognormal_uniform(
  y,
  meanlog = 0,
  sdlog = 1,
  mix = 0.1,
  shift = 0,
  max_uniform = 100
)
```

## Source

Some background, discussion and examples at
<http://www.martinmodrak.cz/2021/04/01/using-brms-to-model-reaction-times-contaminated-with-errors/>

## Arguments

- n:

  the number of values to draw from the RNG

- meanlog:

  the mean of the lognormal component

- sdlog:

  the sd of the lognormal component

- mix:

  the probability of the value comming from the uniform component

- shift:

  the shift the lognormal distribution

- max_uniform:

  the maximum value of the uniform component

- y:

  the observed value

## Details

The mixture of shifted lognormal and uniform can be described as \$\$y_i
= \begin{cases} u_i & \mathrm{if} \quad z_i = 0 \\ s_i + r_i &
\mathrm{if} \quad z_i = 1 \end{cases} \\ u_i \sim Uniform(0, \alpha) \\
\log(r_i) \sim Normal(\mu_i, \sigma) \\ P(z_i = 0) = \theta\$\$

Where θ corresponds to `mix`, α to `max_uniform` and \\s_i\\ to `shift`.
