# bayesfam

This is a collection of custom brms families initially written for our [simulation study ](https://arxiv.org/abs/2210.06927) and [bayesim](https://github.com/sims1253/bayesim), our simulation framework.

They can be used like any other brms family. The only thing to remember is to explicitly hand brm the `stanvar`
value that is part of the custom family object.
```r
library(brms)
library(bayesfam)
data <- list(y = rbetaprime(1000, mu = exp(2.3), phi = 2))
brm(
    y ~ 1 ,
    data = data,
    family = betaprime(),
    stanvars = betaprime()$stanvars, # This is the important part to remember
  )

```

```
Family: betaprime 
  Links: mu = log; phi = identity 
Formula: y ~ 1 
   Data: data (Number of observations: 1000) 
  Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
         total post-warmup draws = 4000

Population-Level Effects: 
          Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
Intercept     2.29      0.02     2.25     2.33 1.00     1329     2061

Family Specific Parameters: 
    Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
phi     2.12      0.18     1.78     2.48 1.00     1275     1828

Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```

