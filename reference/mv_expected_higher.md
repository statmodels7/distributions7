# Expected Higher Derivatives of a Simplex-Valued Family by Sampling

Averages the closed-form observed third or fourth derivatives over a
sample drawn from the family. This is the route the multivariate branch
takes throughout: a quadrature over the simplex has no counterpart to
the one-dimensional split at quantiles, and the Bartlett route would
need the score's own higher derivatives. The result therefore carries
Monte Carlo error of order `nsim^(-1/2)`, and a caller who needs it
reproducible sets a seed first.

## Usage

``` r
mv_expected_higher(distrib, y, theta, order, nsim)
```

## Arguments

- distrib:

  A `DirichletDistrib` or `MultinomialDistrib` object.

- y:

  The observed response, read only for its number of rows.

- theta:

  A named list of parameters, each component a single number.

- order:

  The derivative order, `3L` or `4L`.

- nsim:

  The Monte Carlo sample size, a single positive whole number.

## Value

A named list of numeric vectors of length `n_obs(distrib, y)`, each
constant across observations, keyed as
`deriv_names(distrib@params, order)`.

## Details

Only the two simplex-valued families reach here. The branch is taken on
whether `distrib` inherits `DirichletDistrib`, everything else being
routed to
[`multinomial_higher()`](https://statmodels7.github.io/distributions7/reference/multinomial_higher.md),
so this function must not be called for any third family.

## See also

[`dirichlet_higher()`](https://statmodels7.github.io/distributions7/reference/dirichlet_higher.md)
and
[`multinomial_higher()`](https://statmodels7.github.io/distributions7/reference/multinomial_higher.md)
for the observed derivatives it averages, and
[`distrib_deriv3.DirichletDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.DirichletDistrib.md)
for the method that calls it.

## Examples

``` r
d <- dirichlet_distrib(3)
theta <- list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8)
set.seed(1)
y <- distrib_rng(d, 5, theta)

set.seed(5)
e3 <- distributions7:::mv_expected_higher(d, y, theta, 3L, 3000)

# One value per observation, and every vector constant.
lengths(e3)[1:3]
#> mean_alr1_mean_alr1_mean_alr1 mean_alr1_mean_alr1_mean_alr2 
#>                             5                             5 
#>       mean_alr1_mean_alr1_phi 
#>                             5 
all(vapply(e3, function(z) diff(range(z)) == 0, TRUE))
#> [1] TRUE

# Two seeds give two answers, this being a sample rather than an integral.
set.seed(6)
c(first = e3[[1]][1],
  second = distributions7:::mv_expected_higher(d, y, theta, 3L, 3000)[[1]][1])
#>      first     second 
#> -0.7331362 -0.6983632 
```
