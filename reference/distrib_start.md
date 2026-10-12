# A Starting Value Drawn From the Data

Returns the starting values that
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
begins from, computed from the response where the family has a method
for it. The result is a list of named parameter lists on the
**parameter** scale. The fit tries them in order until two converged
runs agree on the best objective, so the first element is the one tried
first.

## Usage

``` r
distrib_start(distrib, y, n_start = 5L, ...)
```

## Arguments

- distrib:

  An object inheriting from `distrib`.

- y:

  The response: a numeric vector, or an \\n \times p\\ matrix for a
  multivariate family. A method may ignore it, and the base one does.

- n_start:

  How many starting values are wanted, a single positive integer.
  Defaults to 5. A method that returns its own estimate returns one and
  ignores this, which the two multivariate methods and the two
  Poisson-inverse Gaussian ones do.

- ...:

  Passed to methods. The shipped methods do not use it.

## Value

A list of named parameter lists on the parameter scale, each complete:
one component per entry of `distrib@params`, named and ordered as
`distrib@params`. Every value is strictly inside its parameter's bounds,
which the validator treats as open.

## Details

A starting value that ignores the data can be arbitrarily far from the
answer, and the distance decides how much work the fit does. For a
multivariate Gaussian the sample mean and the sample covariance are the
maximum likelihood estimate, so a fit that starts there converges in one
iteration.

The base method draws each parameter from its own domain and never reads
`y`, so a family that registers nothing loses nothing. A family with a
better estimator registers a method: an exact maximum likelihood
estimator where one is known, a method-of-moments estimator otherwise,
or the estimate of a simpler family the harder one contains. The two
univariate methods route through
[`start_from_moments()`](https://statmodels7.github.io/distributions7/reference/start_from_moments.md),
which uses the family's own moment inversion where
[`moment_estimates()`](https://statmodels7.github.io/distributions7/reference/moment_estimates.md)
has one and reads `params_interpretation` where it does not.

## See also

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
which calls this;
[`start_from_moments()`](https://statmodels7.github.io/distributions7/reference/start_from_moments.md)
for the univariate route;
[`moment_estimates()`](https://statmodels7.github.io/distributions7/reference/moment_estimates.md)
for the family-by-family inversions;
[`generate_random_theta()`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.md)
for the draw the base method makes.

## Examples

``` r
set.seed(1)
d <- mvgaussian1_distrib(2)
y <- distrib_rng(d, 200, list(mu1 = 1, mu2 = -1, sigma_log_L1 = 0,
                              sigma_log_L2 = 0, sigma_L2.1 = 0.5))

# An unstructured covariance has a closed-form estimate, so the fit starts
# at the answer and has nothing left to do.
start <- distrib_start(d, y)[[1]]
rbind(start = mv_location(d, start), sample = colMeans(y))
#>             v1         v2
#> start  1.03554 -0.9415925
#> sample 1.03554 -0.9415925

# A univariate family gives one start from the data and the rest at
# random, so a caller requesting several still explores.
set.seed(2)
s <- distrib_start(gaussian1_distrib(), rnorm(200, 900, 170), n_start = 3)
length(s)
#> [1] 3
vapply(s, function(th) unlist(th), numeric(2))
#>           [,1]     [,2]     [,3]
#> mu    900.9564 0.439958 2.663036
#> sigma 184.9586 1.091506 5.666502
```
