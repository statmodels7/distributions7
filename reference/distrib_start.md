# A Starting Value Drawn from the Data

Returns a starting value for
[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
computed from the response rather than guessed.

## Usage

``` r
distrib_start(distrib, y, n_start = 5L, ...)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  The response.

- n_start:

  How many starting values are wanted. A method free to supply only one
  may ignore it.

- ...:

  Passed to methods.

## Value

A list of named parameter lists, on the **parameter** scale.

## Details

A starting value that ignores the data is a starting value that can be
arbitrarily far from the answer, and how far decides whether the fit
takes one step or never arrives. A four-dimensional gaussian fitted to
the iris measurements is the plain case: started at the origin of the
unconstrained scale, which is a unit covariance and a zero mean, Newton
with the expected information spends five hundred iterations and stops
at a log-likelihood of \\-836\\; started at the sample mean and the
sample covariance it converges in one iteration to \\-379.9146\\, which
is the exact maximum. Nothing about the arithmetic changed.

The default method returns random parameters, as before, so a
distribution that says nothing loses nothing. A family that can do
better says so by registering a method: an exact maximum likelihood
estimator where one is known, a method-of-moments estimator otherwise,
or the estimate of a simpler family the harder one contains.

A method may return several starting values, as a list, and
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
will try each; the first is the one it prefers.

## See also

[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
[`generate_random_theta`](https://statmodels7.github.io/distributions7/reference/generate_random_theta.md)

## Examples

``` r
set.seed(1)
d <- mvgaussian_distrib(2)
y <- distrib_rng(d, 200, list(mu1 = 1, mu2 = -1, sigma_log_L1 = 0,
                              sigma_log_L2 = 0, sigma_L2.1 = 0.5))

# the gaussian knows its own maximum likelihood estimate, so the fit starts
# there and has nothing left to do
start <- distrib_start(d, y)[[1]]
mv_location(d, start)
#>         v1         v2 
#>  1.0355396 -0.9415925 
colMeans(y)
#>         v1         v2 
#>  1.0355396 -0.9415925 
```
