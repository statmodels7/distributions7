# Fix Parameters of a Distribution at Known Values

Returns the distribution obtained by holding some parameters of
`distrib` at known values, leaving the others to be supplied and
estimated. It is the only wrapper in this package that REMOVES
parameters, and it derives nothing: the density is the parent's at the
reassembled vector and the derivatives are the parent's components among
the free indices.

## Usage

``` r
fixed(distrib, ...)
```

## Arguments

- distrib:

  The distribution whose parameters are to be fixed, inheriting from
  `continuous_distrib`, `discrete_distrib` or `multivariate_distrib`.

- ...:

  The fixed values, named after the parameters they fix, as in
  `fixed(gaussian1_distrib(), mu = 0)`. Each must be a single finite
  number strictly inside its parameter's domain, and each name must be a
  parameter of `distrib`. A name that is not, a value outside the
  domain, a value that is not a single number, and an empty `...` are
  each rejected with an error saying which condition failed.

## Value

An S7 object of class
[FixedContinuousDistrib](https://statmodels7.github.io/distributions7/reference/FixedContinuousDistrib.md),
[FixedDiscreteDistrib](https://statmodels7.github.io/distributions7/reference/FixedDiscreteDistrib.md)
or
[FixedMultivariateDistrib](https://statmodels7.github.io/distributions7/reference/FixedMultivariateDistrib.md),
matching the parent's branch. Its `params` are the parent's less the
fixed names in the parent's order, `n_params` falls by as many,
`fixed_params` holds the values, and `distrib_name` is the parent's with
the held values in brackets.

## The construction

Splitting the parent's parameters into a fixed part \\\theta_C = c\\ and
a free part \\\theta_F\\, \$\$f\_{\mathrm{fix}}(y; \theta_F) = f(y;
\theta_F, c), \qquad \ell^{(i_1 \cdots i_k)}\_{\mathrm{fix}} =
\ell^{(i_1 \cdots i_k)}, \quad i_1, \dots, i_k \in F.\$\$ `theta`
carries only the free parameters, every generic answers as the parent
does at the full vector, and a derivative is the parent's restricted to
the free indices: a subvector of the score, a submatrix of the Hessian,
sub-arrays at orders three and four. Nothing is recomputed, no
normalizing constant changes, and
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
estimates the free parameters with standard errors and intervals for
them alone.

## What is accepted

Fixed values are single finite numbers, strictly inside the OPEN domain
of their parameter. Fixing a parameter of a distribution that is already
a fixed-parameter wrapper collapses the two into one wrapper around the
original parent. Fixing a WRAPPER's own parameter is allowed and useful:
`fixed(zero_inflated(d), zi = 0.3)` is a zero-inflated model with a
known inflation rate. Fixing every parameter is allowed too and gives a
fully known distribution with an empty parameter set. Calling with no
named value is an error: the result would be the parent unchanged, and
returning it silently would hide a missing argument.

## What it is for

A prior. `fixed(gaussian1_distrib(), mu = 0)` is the ridge penalty with
its scale free, `fixed(laplace2_distrib(), mu = 0)` is the lasso, and
`fixed(mvgaussian1_distrib(p), mu1 = 0, ...)` is what a random effect is
distributed by. `fixed(folded(gaussian1_distrib()), mu = 0)` is the
half-normal.

The per-parameter smoothness declaration travels with the free
parameters, so fixing the location of a Laplace leaves a distribution
whose remaining parameter is smooth.

## Notation

\\f\\ is the parent's density, \\\theta_C = c\\ the fixed parameters,
\\\theta_F\\ the free ones, \\F\\ their index set and \\\ell^{(i_1\cdots
i_k)}\\ a derivative of the log-density in the parameters named.

## See also

[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md),
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md),
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md)
and
[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md)
for the other wrappers, and
[FixedContinuousDistrib](https://statmodels7.github.io/distributions7/reference/FixedContinuousDistrib.md)
for what the methods do.

## Examples

``` r
# A gaussian with a known mean: only sigma remains.
d <- fixed(gaussian1_distrib(), mu = 0)
d@params
#> [1] "sigma"
theta <- list(sigma = 2)
distrib_pdf(d, c(-1, 0, 1), theta)
#> [1] 0.1760327 0.1994711 0.1760327

# The score is the corresponding component of the parent's, unchanged.
full <- distrib_gradient(gaussian1_distrib(), c(-1, 0, 1),
                         list(mu = 0, sigma = 2))
all.equal(distrib_gradient(d, c(-1, 0, 1), theta)$sigma, full$sigma)
#> [1] TRUE

# The lasso prior: a Laplace in its rate, centered at zero.
lasso <- fixed(laplace2_distrib(), mu = 0)
lasso@params
#> [1] "lambda"
b <- c(-1, 0.5, 2)
all.equal(-distrib_pdf(lasso, b, list(lambda = 2), log = TRUE),
          2 * abs(b) - log(2 / 2))
#> [1] TRUE

# A wrapper's own parameter can be held.
fixed(zero_inflated(poisson_distrib()), zi = 0.3)@params
#> [1] "mu"

# Two calls collapse into one wrapper, and fixing everything is legal.
fixed(fixed(gaussian1_distrib(), mu = 0), sigma = 1)@fixed_params
#> $mu
#> [1] 0
#> 
#> $sigma
#> [1] 1
#> 

# Four refusals, each naming the condition that failed.
try(fixed(gaussian1_distrib(), nope = 1))
#> Error : 'nope' is not a parameter of 'gaussian1'. Parameters: mu, sigma.
try(fixed(gaussian1_distrib(), sigma = -1))
#> Error : The value fixing 'sigma' (-1) is outside its open domain (0, Inf).
try(fixed(gaussian1_distrib(), mu = c(0, 1)))
#> Error : The value fixing 'mu' must be a single finite number.
try(fixed(gaussian1_distrib()))
#> Error : fixed() needs at least one named value, as in fixed(d, mu = 0). With
#>   none, the result would be the parent distribution unchanged; returning
#>   it silently would hide a missing argument rather than report it.
```
