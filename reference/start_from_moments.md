# Where a Univariate Family Starts, From the Data

Returns `n_start` starting values whose first is computed from the
response and whose rest are random draws. It is the method that both
`continuous_distrib` and `discrete_distrib` register, so a univariate
family uses it unless the family registers a method of its own, as
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md),
[`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md),
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
and
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
do.

The data-based value comes from
[`moment_estimates()`](https://statmodels7.github.io/distributions7/reference/moment_estimates.md)
where the family has an entry there. The other families fall back to
reading `params_interpretation`: a parameter meaning a location is
started at the sample median, one meaning a spread at the sample
standard deviation or its square, one meaning degrees of freedom at what
the sample kurtosis implies, and a shape, a dispersion or a probability
keeps its draw, those being of order one whatever the data.

## Usage

``` r
start_from_moments(distrib, y, n_start = 5L, ...)
```

## Arguments

- distrib:

  A univariate distribution, supplying `params`, `params_interpretation`
  and `params_bounds`.

- y:

  The response, a numeric vector. Non-finite entries are dropped, and a
  sample with fewer than two usable values gets the random draws alone.

- n_start:

  How many starting values are wanted, a single positive integer.
  Defaults to 5, and a value below 1 is raised to 1.

- ...:

  Unused.

## Value

A list of `n_start` named parameter lists on the parameter scale, the
first from the data where one could be computed.

## Why the base method is not used for univariate families

The base method draws each parameter from its own domain and does not
read `y`. That is adequate while the response is of order one and fails
when it is not: on a response of mean 919 and standard deviation 169 the
draws are of order one, the first Newton step is taken where the
residuals are hundreds of standard deviations wide, and the scale can
run to the largest representable double. A start computed from the data
is on the scale of the data.

## Degrees of freedom from the kurtosis

A \\t\\ of \\\nu\\ degrees has excess kurtosis \\6/(\nu-4)\\, so a
sample kurtosis above 0.05 inverts to \\\nu = 6/\hat\gamma_2 + 4\\,
capped at 100; a sample no heavier-tailed than a Gaussian starts at 30.
A large start keeps a `student_t2`, whose scale parameter is the
standard deviation, away from the degenerate ridge at its lower bound of
\\\nu = 2\\, where the scale tends to infinity as \\\nu\\ falls; a run
started at a small \\\nu\\ can follow the ridge down to the bound.

## Bounds and the remaining starts

Every value is moved strictly inside its parameter's bounds before it is
returned, because a sample median can land exactly on the boundary of a
support and the validator treats bounds as open. Only the first start is
data-based; the rest stay random, so a caller requesting several still
explores, and
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
reaches them when the first fails.

## See also

[`moment_estimates()`](https://statmodels7.github.io/distributions7/reference/moment_estimates.md)
for the family-by-family inversions;
[`distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.md)
for the generic;
[`clamp_to_bounds()`](https://statmodels7.github.io/distributions7/reference/clamp_to_bounds.md)
for the boundary rule of the moment estimates;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
which consumes the list.
