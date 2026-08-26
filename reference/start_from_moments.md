# Where a Univariate Family Starts, From the Data

Returns `n_start` starting values whose first is computed from the
response and whose rest are random draws. It is the method both
`continuous_distrib` and `discrete_distrib` register, so it is what
every univariate family in the package uses.

The data-based value comes from
[`moment_estimates()`](https://statmodels7.github.io/distributions7/reference/moment_estimates.md)
where the family has an entry there, which 37 of the 42 univariate
families do. The other five fall back to reading
`params_interpretation`: a parameter meaning a location is started at
the sample median, one meaning a spread at the sample standard deviation
or its square, one meaning degrees of freedom at what the sample
kurtosis implies, and a shape, a dispersion or a probability keeps its
draw, those being of order one whatever the data.

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

## What this replaced, and why it was necessary

The base method draws each parameter from its own domain and never reads
`y`. That is adequate while the response is of order one and fails
completely when it is not: on a response of mean 919 and standard
deviation 169 the draws are of order one, the first Newton step is taken
where the residuals are hundreds of standard deviations wide, and the
scale runs to the largest representable double. Measured on a Gaussian,
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
recovers \\N(5, 2)\\ and \\N(50, 20)\\ and fails on \\N(500, 200)\\: the
defect was a threshold in the scale of the data, not in the family.

## Degrees of freedom from the kurtosis

A \\t\\ of \\\nu\\ degrees has excess kurtosis \\6/(\nu-4)\\, so a
sample kurtosis above 0.05 inverts to \\\nu = 6/\hat\gamma_2 + 4\\,
capped at 100; a sample no heavier-tailed than a Gaussian starts at 30.
Starting large matters. A `student_t2`, whose scale parameter is the
standard deviation, has a degenerate ridge at its lower bound of \\\nu =
2\\ where the scale runs to infinity as \\\nu\\ falls, and a run started
small slides down it: measured on 610 abdominal circumferences, whose
excess kurtosis is \\-1.09\\, the random draws put \\\nu\\ between 2.8
and 8 and the fit reached the boundary at a log-likelihood of
\\-3688.28\\, while any start of 2.5 or more with the location and scale
at their sample values reaches \\-3600.71\\, the Gaussian limit and the
true maximum.

## The clamp, and why only the first start

Every value is moved strictly inside its parameter's bounds before it is
returned, because a sample median can land exactly on the boundary of a
support and the validator treats bounds as open. Only the first start is
data-based; the rest stay random, so a caller asking for several still
explores, and
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
reaches them when the first fails.

## See also

[`moment_estimates()`](https://statmodels7.github.io/distributions7/reference/moment_estimates.md)
for the family-by-family inversions;
[`distrib_start()`](https://statmodels7.github.io/distributions7/reference/distrib_start.md)
for the generic;
[`clamp_to_bounds()`](https://statmodels7.github.io/distributions7/reference/clamp_to_bounds.md)
for the boundary rule;
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
which consumes the list.
