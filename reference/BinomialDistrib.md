# Binomial Distribution Class

The S7 class of the binomial family: the number of successes in `size`
independent trials, each succeeding with probability \\\mu \in (0, 1)\\.
It inherits from `discrete_distrib`, so expectations over its support
are sums and no derivative with respect to the response is defined.

The class adds one property to the parent's, `size`, and that property
is the reason it is documented separately from
[`bernoulli_distrib()`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md):
**`size` is a constant of the object, not a parameter**. It is not
estimated, it carries no link and no bound, and it does not appear in
`params`. It may be a single number or one value per observation, which
is how grouped binary data with unequal group sizes is described.

The default link is the logit, the **canonical** link here: on its scale
the observed and the expected information coincide.

Build one with
[`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md).
This page documents the raw S7 constructor, which takes the parent's
properties plus `size` and validates none of the relationships between
them.

## Usage

``` r
BinomialDistrib(
  distrib_name = character(0),
  dimension = character(0),
  bounds = integer(0),
  params = character(0),
  params_interpretation = character(0),
  n_params = integer(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0),
  size = integer(0)
)
```

## Arguments

- distrib_name:

  A single character string specifying the name of the distribution
  (e.g., `"student t"`).

- dimension:

  A character string indicating the dimensionality (`"univariate"` or
  `"multivariate"`).

- bounds:

  A numeric vector of length 2 defining the overall support of the
  distribution `c(lower, upper)`.

- params:

  A character vector containing the names of the distribution parameters
  (e.g., `c("mu", "sigma")`).

- params_interpretation:

  A character vector (typically named) providing the statistical
  interpretation of each parameter (e.g., `c(mu = "location")`).

- n_params:

  A numeric value specifying the total number of parameters.

- params_bounds:

  A list of numeric vectors of length 2, specifying the valid
  mathematical domain for each individual parameter.

- link_params:

  A list of link function objects corresponding to each parameter,
  primarily used to map parameters to the unconstrained real line for
  optimization algorithms.

- params_smooth:

  An optional named logical vector flagging, for each parameter, whether
  the log-likelihood is differentiable with respect to it. Defaults to
  all `TRUE` (leave empty). Set an entry to `FALSE` for parameters at
  which the log-likelihood has a kink (e.g. the location of a Laplace
  distribution): the observed Hessian is then degenerate and the
  expected information must be obtained from the score variance rather
  than from \\-\mathbb{E}\[H\]\\ (see
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md)).

- size:

  A numeric vector, the number of trials. A single value applies to
  every observation; a vector of the length of the response gives one
  count of trials per observation. Not a parameter: it is fixed data,
  and the object stores it in the `size` property.

## Value

An S7 object of class `BinomialDistrib`, inheriting from
`discrete_distrib` and from `distrib`. Its properties are the parent's
(`distrib_name`, `dimension`, `bounds`, `params`,
`params_interpretation`, `n_params`, `params_bounds`, `link_params`,
`params_smooth`) plus `size`. For an object built by
[`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md)
they hold `"binomial"`, `"univariate"`, `c(0, max(size))`, `"mu"`,
`c(mu = "probability")`, `1`, the domain \\(0, 1)\\, the one link, and
the `size` given.

## Methods

Registered on this class:
[`distrib_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.BinomialDistrib.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.BinomialDistrib.md),
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.BinomialDistrib.md),
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.BinomialDistrib.md),
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.BinomialDistrib.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.BinomialDistrib.md),
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.BinomialDistrib.md),
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.BinomialDistrib.md),
[`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.BinomialDistrib.md)

Everything else is inherited from
[`discrete_distrib()`](https://statmodels7.github.io/distributions7/reference/discrete_distrib.md).

## See also

[`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md)
to build one;
[`bernoulli_distrib()`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md),
the case `size = 1`;
[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md)
and
[`betabinom2_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md)
when the probability varies between groups;
[`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md),
the limit of many trials at a small probability.

## Examples

``` r
d <- binomial_distrib(size = 10)
S7::S7_inherits(d, discrete_distrib)
#> [1] TRUE

# size is a property of the object; only mu is a parameter.
d@size
#> [1] 10
d@params
#> [1] "mu"
d@bounds
#> [1]  0 10

# The mean is n p and the variance n p (1-p).
th <- list(mu = 0.3)
c(mean = mean(d, th), var = variance(d, th), n_p = 10 * 0.3)
#> mean  var  n_p 
#>  3.0  2.1  3.0 

# size = 1 is the Bernoulli, exactly.
all.equal(distrib_pdf(binomial_distrib(size = 1), c(0, 1), th),
          distrib_pdf(bernoulli_distrib(), c(0, 1), th))
#> [1] TRUE
```
