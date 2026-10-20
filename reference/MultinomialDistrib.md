# Multinomial Distribution Class

This is the S7 class of the multinomial family on \\p\\ categories and
\\n\\ trials, parametrized by a probability vector carried on a
`parameters7` simplex. It inherits from `multivariate_distrib`. Its
support is a finite set of points: every vector of non-negative integers
summing to \\n\\.

The class carries two properties beyond the parent's: `size`, the number
of trials, fixed at construction; and `param`, the simplex on which the
probabilities lie, whose free names become the family's own parameter
names prefixed by `probs_`. Build one with
[`multinomial_distrib()`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md),
which validates `n_dim`, `size` and the simplex and fills the properties
in. This page documents the raw S7 constructor, which runs only the
validators of `multivariate_distrib` and `distrib` and does not check
`size` or `param`.

## Usage

``` r
MultinomialDistrib(
  distrib_name = character(0),
  dimension = character(0),
  bounds = integer(0),
  params = character(0),
  params_interpretation = character(0),
  n_params = integer(0),
  params_bounds = list(),
  link_params = list(),
  params_smooth = logical(0),
  n_dim = integer(0),
  size = integer(0),
  param = NULL
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

- n_dim:

  The dimension \\p\\ of one observation, a single positive integer. The
  validator rejects any other value, and also a `dimension` property
  other than `"multivariate"`.

- size:

  The number of trials \\n\\, a single positive integer stored as a
  numeric. It belongs to the object, so an object cannot be reused
  across data sets whose trial counts differ.

- param:

  The `parameters7`
  [`parameters7::simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.html)
  carrying the probabilities. Its `n_free` is \\p-1\\, one fewer than
  the number of categories.

## Value

An S7 object of class `MultinomialDistrib`, inheriting from
`multivariate_distrib` and from `distrib`. Beyond `size` and `param` its
properties are the parent's. For an object built by
[`multinomial_distrib()`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md),
`params` is the simplex's free names prefixed by `probs_`, `n_params` is
\\p-1\\, `bounds` is `c(0, size)`, and every parameter carries an
identity link, because its free value is already unconstrained. The
family has no dispersion parameter.

## Methods

The following methods are registered on this class, and each is linked
to the page that documents it:

- the distribution:
  [`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.MultinomialDistrib.md)
  and
  [`distrib_rng()`](https://statmodels7.github.io/distributions7/reference/distrib_rng.MultinomialDistrib.md);

- the derivatives in the parameters:
  [`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.MultinomialDistrib.md),
  [`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.MultinomialDistrib.md),
  [`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.MultinomialDistrib.md),
  [`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.MultinomialDistrib.md)
  and
  [`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.MultinomialDistrib.md);

- the moments:
  [`mean()`](https://statmodels7.github.io/distributions7/reference/mean.MultinomialDistrib.md)
  and
  [`variance()`](https://statmodels7.github.io/distributions7/reference/variance.MultinomialDistrib.md);

- the multivariate summaries:
  [`mv_location()`](https://statmodels7.github.io/distributions7/reference/mv_location.MultinomialDistrib.md),
  [`mv_sigma()`](https://statmodels7.github.io/distributions7/reference/mv_sigma.MultinomialDistrib.md),
  [`mv_marginal()`](https://statmodels7.github.io/distributions7/reference/mv_marginal.MultinomialDistrib.md)
  and
  [`mv_support()`](https://statmodels7.github.io/distributions7/reference/mv_support.MultinomialDistrib.md).

Every other generic is inherited from
[`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md),
whose methods signal an error for the distribution function, the
quantile function, the response gradient and Hessian, the mixed response
derivative, the skewness and the kurtosis.

## See also

[`multinomial_distrib()`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md)
to build one;
[`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md),
which is a coordinate's marginal and the two-category case;
[`dirichlet_distrib()`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md)
for the continuous family on the same simplex, which is conjugate to
this one;
[`parameters7::simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.html)
for the chart.

## Examples

``` r
d <- multinomial_distrib(3, size = 5)
S7::S7_inherits(d, multivariate_distrib)
#> [1] TRUE

# Three categories and five trials give two free probabilities.
c(categories = d@n_dim, trials = d@size, parameters = d@n_params)
#> categories     trials parameters 
#>          3          5          2 
d@params
#> [1] "probs_alr1" "probs_alr2"

# The support is finite and enumerable.
nrow(mv_support(d, NULL))
#> [1] 21
```
