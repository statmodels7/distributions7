# Apply a Variable Transformation to a Distribution Object

Builds the distribution of \\Y = g(X)\\, where \\X\\ follows an existing
CONTINUOUS distribution and \\g\\ is a bijective transformation
described by a
[`transformer()`](https://statmodels7.github.io/distributions7/reference/transformer.md).
The result carries exactly the parent's parameters: the transformation
adds none and removes none, and the density is the change-of-variables
formula \$\$f_Y(y) = f_X(g^{-1}(y)) \cdot \left\lvert\frac{d}{dy}
g^{-1}(y)\right\rvert.\$\$

## Usage

``` r
transformation(distrib, transformer, new_name = NULL)
```

## Arguments

- distrib:

  An object inheriting from `continuous_distrib`. A discrete one is
  rejected.

- transformer:

  A
  [`transformer()`](https://statmodels7.github.io/distributions7/reference/transformer.md)
  object, such as
  [`log_transform()`](https://statmodels7.github.io/distributions7/reference/log_transform.md),
  `bc_transform(0.5)` or `affine_transform(1, 2)`. It is rejected where
  its own `valid_support` returns `FALSE` for `distrib@bounds`.

- new_name:

  A single non-empty string naming the result, or `NULL`, the default,
  to compose the parent's name with the transformer's as `"g(parent)"`.

## Value

An S7 object of class
[TransformedDistrib](https://statmodels7.github.io/distributions7/reference/TransformedDistrib.md),
carrying `parent_distrib` and `transformer`. Its `params`,
`params_bounds` and `link_params` are the parent's unchanged; `bounds`
is the transformer's image of the parent's; `params_interpretation`
names the parent's scale in brackets; and `params_smooth` is the
parent's, because a kink in the parameters survives a change of
variables in the response untouched.

## What comes from where

The distribution function, the quantile function and the generator are
obtained by mapping through \\g\\, and for a decreasing transformation
the tails are swapped. Since \\g\\ does not depend on the parameters,
the score, the observed Hessian and the expected Hessian COINCIDE with
the parent's, evaluated at \\x = g^{-1}(y)\\. The function does not
recompute them and loses no accuracy. The moments are available
numerically through
[`moment()`](https://statmodels7.github.io/distributions7/reference/moment.md).

## What is rejected

The function rejects a parent that is not continuous (discrete or
multivariate), since a change of variables of one coordinate needs a
continuous density; a `transformer` argument that is not a `transformer`
object; a transformer that is not valid on the parent's support, as its
`valid_support` function determines; and a `new_name` that is not a
single non-empty string. A map that is not injective, such as the
absolute value, cannot be a transformer at all and is
[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md)
instead.

## Naming the result

Several standard families are a transformation of another family that
the package already provides, and `new_name` sets the name that the
result carries, in place of the default name that records the recipe
that produced it. The reciprocal of a gamma is an inverse gamma, the
exponential of a logistic is a log-logistic, and the exponential of an
exponential is a Pareto. Only the printed name changes, and the
distribution does not depend on it.

## Notation

\\g\\ is the transformation, \\f_X\\ the parent's density, \\f_Y\\ the
transformed one and \\\theta\\ the parameters shared by both.

## See also

[`transformer()`](https://statmodels7.github.io/distributions7/reference/transformer.md)
for the map and the twelve ready-made ones,
[TransformedDistrib](https://statmodels7.github.io/distributions7/reference/TransformedDistrib.md)
for the class,
[`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md)
for a map that is not injective, and
[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md),
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md),
[`zero_inflated()`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md)
and
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
for the other wrappers.

## Examples

``` r
# A lognormal built by transformation is equal to lognormal1_distrib() with
# its second parameter read as a variance.
logn <- transformation(gaussian1_distrib(), exp_transform())
logn@params
#> [1] "mu"    "sigma"
c(built = distrib_pdf(logn, 2, list(mu = 0, sigma = 1)), dlnorm = dlnorm(2))
#>    built   dlnorm 
#> 0.156874 0.156874 

# The reciprocal of a gamma is an inverse gamma, and new_name gives the
# result that name.
ig <- transformation(gamma2_distrib(), inverse_transform(),
                     new_name = "inverse gamma")
ig@distrib_name
#> [1] "inverse gamma"

# The derivatives are the parent's at the preimage, so a fit of the
# transformed family is the parent's fit of the transformed data.
set.seed(1)
y <- distrib_rng(logn, 500, list(mu = 0, sigma = 1))
rbind(transformed = coef(fit_distrib(logn, y)),
      parent_on_logs = coef(fit_distrib(gaussian1_distrib(), log(y))))
#>                        mu    sigma
#> transformed    0.02264409 1.010916
#> parent_on_logs 0.02264409 1.010916

# Three errors follow, and each message names the condition that failed.
try(transformation(poisson_distrib(), log_transform()))
#> Error : transformation() currently supports only continuous distributions.
try(transformation(gaussian1_distrib(), log_transform()))
#> Error : The 'log' transformation is not valid for the support of the 'gaussian1' distribution.
try(transformation(gaussian1_distrib(), exp_transform(), new_name = ""))
#> Error : 'new_name' must be NULL or a single non-empty string.
```
