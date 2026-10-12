# Replace the Derivatives of a Truncated Discrete Part Near Degeneracy

Overwrites, in a result of
[`distrib_gradient()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md),
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
or of the third to fifth derivatives of a discrete
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
family or of a discrete
[`zero_adjusted()`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
family, the components in the parent's parameters at the observations
where the retained mass sits at the first retained point, \\f(L)/Z \>
0.9\\, with the values of
[`head_derivs()`](https://statmodels7.github.io/distributions7/reference/head_derivs.md).
The components in the zero-adjustment parameter are left as they are.

## Usage

``` r
head_override(distrib, y, theta, order, res)
```

## Arguments

- distrib:

  A `TruncatedDiscreteDistrib` or `ZeroAdjustedDiscreteDistrib`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters, aligned.

- order:

  The derivative order of `res`, 1 to 5.

- res:

  The result to correct, a named list of components.

## Value

`res`, with the components replaced where the condition holds.

## See also

[`head_derivs()`](https://statmodels7.github.io/distributions7/reference/head_derivs.md).

## Examples

``` r
d <- truncated(poisson_distrib(), lower = 1)
distrib_hessian(d, 1, list(mu = 1e-6))$mu_mu
#> [1] -0.08333333
```
