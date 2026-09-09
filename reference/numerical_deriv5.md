# Fifth-Order Derivatives by One Central Difference

The fifth-order derivative components, obtained by differentiating the
fourth-order ones once in each parameter.

## Usage

``` r
numerical_deriv5(
  distrib,
  y,
  theta,
  scale = c("parameter", "link"),
  accuracy = 2L
)
```

## Arguments

- distrib:

  A distribution object.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters, aligned by the generic.

- scale:

  `"parameter"` or `"link"`; the scale the fourth order is read on and
  the scale the perturbation is applied on.

- accuracy:

  The order of accuracy of the central rule, passed to
  [`numericals7::fd_offsets()`](https://statmodels7.github.io/numericals7/reference/fd_offsets.html)
  and
  [`numericals7::fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.html).
  The default 2 is the three-point rule; 4 costs four evaluations of the
  fourth order per parameter instead of two.

## Value

A named list of fifth-derivative component vectors, each of length
`length(y)`, keyed lexicographically as
[`deriv_names(distrib@params, 5)`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
gives them.

## Details

One central difference of the analytic fourth order, and never a
difference of a difference. That is the rule the rest of this package's
numerical surface obeys and it is worth restating here, because the
sibling
[`numerical_deriv4()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv4.md)
does **not** obey it: it differences the analytic Hessian twice, which
it can afford because a second difference of an exact quantity is still
only two orders removed from one. A fifth order built the same way –
differencing the analytic third three times – would not be.

The component of a non-decreasing five-tuple is the four-tuple that
remains when its last index is dropped, differentiated in the parameter
that index names. The tuple is non-decreasing, so dropping the last
index leaves a valid order-4 key, and mixed partials commute, so which
index is dropped does not matter. Components are grouped by that last
index, so the whole order costs `2p` evaluations of the fourth rather
than two per component.

## The stencil, and the step

The nodes, the weights and the step are all numericals7's, so the rule
lives in one place and raising `accuracy` is an argument rather than new
arithmetic.
[`numericals7::fd_derivative()`](https://statmodels7.github.io/numericals7/reference/fd_derivative.html)
is not called directly for the reason
[`numDeriv_grad()`](https://statmodels7.github.io/distributions7/reference/numDeriv_grad.md)
already records: its `f` maps a vector of points to the values at those
points, while this one reads a whole named list at each node and keeps
every component that shares the differentiated index. A node of zero
weight is skipped, so the default central rule costs two evaluations of
the fourth order per parameter and not three.

[`numericals7::fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.html)
at order 1 is \\\epsilon^{1/3}\\ scaled by the magnitude of the
evaluation point and shrunk to keep the stencil inside the parameter's
domain. It is the right rule here because the quantity being differenced
is evaluated to machine precision. The same rule would be a thousand
times too small in statmodels7's outer stencil, where what limits the
step is the reproducibility of refitting the mode rather than the
rounding of an arithmetic expression. The two cases look alike and take
opposite answers.

On the link scale the perturbation is applied to \\\eta\\ and the fourth
order is asked for on the link scale as well, so the result is the
derivative of an analytic link-scale quantity. Nothing here needs the
fifth derivative of a link, nor an order-5 entry in
[`bell_partial()`](https://statmodels7.github.io/distributions7/reference/bell_partial.md),
both of which the chain rule route would have required.

## See also

[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md),
the generic;
[`numerical_deriv4()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv4.md),
the order below and the one shape not to copy.

## Examples

``` r
numerical_deriv5(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#> $mu_mu_mu_mu_mu
#> [1] 0 0 0
#> 
#> $mu_mu_mu_mu_sigma
#> [1] 0 0 0
#> 
#> $mu_mu_mu_sigma_sigma
#> [1] 0 0 0
#> 
#> $mu_mu_sigma_sigma_sigma
#> [1] 24 24 24
#> 
#> $mu_sigma_sigma_sigma_sigma
#> [1] -120    0  120
#> 
#> $sigma_sigma_sigma_sigma_sigma
#> [1] 336 -24 336
#> 
```
