# Expected Hessian of a Reparametrized Distribution

The same congruence applied to the parent's expected information,
\\\mathbb{E}\[\ell^{(ab)}\] = \sum\_{ij}\mathbb{E}\[\ell^{(ij)}\]h^i_a
h^j_b\\. The second-derivative term of the map drops out because the
score has mean zero, so a parent with an exact expected information
gives an exact one here.

## Usage

``` r
reparam_expected_hessian(
  distrib,
  y,
  theta,
  scale = c("parameter", "link"),
  approx = c("bartlett", "integrate", "mc", "opg"),
  nsim = 10000,
  ...
)
```

## Arguments

- distrib:

  A reparametrized distribution.

- y:

  The response.

- theta:

  A named list of the new parameters.

- scale:

  Either `"parameter"` or `"link"`; handled by the generic.

- approx:

  Passed to the parent when it has no closed form.

- nsim:

  Passed to the parent.

- ...:

  Unused.

## Value

A named list of expected second derivatives.

## See also

[`reparametrize`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
