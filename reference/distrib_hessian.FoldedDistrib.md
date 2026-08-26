# Folded Observed Hessian

Computes the second derivatives of the folded log-density by the
moment-to-cumulant relation applied to the block ratios \$\$\frac{d^B
L}{L} = w \frac{d^B f(x)}{f(x)} + (1-w) \frac{d^B f(-x)}{f(-x)},\$\$
which at second order is the familiar mixture Hessian: the weighted
average of the two components' second-order ratios, less the square of
the weighted average of their first-order ones.

## Arguments

- distrib:

  A `FoldedDistrib` object, from
  [`folded()`](https://statmodels7.github.io/distributions7/reference/folded.md).

- y:

  A numeric vector of observations, non-negative.

- theta:

  A named list of the parent's parameters.

- scale:

  One of `"parameter"` (the default) or `"link"`, handled by the generic
  before dispatch.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors, one per unordered pair of parameters,
keyed as
[`hess_names(distrib@params)`](https://statmodels7.github.io/distributions7/reference/hess_names.md).

## Details

No formula is written out for the fold.
[`fold_ratio()`](https://statmodels7.github.io/distributions7/reference/fold_ratio.md)
supplies the ratios and
[`log_deriv()`](https://statmodels7.github.io/distributions7/reference/log_deriv.md)
turns them into derivatives of \\\log L\\, the same two functions the
third and fourth orders use. Each ratio is a complete Bell polynomial in
the parent's own log-derivatives, so a parent with closed forms gives
the folded family closed forms.

## Notation

\\f\\ is the parent's density, \\L\\ the folded one, \\B\\ a multiset of
parameter names, \\w\\ the weight of the positive preimage and \\\ell\\
the folded log-density.

## See also

[`distrib_gradient.FoldedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.FoldedDistrib.md)
for the first order,
[`distrib_deriv3.FoldedDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.FoldedDistrib.md)
for the third,
[`fold_ratio()`](https://statmodels7.github.io/distributions7/reference/fold_ratio.md)
for the ratios, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- folded(gaussian1_distrib())
theta <- list(mu = 0.5, sigma = 1.2)
set.seed(2)
y <- distrib_rng(d, 30, theta)

H <- distrib_hessian(d, y, theta)
vapply(H, sum, numeric(1))
#>       mu_mu sigma_sigma    mu_sigma 
#>  -0.2363048 -53.5849052 -24.9136187 

# Against a numerical Hessian of the folded log-likelihood.
ll <- function(v) {
  t2 <- as.list(v); names(t2) <- d@params
  sum(distrib_pdf(d, y, t2, log = TRUE))
}
Hn <- numDeriv::hessian(ll, unlist(theta))
ref <- vapply(distributions7:::hess_pairs(d@params),
              function(q) Hn[match(q[1], d@params), match(q[2], d@params)],
              numeric(1))
max(abs(vapply(H, sum, numeric(1)) - ref))
#> [1] 5.593428e-10
```
