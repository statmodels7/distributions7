# Multinomial Observed Hessian

Computes the distinct second derivatives of the multinomial log-mass
with respect to the simplex's free values, one value per observation, in
closed form. With \\A = \partial p/\partial\eta\\ and \\B\_{\cdot,kl} =
\partial^2 p/\partial\eta_k\partial\eta_l\\, \$\$\ell^{(\eta_k\eta_l)} =
\sum\_{j=1}^{p}\left(\dfrac{y_j}{p_j}B\_{j,kl} -
\dfrac{y_j}{p_j^2}A\_{jk}A\_{jl}\right).\$\$ Both terms carry the data,
so the observed matrix differs from its expectation at every
observation.

## Arguments

- distrib:

  A `MultinomialDistrib` object, from
  [`multinomial_distrib()`](https://statmodels7.github.io/distributions7/reference/multinomial_distrib.md).

- y:

  A numeric matrix with one row per observation and \\p\\ columns, each
  row on the support. A single observation may be given as a plain
  vector.

- theta:

  A named list of the simplex's free values on the parameter scale, each
  of length 1.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of numeric vectors, one per distinct entry of the symmetric
matrix and named as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
names them, each of length `nrow(y)`. For \\p\\ categories there are
\\p(p-1)/2\\ of them.

## Notation

\\\ell\\ is the log-mass of one observation, \\p\\ the probability
vector, \\\eta\\ the free vector of the simplex chart, \\A\\ its
Jacobian and \\B\\ its second-derivative arrays.

## See also

[`distrib_gradient.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.MultinomialDistrib.md)
for the score,
[`distrib_expected_hessian.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.MultinomialDistrib.md)
for the expectation of this quantity,
[`distrib_deriv3.MultinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.MultinomialDistrib.md)
for the order above, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- multinomial_distrib(3, size = 5)
th <- list(probs_alr1 = 0.3, probs_alr2 = -0.2)
set.seed(1)
Y <- distrib_rng(d, 4, th)
h <- distrib_hessian(d, Y, th)
names(h)
#> [1] "probs_alr1_probs_alr1" "probs_alr2_probs_alr2" "probs_alr1_probs_alr2"

# numDeriv on the summed log-mass reproduces the summed matrix.
fn <- function(v)
  sum(distrib_pdf(d, Y, as.list(stats::setNames(v, d@params)), log = TRUE))
H <- numDeriv::hessian(fn, unlist(th))
rbind(numeric = c(H[1, 1], H[2, 2], H[1, 2]),
      closed = vapply(h, sum, numeric(1)))
#>         probs_alr1_probs_alr1 probs_alr2_probs_alr2 probs_alr1_probs_alr2
#> numeric             -4.890517             -3.832489              2.201545
#> closed              -4.890517             -3.832489              2.201545

# The mass-weighted sum over the support is the expected Hessian, exactly.
supp <- mv_support(d, th)
w <- distrib_pdf(d, supp, th)
rbind(summed = vapply(distrib_hessian(d, supp, th),
                      function(v) sum(w * v), numeric(1)),
      expected = vapply(distrib_expected_hessian(d, Y, th),
                        function(v) v[1], numeric(1)))
#>          probs_alr1_probs_alr1 probs_alr2_probs_alr2 probs_alr1_probs_alr2
#> summed               -1.222629            -0.9581222             0.5503861
#> expected             -1.222629            -0.9581222             0.5503861
```
