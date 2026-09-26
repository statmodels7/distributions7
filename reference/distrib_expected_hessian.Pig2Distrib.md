# Poisson-Inverse Gaussian Expected Information

The expected information of the Poisson-inverse Gaussian and its first
and second derivatives in the parameters, computed exactly by one pass
over the support for each observation.

## Arguments

- distrib:

  A `Pig1Distrib` or `Pig2Distrib` object.

- y:

  A numeric vector, read for its length.

- theta:

  A named list with `mu` and `sigma` (pig1) or `alpha` (pig2).

- scale:

  `"parameter"` or `"link"`.

- approx, nsim:

  Unused.

- ...:

  Unused.

- threads:

  A single positive integer, the kernel's thread count.

## Value

A named list keyed as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
or
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md).

## Details

For
[`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md)
the expectations are sums over the support of the mass times products of
the observed derivatives, written through the second Bartlett identity
with the measure moving: \$\$\mathbb{E}\[\ell\_{ab}\] =
-\mathbb{E}\[\ell_a\ell_b\],\$\$
\$\$\partial_c\\\mathbb{E}\[\ell\_{ab}\] =
-\mathbb{E}\[\ell\_{ac}\ell_b + \ell_a\ell\_{bc} +
\ell_a\ell_b\ell_c\],\$\$ and the corresponding expression at the second
order. The forms in \\\mathbb{E}\[\ell\_{ab}\]\\ itself sum terms whose
mean is of a higher order in \\1/\alpha\\ than the terms, and lose their
digits in the Poisson limit; the forms in the score do not.

Each observed derivative reads \\\log S_y(\alpha)\\ and the moments of
the finite sum \\S_y\\. Computed afresh at each \\y\\ they cost \\y\\
terms apiece; here they come from the Bessel recurrence \\S\_{y+1} =
S\_{y-1} + \\(2y-1)/\alpha\\S_y\\, differentiated in \\\alpha\\ and
written on positive quantities, so a whole pass costs a number of terms
linear in the length of the support. The pass stops past the mean once
the geometric bound on the remaining mass, the mass at \\y\\ times
\\2\sigma\mu\\, is below \\10^{-17}\\ of the mass accumulated.

For
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
the recurrence runs in \\w = 1/(2\alpha) = \sigma/(2\sqrt{1 +
2\sigma\mu})\\ instead, \\S\_{y+1} = S\_{y-1} + 2(2y-1)\\w\\S_y\\,
because the row of that parametrization reads \\\log S_y\\ and its
derivatives in \\w\\: all of them are finite as \\\sigma \to 0\\, where
those in \\\alpha\\ carry \\\sigma^{-k}\\ and cancel.

The `approx` and `nsim` arguments are accepted and ignored.

## See also

[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md);
`pig1_expected_cpp()` and `pig2_expected_cpp()` for the kernels.

## Examples

``` r
d <- pig2_distrib()
th <- list(mu = 3, alpha = 2)
distrib_expected_hessian(d, 1, th)
#> $mu_mu
#> [1] -0.05598324
#> 
#> $alpha_alpha
#> [1] -0.1167547
#> 
#> $mu_alpha
#> [1] 9.013046e-17
#> 

# Against the outer product of the scores summed over the support.
y <- 0:400
f <- distrib_pdf(d, y, th)
g <- distrib_gradient(d, y, th)
c(kernel = distrib_expected_hessian(d, 1, th)$alpha_alpha,
  sum = -sum(f * g$alpha^2))
#>     kernel        sum 
#> -0.1167547 -0.1167547 
```
