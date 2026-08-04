# Fisher Scoring, With Its Own Settings

Returns a specification of Fisher scoring for
[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
carrying how the expected information is to be obtained when the
distribution has no closed form for it.

## Usage

``` r
fisher_scoring(
  approx = c("bartlett", "integrate", "mc", "opg"),
  nsim = 10000,
  criterion = NULL,
  maxit = NULL
)
```

## Arguments

- approx:

  How the expectation is approximated when the distribution has no
  closed-form expected information: `"bartlett"` (the default, the outer
  product of the score, equivalently `"opg"`), `"integrate"` for
  quadrature of the observed information, or `"mc"` for Monte Carlo. See
  [`expected_derivative_methods`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

- nsim:

  Monte Carlo sample size, used when `approx = "mc"`.

- criterion:

  A stopping rule from optimizers7, or `NULL` to use the one
  [`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
  builds from its `tol`.

- maxit:

  An iteration limit, or `NULL` to use the fit's.

## Value

An object of class
[`FisherScoring`](https://statmodels7.github.io/distributions7/reference/FisherScoring.md).

## Details

[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
takes one argument saying how to optimise, and it takes either an
optimiser of optimizers7 or this. The two are the same kind of thing
said in the same place:

|  |  |
|----|----|
| `method = fisher_scoring()` | Newton's method with the **expected** information |
| `method = optimizers7::newton()` | Newton's method with the **observed** Hessian |
| `method = optimizers7::lbfgs()` | whatever that optimiser does |

Fisher scoring is not a separate algorithm, which is why it has no
implementation of its own: it is a Newton step with one matrix replaced
by another. What it does need, and an optimiser cannot carry, is a
statement of how that matrix is to be obtained when the family does not
supply it in closed form — and that is what this object holds. A family
that does supply one ignores `approx` entirely, and
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
refuses the argument in that case rather than accepting something it
will not use.

`criterion` and `maxit` default to `NULL`, meaning the values
[`fit_distrib()`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md)
was given through its own `tol` and `maxit`. Setting them here overrides
those, in the way an optimiser object's own settings already do.

## See also

[`fit_distrib`](https://statmodels7.github.io/distributions7/reference/fit_distrib.md),
[`expected_derivative_methods`](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)

## Examples

``` r
set.seed(1)
d <- gaussian_distrib()
y <- distrib_rng(d, 200, list(mu = 1, sigma = 2))

# the default, and the same thing said explicitly
coef(fit_distrib(d, y))
#>       mu    sigma 
#> 1.071079 1.853543 
coef(fit_distrib(d, y, method = fisher_scoring()))
#>       mu    sigma 
#> 1.071079 1.853543 

# A family whose expected information has no closed form takes a strategy.
# The same argument on a family that HAS one is refused rather than
# silently ignored.
sn <- skewnormal_distrib()
set.seed(2)
ys <- distrib_rng(sn, 300, list(mu = 0, sigma = 1, alpha = 3))
coef(fit_distrib(sn, ys, method = fisher_scoring(approx = "opg")))
#>         mu      sigma      alpha 
#> 0.09063071 1.01986793 3.15712508 

try(fit_distrib(d, y, method = fisher_scoring(approx = "mc")))
#> Error : 'gaussian' computes its expected information in closed form, so the 'approx'
#>   of fisher_scoring() would be ignored. Use fisher_scoring() with no
#>   arguments: the fit will take the exact expression.
```
