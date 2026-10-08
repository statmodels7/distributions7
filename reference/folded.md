# Fold a Distribution at Zero

Wraps a continuous distribution into the distribution of the absolute
value of its variable. The result has density \$\$L(x; \theta) = f(x;
\theta) + f(-x; \theta), \qquad x \ge 0,\$\$ distribution function
\\F(x) - F(-x)\\, and exactly the parent's parameters: folding adds none
and removes none, as
[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md)
does not. The half-normal is
`fixed(folded(gaussian1_distrib()), mu = 0)`, and the folded normal
proper is `folded(gaussian1_distrib())`.

## Usage

``` r
folded(distrib)
```

## Arguments

- distrib:

  A `continuous_distrib` object whose support reaches below zero and
  which declares no atom. A discrete distribution, a distribution
  supported on \\\[0, \infty)\\, an already folded distribution, and a
  zero-adjusted continuous one are each rejected with an error saying
  which condition failed.

## Value

An S7 object of class
[FoldedDistrib](https://statmodels7.github.io/distributions7/reference/FoldedDistrib.md)
carrying `parent_distrib`. Its `params`, `params_bounds` and
`link_params` are the parent's unchanged; `bounds` is
`c(0, max(abs(parent bounds)))`; `distrib_name` is `"folded "` followed
by the parent's; and `params_smooth` is the parent's smoothness.

## A fold is not a change of variable

The map \\y \mapsto \|y\|\\ is TWO TO ONE and has no inverse, so it
cannot be a
[`transformer()`](https://statmodels7.github.io/distributions7/reference/transformer.md).
Instead of carrying a density through a Jacobian it ADDS the two
preimages. Every point above zero has two of them and zero has one,
which is why the constructor is careful about atoms.

## Where the derivatives come from

Every one comes from the wrapper machinery unchanged. Writing \\w =
f(x)/L(x)\\ for the weight of the positive preimage, the ratios that
machinery consumes are \$\$\frac{d^B L}{L} = w \frac{d^B f(x)}{f(x)} +
(1-w) \frac{d^B f(-x)}{f(-x)},\$\$ each term a complete Bell polynomial
in the parent's own log-derivatives, one at \\+x\\ and one at \\-x\\;
\\\log L\\ then follows by the moment-to-cumulant relation. At first
order this reads \\w s(x) + (1-w) s(-x)\\, the score of a two-component
mixture, which is what a fold is. A parent analytic to fourth order
gives a folded family analytic to fourth order.

## What is rejected

A parent that does not reach below zero folds to itself, so the call
would be a no-op; returning it unchanged would hide a mistaken call
rather than report it, and the same check makes `folded()` of a folded
distribution an error. A parent with an atom is rejected too: zero is
its own preimage while every other point has two, so an atom at zero
would be counted twice by the sum above and one elsewhere would be moved
onto its reflection. A discrete parent is rejected outright.

## The sign of a symmetric parent's location is not identified

When the parent is symmetric about its location, \\f(-x; \mu) = f(x;
-\mu)\\, so the two terms of \\L\\ merely swap and the likelihood is
EXACTLY even in \\\mu\\. A fit converges to \\+\hat\mu\\ or \\-\hat\mu\\
according to where it started, at the same maximized value to every
digit: measured on 400 draws at \\\mu = 1.5\\, two fits started at \\\pm
1\\ reach \\\pm 1.50494\\ with the same log-likelihood of \\-511.2516\\.
This is a property of the model, not of the fitting, and it is not
rejected, the folded normal being a standard family; what is estimable
is \\\|\mu\|\\ together with the remaining parameters. Holding the
location at zero removes the question and gives the half-normal. A
parent that is not symmetric about its location, such as
[`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md),
has no such invariance and its sign is identified.

## Notation

\\f\\ is the parent's density, \\F\\ its distribution function, \\L\\
the folded density, \\w\\ the weight of the positive preimage, \\s\\ the
parent's score and \\B\\ a multiset of parameter names.

## See also

[`truncated()`](https://statmodels7.github.io/distributions7/reference/truncated.md),
the other wrapper that adds no parameter,
[`fixed()`](https://statmodels7.github.io/distributions7/reference/fixed.md),
which gives the half-normal,
[`transformation()`](https://statmodels7.github.io/distributions7/reference/transformation.md)
for a map that IS one to one, and
[FoldedDistrib](https://statmodels7.github.io/distributions7/reference/FoldedDistrib.md)
for the class.

## Examples

``` r
d <- folded(gaussian1_distrib())
theta <- list(mu = 0.5, sigma = 1)
distrib_pdf(d, c(0, 0.5, 2), theta)
#> [1] 0.7041307 0.6409130 0.1470459

# The half-normal: a folded gaussian with its location held at zero, whose
# density is twice the gaussian's and whose mean is sigma sqrt(2 / pi).
hn <- fixed(folded(gaussian1_distrib()), mu = 0)
hn@params
#> [1] "sigma"
all.equal(distrib_pdf(hn, c(0.5, 1), list(sigma = 2)),
          2 * dnorm(c(0.5, 1), 0, 2))
#> [1] TRUE
c(mean = mean(hn, list(sigma = 2)), theory = 2 * sqrt(2 / pi))
#>     mean   theory 
#> 1.595769 1.595769 

# The sign of a symmetric parent's location is not identified: two fits
# started either side reach the same maximum at opposite signs.
set.seed(3)
z <- distrib_rng(d, 400, list(mu = 1.5, sigma = 1))
f1 <- fit_distrib(d, z, start = list(mu = 1, sigma = 1))
f2 <- fit_distrib(d, z, start = list(mu = -1, sigma = 1))
rbind(from_plus = c(coef(f1), logLik = as.numeric(logLik(f1))),
      from_minus = c(coef(f2), logLik = as.numeric(logLik(f2))))
#>                   mu    sigma    logLik
#> from_plus   1.504939 1.041631 -511.2516
#> from_minus -1.504939 1.041631 -511.2516

# A parent that is not symmetric about its location has no such invariance.
sn <- folded(skewnormal1_distrib())
set.seed(4)
zs <- distrib_rng(sn, 300, list(mu = 1, sigma = 1, alpha = 3))
lf <- function(m)
  sum(distrib_pdf(sn, zs, list(mu = m, sigma = 1, alpha = 3), log = TRUE))
c(at_plus_1 = lf(1), at_minus_1 = lf(-1))
#>  at_plus_1 at_minus_1 
#>  -266.6752 -1044.9657 

# Three rejections, each naming the condition that failed.
try(folded(gamma1_distrib()))
#> Error : 'gamma1' is supported on [0, Inf], which the absolute value leaves alone.
#>   folded() would return the same distribution, so the call is a mistake
#>   rather than a no-op and is reported as one.
try(folded(folded(gaussian1_distrib())))
#> Error : 'folded gaussian1' is supported on [0, Inf], which the absolute value leaves alone.
#>   folded() would return the same distribution, so the call is a mistake
#>   rather than a no-op and is reported as one.
try(folded(poisson_distrib()))
#> Error : folded() supports continuous distributions only.
```
