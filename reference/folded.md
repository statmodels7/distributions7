# Fold a Distribution at Zero

Wraps a continuous distribution into the distribution of the absolute
value of its variable.

## Usage

``` r
folded(distrib)
```

## Arguments

- distrib:

  A `continuous_distrib` object whose support reaches below zero.

## Value

An S7 object of class
[`FoldedDistrib`](https://statmodels7.github.io/distributions7/reference/FoldedDistrib.md).

## Details

Folding is not a change of variable. The map \\y \mapsto \|y\|\\ is two
to one and has no inverse, so it cannot be a
[`transformer`](https://statmodels7.github.io/distributions7/reference/transformer.md):
instead of carrying a density through a Jacobian it adds the two
preimages, \$\$L(x; \theta) = f(x; \theta) + f(-x; \theta), \qquad x \ge
0,\$\$ with distribution function \\F(x) - F(-x)\\. No parameter is
added and none is removed, as with
[`truncated`](https://statmodels7.github.io/distributions7/reference/truncated.md).

Every derivative comes from the wrapper machinery unchanged. Writing \\w
= f(x)/L(x)\\ for the weight of the positive preimage, the ratios that
machinery consumes are \$\$\dfrac{d^B L}{L} = w \dfrac{d^B f(x)}{f(x)} +
(1-w) \dfrac{d^B f(-x)}{f(-x)},\$\$ each term a complete Bell polynomial
in the parent's own log-derivatives, one evaluated at \\+x\\ and one at
\\-x\\; \\\log L\\ then follows from the moment-to-cumulant relation. At
first order this is \\w s(x) + (1-w) s(-x)\\, the score of a
two-component mixture, which is what a fold is.

The parent must reach below zero. A distribution already supported on
the non-negative half line folds to itself, and returning it unchanged
would hide a mistaken call rather than report it; the same check makes
`folded()` of a folded distribution an error.

A parent with an atom is rejected as well. The point zero is its own
preimage while every other point has two, so an atom at zero would be
counted twice by the sum above and an atom elsewhere would be moved onto
its reflection.

**The half-normal** is `fixed(folded(gaussian1_distrib()), mu = 0)`, and
the folded normal proper is `folded(gaussian1_distrib())`.

**The sign of a symmetric parent's location is not identified.** When
the parent is symmetric about its location, \\f(-x; \mu) = f(x; -\mu)\\,
so the two terms of \\L\\ merely swap and the likelihood is *exactly*
even in \\\mu\\: a fit converges to \\+\hat\mu\\ or \\-\hat\mu\\
according to where it started, at the same maximized value to every
digit. This is a property of the model rather than of the fitting, and
it is not rejected, the folded normal being a standard family; what is
estimable is \\\|\mu\|\\ together with the remaining parameters. Holding
the location at zero removes the question and gives the half-normal. A
parent that is not symmetric about its location, such as
[`skewnormal1_distrib`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md),
has no such invariance and its sign is identified.

## See also

[`truncated`](https://statmodels7.github.io/distributions7/reference/truncated.md),
[`fixed`](https://statmodels7.github.io/distributions7/reference/fixed.md),
[`transformation`](https://statmodels7.github.io/distributions7/reference/transformation.md)

## Examples

``` r
d <- folded(gaussian1_distrib())
theta <- list(mu = 0.5, sigma = 1)
distrib_pdf(d, c(0, 0.5, 2), theta)
#> [1] 0.7041307 0.6409130 0.1470459

# the half-normal: a folded gaussian with its location held at zero
hn <- fixed(folded(gaussian1_distrib()), mu = 0)
hn@params
#> [1] "sigma"
distrib_pdf(hn, c(0.5, 1), list(sigma = 2))
#> [1] 0.3866681 0.3520653
```
