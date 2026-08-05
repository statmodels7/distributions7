# Generalized Gamma Distribution Object

Creates a distribution object for the generalized gamma distribution in
Stacy's form, with a scale \\a\\ and two shapes \\d\\ and \\p\\.

## Usage

``` r
gengamma1_distrib(
  link_a = log_link(),
  link_d = log_link(),
  link_p = log_link()
)
```

## Arguments

- link_a:

  A link function object for \\a\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

- link_d:

  A link function object for \\d\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

- link_p:

  A link function object for \\p\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

## Value

An S7 object of class `GenGamma1Distrib`.

## Details

The flexible family for a positive response, and the one that makes a
choice between the gamma, the Weibull and the lognormal something to
estimate rather than to assume.

**Density:** \$\$f(y) =
\dfrac{p}{a^{d}\\\Gamma(d/p)}\\y^{d-1}e^{-(y/a)^{p}}\$\$

**What it nests**, which Stacy's parametrization is chosen to make
visible:

- \\p = 1\\ is the
  [gamma](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md)
  with shape \\d\\ and scale \\a\\;

- \\d = p\\ is the
  [Weibull](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)
  with shape \\p\\ and scale \\a\\;

- \\d = p = 1\\ is the
  [exponential](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md);

- \\p \to 0\\ with \\d/p\\ held large approaches the
  [lognormal](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md).

The first three are exact and testable; the fourth is a limit and is not
reached at any admissible value.

**Score and information.** Writing \\w = (y/a)^{p}\\, \\L = \log(y/a)\\
and \\k = d/p\\, \$\$\dfrac{\partial\ell}{\partial a} = \dfrac{pw-d}{a},
\qquad \dfrac{\partial\ell}{\partial d} = L - \dfrac{\psi(k)}{p}, \qquad
\dfrac{\partial\ell}{\partial p} = \dfrac{1}{p} +
\dfrac{d\psi(k)}{p^{2}} - wL\$\$ and the expected information is
**closed form**: \\u = w\\ is Gamma with shape \\k\\ and unit rate, so
every expectation the Hessian needs is one of \\\mathbb{E}\[u\] = k\\,
\\\mathbb{E}\[u\log u\] = k\psi(k+1)\\ and \\\mathbb{E}\[u(\log u)^{2}\]
= k\\\psi(k+1)^2 + \psi'(k+1)\\\\. That is the same device the Weibull
and the Gumbel use, where the corresponding variable is standard
exponential.

**Moments:** \\\mathbb{E}\[Y^{r}\] =
a^{r}\Gamma\\(d+r)/p\\/\Gamma(d/p)\\, finite for every \\r \> -d\\.

**Parameter domains:** all three are positive.

The three parameters are **weakly identified together** on small
samples: \\d\\ and \\p\\ enter the density largely through their ratio,
and the profile likelihood in that direction is flat. A fit of all three
wants several hundred observations, and holding one with
[`fixed`](https://statmodels7.github.io/distributions7/reference/fixed.md)
is often the better model.

## References

Stacy, E. W. (1962). A generalization of the gamma distribution. *Annals
of Mathematical Statistics* 33, 1187-1192.

## See also

[`gamma2_distrib`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md),
[`weibull1_distrib`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md),
[`lognormal1_distrib`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md),
[`exponential_distrib`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md)

## Examples

``` r
d <- gengamma1_distrib()
d@params
#> [1] "a" "d" "p"

theta <- list(a = 2, d = 3, p = 1.5)
distrib_pdf(d, c(0.5, 2, 5), theta)
#> [1] 0.04136704 0.27590958 0.08999981

# p = 1 is the gamma with shape d and scale a
max(abs(distrib_pdf(d, c(0.5, 2, 5), list(a = 2, d = 3, p = 1)) -
        dgamma(c(0.5, 2, 5), shape = 3, scale = 2)))
#> [1] 2.775558e-17

# d = p is the Weibull with shape p and scale a
max(abs(distrib_pdf(d, c(0.5, 2, 5), list(a = 2, d = 1.5, p = 1.5)) -
        dweibull(c(0.5, 2, 5), shape = 1.5, scale = 2)))
#> [1] 5.551115e-17
```
