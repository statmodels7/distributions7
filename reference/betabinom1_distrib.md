# Beta-Binomial Distribution Object

Creates a distribution object for the beta-binomial distribution,
parametrized by the mean proportion \\\mu\\ and a dispersion parameter
\\\sigma\\.

## Usage

``` r
betabinom1_distrib(size, link_mu = logit_link(), link_sigma = log_link())
```

## Arguments

- size:

  The number of trials \\n\\. A constant of the distribution rather than
  a parameter, as for
  [`binomial_distrib`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md).

- link_mu:

  A link function object for \\\mu\\. Defaults to
  [`logit_link`](https://statmodels7.github.io/linkfunctions7/reference/logit_link.html),
  the mean being a proportion.

- link_sigma:

  A link function object for \\\sigma\\. Defaults to
  [`log_link`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  to ensure positivity.

## Value

An S7 object of class `BetaBinom1Distrib`.

## Details

The family is the binomial with its success probability drawn from a
Beta, which is what makes it the natural model for a proportion whose
trials are not independent. Writing \\\alpha = \mu/\sigma\\ and \\\beta
= (1-\mu)/\sigma\\, the mean is \\n\mu\\ and \$\$\operatorname{Var}(Y) =
n\mu(1-\mu) \left(1 + (n-1)\dfrac{\sigma}{1+\sigma}\right),\$\$ so the
family is overdispersed relative to the binomial at every \\\sigma \>
0\\ and approaches it as \\\sigma \to 0\\. The intraclass correlation is
\\\sigma/(1+\sigma)\\.

**Probability mass function:** \$\$P(Y = y) = \binom{n}{y}
\dfrac{B(y+\alpha,\\ n-y+\beta)}{B(\alpha, \beta)}\$\$

**Score and information.** The parameters enter only through the two
shapes, where every derivative is a difference of polygammas, so the
score is that difference carried through the chain rule of \\(\alpha,
\beta) = (\mu/\sigma, (1-\mu)/\sigma)\\. The expected information is an
**exact sum** over the finite support rather than a quadrature, which is
what a bounded count buys.

**Parameter domains:**

- \\\mu \in (0, 1)\\

- \\\sigma \in (0, +\infty)\\

The family is **not** reachable from anything already in the package: it
is neither a binomial with a parameter held fixed nor a wrapper over
one, the mixing being over the success probability rather than over the
outcome.

Third and fourth derivatives come from the numerical fallback. On a
finite support that fallback differences an exact mass function, so its
accuracy is the usual `1e-8` rather than the poorer figure a quadrature
would give.

## See also

[`binomial_distrib`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md),
[`beta1_distrib`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md),
[`negbin2_distrib`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)

## Examples

``` r
d <- betabinom1_distrib(size = 10)
d@params
#> [1] "mu"    "sigma"

theta <- list(mu = 0.3, sigma = 0.5)
distrib_pdf(d, 0:10, theta)
#>  [1] 0.26446066 0.15257346 0.11686478 0.09645982 0.08212120 0.07082953
#>  [7] 0.06121071 0.05246632 0.04397912 0.03502041 0.02401400
c(mean = mean(d, theta), variance = variance(d, theta))
#>     mean variance 
#>      3.0      8.4 

# overdispersed relative to the binomial with the same mean
variance(binomial_distrib(size = 10), list(mu = 0.3))
#> [1] 2.1
```
