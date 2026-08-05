# Multinomial Distribution Object

Creates a distribution object for the multinomial distribution,
parametrized by a probability vector on the simplex.

## Usage

``` r
multinomial_distrib(n_dim, size, probs = parameters7::simplex(n_dim))
```

## Arguments

- n_dim:

  The number of categories \\p\\.

- size:

  The number of trials \\n\\. A constant of the distribution rather than
  a parameter, as for
  [`binomial_distrib`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md).

- probs:

  A parameters7
  [`simplex`](https://statmodels7.github.io/parameters7/reference/simplex.html)
  of the same dimension. Defaults to `parameters7::simplex(n_dim)`.

## Value

An S7 object of class `MultinomialDistrib`.

## Details

The first family here that is multivariate and **discrete**, so its
support is a finite set of points — every vector of non-negative
integers summing to \\n\\ — rather than a region.
[`mv_support`](https://statmodels7.github.io/distributions7/reference/mv_support.md)
enumerates them, which is what lets an expectation be an exact sum and
the validator check the total mass by addition rather than by sampling.

The probabilities are carried by a parameters7 simplex and flattened
into scalars with identity links, exactly as a covariance is for the
multivariate gaussian: the constraint that they be positive and sum to
one lives in the parameter, where a scalar link could not express it.

**Probability mass function:** \$\$P(Y=y) = \dfrac{n!}{\prod_j
y_j!}\prod_j p_j^{y_j}\$\$

**Score and information.** With \\A = \partial p/\partial\eta\\,
\$\$\dfrac{\partial\ell}{\partial\eta_k} = \sum_j
\dfrac{y_j}{p_j}A\_{jk}, \qquad \mathbb{E}\[\ell^{(\eta_k\eta_l)}\] =
-n\sum_j \dfrac{A\_{jk}A\_{jl}}{p_j}\$\$ The expected form is closed
because \\\mathbb{E}\[y_j\] = np_j\\ turns the second-derivative term
into \\n\sum_j B\_{j,kl}\\, which vanishes: the probabilities sum to
one, so every derivative of their sum is zero.

**Moments:** mean \\np\\ and \\\operatorname{Cov}(Y_i,Y_j) =
n(\delta\_{ij}p_i - p_ip_j)\\, singular by construction.

**The marginals are binomial**, coordinate \\j\\ being
\\\mathrm{Binomial}(n, p_j)\\ with the other categories collapsed into a
single failure.

## See also

[`binomial_distrib`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md),
[`dirichlet_distrib`](https://statmodels7.github.io/distributions7/reference/dirichlet_distrib.md)
for the family that is conjugate to it,
[`simplex`](https://statmodels7.github.io/parameters7/reference/simplex.html)

## Examples

``` r
d <- multinomial_distrib(3, size = 5)
d@params
#> [1] "probs_alr1" "probs_alr2"

theta <- as.list(stats::setNames(c(0.3, -0.2), d@params))
mv_location(d, theta)
#> [1] 2.130063 1.291948 1.577989

# the support is a finite set of points, so the mass sums exactly
supp <- mv_support(d, theta)
c(points = nrow(supp), mass = sum(distrib_pdf(d, supp, theta)))
#> points   mass 
#>     21      1 
```
