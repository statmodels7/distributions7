# Beta-Binomial Distribution in Its Shapes

Creates a beta-binomial distribution object in its canonical
parametrization, the two beta shapes.

## Usage

``` r
betabinom2_distrib(size, link_alpha = log_link(), link_beta = log_link())
```

## Arguments

- size:

  The number of trials, a constant of the distribution.

- link_alpha:

  Link function for \\\alpha\\. Defaults to the log.

- link_beta:

  Link function for \\\beta\\. Defaults to the log.

## Value

An S7 object of class
[`BetaBinom2Distrib`](https://statmodels7.github.io/distributions7/reference/BetaBinom2Distrib.md).

## Details

The same law as
[`betabinom1_distrib`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md),
which carries a mean proportion and a dispersion: \\\alpha =
\mu/\sigma\\ and \\\beta = (1-\mu)/\sigma\\.

Every derivative is a sum of polygamma functions, since the log-mass is
a sum of log-gamma terms and differentiating it \\k\\ times replaces
each by \\\psi^{(k-1)}\\. The expectations are **exact finite sums**
over \\\\0, \dots, n\\\\ rather than quadratures.

## See also

[`betabinom1_distrib`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md)

## Examples

``` r
d <- betabinom2_distrib(size = 10)
theta <- list(alpha = 2, beta = 3)
distrib_pdf(d, 0:10, theta)
#>  [1] 0.06593407 0.10989011 0.13486513 0.14385614 0.13986014 0.12587413
#>  [7] 0.10489510 0.07992008 0.05394605 0.02997003 0.01098901
c(mean = mean(d, theta), variance = variance(d, theta))
#>     mean variance 
#>        4        6 
```
