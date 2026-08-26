# Complete Bell Polynomial in the Parent's Log-Derivatives

Computes \\d^I f / f\\ from the derivatives of \\\log f\\, as the sum
over set partitions \$\$\frac{d^I f}{f} = \sum\_{\pi} \prod\_{B \in \pi}
\ell^{(B)},\$\$ where \\\pi\\ runs over the partitions of the
multi-index \\I\\ and \\\ell^{(B)}\\ is the derivative of the
log-density in the parameters that block names.

## Usage

``` r
bell_f_ratio(idx, ell)
```

## Arguments

- idx:

  A character vector of parameter names, with repetition, naming the
  multi-index \\I\\: `c("mu", "mu", "sigma")` is
  \\\partial^3/\partial\mu^2\partial\sigma\\. Its length is the order.

- ell:

  A function of one block, returning \\\ell^{(B)}\\ for the parameters
  that block names, as a numeric vector. Called once per block of every
  partition, so a caller with an expensive parent memoizes it.

## Value

A numeric vector, the length of whatever `ell` returns: the ratio \\d^I
f / f\\ at each observation.

## Details

This is the Bartlett lemma read backwards. The identity is normally used
to eliminate a derivative; here it builds one, and that is what every
wrapper needs: each wrapper's log-likelihood is the parent's log-density
plus, or in place of, \\\log L\\ for some \\\theta\\-dependent \\L\\,
and the orders above two are assembled from this sum together with
[`log_deriv()`](https://statmodels7.github.io/distributions7/reference/log_deriv.md).

At order 1 it returns \\\ell^{(i)}\\ and at order 2 \\\ell^{(ij)} +
\ell^{(i)}\ell^{(j)}\\, the ordinary relation between the second
derivative of a density and of its logarithm. Those two cases reproduce
the hand-written closed forms exactly, and that agreement is the license
for the orders that have nothing to compare against.

## See also

[`log_deriv()`](https://statmodels7.github.io/distributions7/reference/log_deriv.md),
the companion identity for a logarithm;
[`index_partitions()`](https://statmodels7.github.io/distributions7/reference/index_partitions.md),
which supplies the partitions;
[`numericals7::set_partitions()`](https://statmodels7.github.io/numericals7/reference/set_partitions.html)
for the enumeration itself.
