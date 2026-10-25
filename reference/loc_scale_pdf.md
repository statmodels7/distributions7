# The Density at the Nodes of the Location-Scale Rule

Returns the density at the nodes of the rule. For every family but the
skew t this is
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md).
For a skew t with \\\nu \le 10\\ at \\\|z\| \> 10\\ the density is
formed as a product, \$\$f(z) = \frac{2}{\sigma}\\ C\\ \nu^{(\nu+1)/2}\\
\|z\|^{-\nu}\\ \|z\|^{-1}\\ (1 + \nu/z^2)^{-(\nu+1)/2}\\
T\_{\nu+1}(w),\$\$ with \\C = \Gamma((\nu+1)/2) /
(\sqrt{\nu\pi}\\\Gamma(\nu/2))\\ and \\w\\ as in
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md).
The power of \\\|z\|\\ is taken directly, so that the relative error of
the density does not grow with \\\|\log f\|\\, which the exponential of
the logarithm does at about the unit roundoff times \\\|\log f\|\\.

## Usage

``` r
loc_scale_pdf(distrib, y, th)
```

## Arguments

- distrib:

  A distribution object of a location-scale family.

- y:

  The nodes, shifted by the center of the rule.

- th:

  The parameters at each node, a named list.

## Value

A numeric vector.
