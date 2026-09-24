# The Beta-Binomial's Expected-Information Derivatives in Mean and Dispersion

[`betabinom1_distrib()`](https://statmodels7.github.io/distributions7/reference/betabinom1_distrib.md)
is the shapes \\a = \mu/\sigma\\ and \\b = (1-\mu)/\sigma\\ under
another name, so its derivatives are the shapes' ones, summed exactly
over the support by the compiled kernel, and carried across by
[`dexpected_chain()`](https://statmodels7.github.io/distributions7/reference/dexpected_chain.md)
with the map's partials written out: \\a\_\mu = 1/\sigma\\, \\a\_\sigma
= -\mu/\sigma^2\\, \\a\_{\mu\sigma} = -1/\sigma^2\\, \\a\_{\sigma\sigma}
= 2\mu/\sigma^3\\, \\a\_{\mu\sigma\sigma} = 2/\sigma^3\\,
\\a\_{\sigma\sigma\sigma} = -6\mu/\sigma^4\\, and for \\b\\ the same
with \\\mu\\ replaced by \\1-\mu\\ and every derivative in \\\mu\\
changing sign; every derivative taking \\\mu\\ twice is zero.

## Usage

``` r
betabinom1_dexpected(y, mu, sigma, size, order, threads = 1L)
```

## Arguments

- y:

  The response, read for its length.

- mu, sigma:

  The mean proportion and the dispersion.

- size:

  The number of trials.

- order:

  `1L` or `2L`.

- threads:

  Passed to the kernel.

## Value

A named list on the parameter scale, keyed as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
at order 1 and as
[`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
followed by
[`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md)
at order 2.

## See also

[`support_dexpected()`](https://statmodels7.github.io/distributions7/reference/support_dexpected.md),
[`dexpected_chain()`](https://statmodels7.github.io/distributions7/reference/dexpected_chain.md)
