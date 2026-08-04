# The Quantities a Skew t's Moments Are Built From

Returns \\\delta = \alpha/\sqrt{1+\alpha^2}\\, the constant \\b\_\nu =
\sqrt{\nu/\pi}\\\Gamma\\(\nu-1)/2\\/\Gamma(\nu/2)\\, the mean \\\mu_z =
\delta b\_\nu\\ and the variance \\\sigma_z^2\\ of the standardised
variable.

## Usage

``` r
skewt_moment_pieces(alpha, nu)
```

## Arguments

- alpha:

  The shape parameter.

- nu:

  The degrees of freedom.

## Value

A list with `delta`, `bnu`, `mz` and `vz`.

## Details

\\b\_\nu\\ is finite only for \\\nu \> 1\\ and \\\sigma_z^2\\ only for
\\\nu \> 2\\; each is `NaN` otherwise, and the moments that use it
inherit that.
