# One Minus the Squared Skewness Parameter

Returns \\1 - \delta^2\\, the quantity the shape of the direct
parametrization divides by. Written in the skewness alone it is \$\$1 -
\delta^2 = 1 -
\left(\dfrac{\lvert\gamma_1\rvert}{\gamma\_{\max}}\right)^{2/3},\$\$
with \\\gamma\_{\max}\\ the ceiling of
[`sn_max_skew()`](https://statmodels7.github.io/distributions7/reference/sn_max_skew.md),
and it is evaluated through `log` and `expm1`, which leaves no
cancellation in it.

## Usage

``` r
sn_one_minus_delta2(gamma1, s)
```

## Arguments

- gamma1:

  The skewness, a numeric vector. Nothing is validated here: the result
  is positive strictly inside \\(-\gamma\_{\max}, \gamma\_{\max})\\,
  zero at either bound and negative outside them.

- s:

  The sign of `gamma1`, \\\pm 1\\, taken by the caller from its plain
  value.

## Value

A numeric vector of the length of `gamma1`.

## Details

The identity is exact rather than an expansion near the ceiling. With
\\c = (2\lvert\gamma_1\rvert/(4-\pi))^{1/3}\\ the map forms \\\mu_z =
c/\sqrt{1+c^2}\\ and \\\delta = \mu_z/b\\, and \\c^2\\ equals
\\(\lvert\gamma_1\rvert/\gamma\_{\max})^{2/3}\\b^2/(1-b^2)\\, so the
ratio to the ceiling carries the whole of the degeneracy.

The two forms this replaces differ from it only at the top of the range,
and both subtract two nearly equal numbers there. `1 - (mu_z/b)^2`
reaches exactly zero at the largest skewness
[`linkfunctions7::bounded_link()`](https://statmodels7.github.io/linkfunctions7/reference/bounded_link.html)
can produce, and the \\b^2 + (b^2-1)c^2\\ of
[`md_skewnormal2()`](https://statmodels7.github.io/distributions7/reference/reparam_map_derivs.md)
reaches \\-1.11\times 10^{-16}\\. At that skewness, one unit in the last
place inside the bound, this returns \\9.42\times 10^{-17}\\, and the
shape that follows from it is \\1.36\times 10^{8}\\ and finite.

## See also

[`sn_cp_to_dp()`](https://statmodels7.github.io/distributions7/reference/sn_cp_to_dp.md)
and
[`md_skewnormal2()`](https://statmodels7.github.io/distributions7/reference/reparam_map_derivs.md),
the two places that divide by it, and
[`sn_max_skew()`](https://statmodels7.github.io/distributions7/reference/sn_max_skew.md)
for the ceiling it is written against.

## Examples

``` r
g <- distributions7:::sn_max_skew() * c(0.5, 0.9, 0.999, 1 - 1e-15)
distributions7:::sn_one_minus_delta2(g, 1)
#> [1] 3.700395e-01 6.783025e-02 6.667778e-04 6.661338e-16

# The expression it replaces has already lost the last of those.
b <- distributions7:::sn_b()
cc <- (2 * g / (4 - pi))^(1 / 3)
1 - (cc / sqrt(1 + cc^2) / b)^2
#> [1] 0.1759033670 0.0257605658 0.0002423968 0.0000000000
```
