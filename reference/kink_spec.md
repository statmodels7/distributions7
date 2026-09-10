# The Composition That Carries a Kink

Records the one non-smooth piece of a log-density, as the composition
\\c(\theta)\\\phi(v(y, \theta))\\ in which \\\phi\\ is not smooth at the
origin and everything else is. A family declares one through
[`kink_decomposition()`](https://statmodels7.github.io/distributions7/reference/kink_decomposition.md);
[`params_order()`](https://statmodels7.github.io/distributions7/reference/params_order.md)
reads the order of differentiability off it.

## Usage

``` r
kink_spec(
  phi = character(0),
  v = function() NULL,
  dv = function() NULL,
  coef = function() NULL
)
```

## Arguments

- phi:

  A single string, one of `"abs"`, `"hinge"`, `"hinge2"`, `"hinge3"` or
  `"step"`.

- v:

  A function of `(y, theta)` returning the argument of `phi`, one value
  per observation.

- dv:

  A function of `(y, theta)` returning a named list with one entry per
  parameter, \\\partial v / \partial \theta_p\\. An entry that is
  identically zero says the kink does not move with that parameter.

- coef:

  A function of `theta` returning the multiplier \\c(\theta)\\ standing
  in front of \\\phi\\.

## Value

An object of class `kink_spec` carrying the four fields.

## Details

The five compositions and the order each one leaves, measured as the
first derivative that jumps across \\v = 0\\:

|            |                     |       |
|------------|---------------------|-------|
| `phi`      | \\\phi(v)\\         | order |
| `"abs"`    | \\\lvert v \rvert\\ | 0     |
| `"hinge"`  | \\(v)\_+\\          | 0     |
| `"hinge2"` | \\(v)\_+^2\\        | 1     |
| `"hinge3"` | \\(v)\_+^3\\        | 2     |
| `"step"`   | \\1\\v \> 0\\\\     | -1    |

The others reduce to the absolute value, \\(v)\_+ = (v + \lvert v
\rvert)/2\\ and \\1\\v\>0\\ = (1 + \mathrm{sign}(v))/2\\, so a smoothing
of \\\lvert \cdot \rvert\\ covers all of them by composition.

\\v\\ is not required to be \\y - \mu\\. A Huber likelihood puts its
kink at \\v = \lvert y - \mu \rvert - k\sigma\\, which moves with the
location, the scale and the cut-off together, and the set of non-smooth
parameters is read off \\\partial v / \partial \theta_p \ne 0\\ rather
than assumed.

## See also

[`kink_decomposition()`](https://statmodels7.github.io/distributions7/reference/kink_decomposition.md)
for the generic a family answers,
[`params_order()`](https://statmodels7.github.io/distributions7/reference/params_order.md)
for the order deduced from it, and
[`check_kink()`](https://statmodels7.github.io/distributions7/reference/check_kink.md)
for the validator.

## Examples

``` r
# the Laplace: l = -log(2 sigma) - |y - mu| / sigma
k <- kink_spec(
  phi  = "abs",
  v    = function(y, theta) y - theta$mu,
  dv   = function(y, theta) list(mu = -1, sigma = 0),
  coef = function(theta) -1 / theta$sigma)
k
#> <kink_spec> c(theta) * phi(v),  phi = abs,  order 0
kink_order(k@phi)
#> [1] 0
```
