# Finite-Difference Steps That Respect the Support

Returns the step for a central difference in the RESPONSE: `h_rel`
scaled by \\\max(1, \|y\|)\\, then shrunk so that \\y \pm h\\ stays
strictly inside the distribution's support. It is the response
counterpart of
[`fd_steps()`](https://statmodels7.github.io/distributions7/reference/fd_steps.md),
and both numerical response derivatives take their step from it.

## Usage

``` r
fd_steps_y(y, bounds, h_rel)
```

## Arguments

- y:

  A numeric vector of observations.

- bounds:

  A numeric vector of length two, the distribution's support. An
  infinite endpoint imposes no clamp on that side.

- h_rel:

  The relative step size, a single positive number. The callers pass
  \\\varepsilon^{1/3}\\ at first order and \\\varepsilon^{1/4}\\ at
  second.

## Value

A numeric vector of steps, as long as `y`, every entry positive.

## Details

The scaling by \\\max(1, \|y\|)\\ makes the step relative where the
response is large and absolute where it is small, so a value near zero
is not differenced with a step below the resolution of a double.

The clamp is what the support requires. A gamma observation at \\y =
10^{-3}\\ differenced with the default step of \\6 \times 10^{-6}\\
needs no help, but one at \\y = 10^{-8}\\ would be evaluated at a
negative point, where the density is not defined and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
returns `-Inf`. The factor 0.49 leaves the step under half the distance
to the bound, so both evaluation points stay inside with room to spare:
at \\y = 10^{-8}\\ the step becomes \\4.9\times 10^{-9}\\ and the left
point \\5.1\times 10^{-9}\\.

## Notation

\\\ell\\ is the log-density of one observation, \\y\\ the response,
\\h\\ the finite-difference step and \\\varepsilon\\ the machine
epsilon, `.Machine$double.eps`.

## See also

[`fd_steps()`](https://statmodels7.github.io/distributions7/reference/fd_steps.md)
for the parameter counterpart, and
[`numerical_grad_y()`](https://statmodels7.github.io/distributions7/reference/numerical_grad_y.md)
and
[`numerical_hess_y()`](https://statmodels7.github.io/distributions7/reference/numerical_hess_y.md),
its two callers.

## Examples

``` r
# On the whole line the step is h_rel times max(1, |y|).
h_rel <- .Machine$double.eps^(1 / 3)
distributions7:::fd_steps_y(c(-100, 0, 0.5, 100), c(-Inf, Inf), h_rel)
#> [1] 6.055454e-04 6.055454e-06 6.055454e-06 6.055454e-04
h_rel * pmax(1, abs(c(-100, 0, 0.5, 100)))
#> [1] 6.055454e-04 6.055454e-06 6.055454e-06 6.055454e-04

# Near a bound the step is cut to 0.49 of the distance to it, so the left
# evaluation point stays inside the support.
y <- c(1e-8, 1e-3, 1)
h <- distributions7:::fd_steps_y(y, c(0, Inf), h_rel)
rbind(step = h, left_point = y - h)
#>               [,1]         [,2]         [,3]
#> step       4.9e-09 6.055454e-06 6.055454e-06
#> left_point 5.1e-09 9.939445e-04 9.999939e-01

# The two steps the callers use differ by a factor of twenty.
c(first_order = .Machine$double.eps^(1 / 3),
  second_order = .Machine$double.eps^(1 / 4))
#>  first_order second_order 
#> 6.055454e-06 1.220703e-04 
```
