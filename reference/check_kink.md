# Check a Declared Kink Against the Family

Holds a
[`kink_decomposition()`](https://statmodels7.github.io/distributions7/reference/kink_decomposition.md)
to the family that declares it. The declaration says where the kink is,
what sits on top of it and what stands in front; the family's own
density and score are written independently of it, so the two can be
compared.

## Usage

``` r
check_kink(distrib, theta = NULL, tol = 1e-06, verbose = TRUE)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- theta:

  Optional named list of parameter values. Defaults to a probe inside
  every interval.

- tol:

  Relative tolerance for the two numerical comparisons.

- verbose:

  Print the table.

## Value

Invisibly, a named logical vector with one entry per check, `NA` where a
check does not apply.

## Details

Three checks, at probe parameters inside every interval:

- `dv`:

  each \\\partial v/\partial \theta_p\\ against a central difference of
  \\v\\ itself. This isolates `dv`, which nothing else reads.

- `jump`:

  the jump of the score across the kink against the one the declaration
  predicts. Writing \\\Delta\phi'\\ for the jump of \\\phi'\\ at the
  origin – 2 for \\\lvert v \rvert\\ and 1 for \\(v)\_+\\ – the score of
  \\c(\theta)\phi(v)\\ jumps by \$\$c(\theta)\\ \Delta\phi' \\ \partial
  v/\partial \theta_p\$\$ as \\y\\ crosses the surface \\v = 0\\. This
  is the check that matters: it reads `phi`, `v`, `dv` and `coef` at
  once, and a wrong `v` fails it.

- `smooth`:

  the score of every parameter the declaration calls smooth really is
  continuous across the same surface. Without this the check would pass
  a declaration that named too few parameters.

A composition whose \\\Delta\phi'\\ is zero – \\(v)\_+^2\\ and
\\(v)\_+^3\\, where the second and third derivative jump instead – is
reported as not checked rather than checked against zero, which any
declaration would satisfy.

## See also

[`kink_decomposition()`](https://statmodels7.github.io/distributions7/reference/kink_decomposition.md),
[`kink_spec()`](https://statmodels7.github.io/distributions7/reference/kink_spec.md),
[`params_order()`](https://statmodels7.github.io/distributions7/reference/params_order.md),
[`check_distrib()`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
for the family's own battery.

## Examples

``` r
check_kink(laplace_distrib())
#> check_kink: laplace,  phi = abs,  order 0
#>   dv against a difference of v      [PASSED]  worst 4.55e-12
#>   the jump of the score             [PASSED]  worst 0.00e+00
#>   the smooth parameters             [PASSED]
#>       mu         measured            2   declared            2
#>       sigma      measured            0   declared            0
check_kink(laplace2_distrib())
#> check_kink: laplace2,  phi = abs,  order 0
#>   dv against a difference of v      [PASSED]  worst 4.55e-12
#>   the jump of the score             [PASSED]  worst 0.00e+00
#>   the smooth parameters             [PASSED]
#>       mu         measured            2   declared            2
#>       lambda     measured            0   declared            0
```
