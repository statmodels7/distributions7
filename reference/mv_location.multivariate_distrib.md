# Location of a Family Without a Method

Signals an error. The first \\p\\ parameters of a family are not a
location in general, so the class has no default. The shipped families
register their own methods: the elliptical families
[`mv_leading_location()`](https://statmodels7.github.io/distributions7/reference/mv_leading_location.md),
the Dirichlet and the multinomial their mean vector. A family written
elsewhere registers its own method.

## Arguments

- distrib:

  A
  [`multivariate_distrib()`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object with no method of its own.

- theta:

  A named list of parameters. Not read: the error is signalled first.

## Value

Does not return: it signals an error naming the family.

## See also

[`mv_location()`](https://statmodels7.github.io/distributions7/reference/mv_location.md)
for the generic,
[`mv_leading_location()`](https://statmodels7.github.io/distributions7/reference/mv_leading_location.md)
for the implementation that the elliptical families use, and
[`mv_location.MvGaussianDistrib()`](https://statmodels7.github.io/distributions7/reference/mv_location.MvGaussianDistrib.md)
for a method.

## Examples

``` r
# Every shipped family registers a method, so the method of the class is
# called directly.
base <- S7::method(mv_location, multivariate_distrib)
try(base(dirichlet_distrib(3),
         list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8)))
#> Error : mv_location() is not defined for 'dirichlet [3d, mean=simplex]': this family has no location parameter. A family that has one registers a method.

# The Dirichlet's own method returns its mean vector.
mv_location(dirichlet_distrib(3),
            list(mean_alr1 = 0.3, mean_alr2 = -0.2, phi = 8))
#> [1] 0.4260125 0.2583897 0.3155978
```
