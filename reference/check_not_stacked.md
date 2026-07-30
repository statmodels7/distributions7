# Refuse to Stack Two Zero Parameters

Rejects an attempt to wrap a distribution that already models the
probability of a zero, and rejects a parameter name the parent has
already used.

## Usage

``` r
check_not_stacked(distrib, fun, param)
```

## Arguments

- distrib:

  The parent distribution being wrapped.

- fun:

  The calling constructor's name, used in the message.

- param:

  The name of the parameter the wrapper wants to add.

## Value

Invisibly `NULL`; raises an error if either condition fails.

## Details

Two zero parameters cannot both be identified. Zero-truncating a
distribution that already has one removes it from the likelihood
entirely – the factor cancels between the numerator and the truncation
constant, leaving its score identically zero – and mixing a further
point mass in only ever shifts the total mass at zero, which one
parameter already describes.

The distributions this rejects are perfectly well-defined; they are
simply not estimable. That is why the constructor is the place to catch
it: nothing goes wrong at run time. The pmf still sums to one,
[`check_distrib`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
still passes, and the fit still converges – to an arbitrary point of a
flat ridge.

## See also

[`zero_inflated`](https://statmodels7.github.io/distributions7/reference/zero_inflated.md),
[`zero_adjusted`](https://statmodels7.github.io/distributions7/reference/zero_adjusted.md)
