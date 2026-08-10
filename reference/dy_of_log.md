# The k-th Response Derivative of a Function of the Log Response

Evaluates \\\partial^{k} g(\log y)/\partial y^{k}\\ from the derivatives
of \\g\\, by \\y^{-k}\sum\_{j} s(k, j)\\ g^{(j)}(\log y)\\ with \\s(k,
j)\\ the signed Stirling numbers of the first kind.

## Usage

``` r
dy_of_log(gd, y, k)
```

## Arguments

- gd:

  A list whose \\j\\-th element is \\g^{(j)}(\log y)\\.

- y:

  The response.

- k:

  The derivative order, 1 to 4.

## Value

A numeric vector the length of `y`.
