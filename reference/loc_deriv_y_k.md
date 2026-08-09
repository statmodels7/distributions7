# The Response Derivative of a Location Family

Builds the order-`k` response derivative of a family whose response
enters only as \\y - \mu\\, as \\(-1)^{k}\\ times the pure derivative in
the location.

## Usage

``` r
loc_deriv_y_k(order)
```

## Arguments

- order:

  The derivative order, 3 or 4.

## Value

A function suitable for registering as an S7 method.

## Details

The identity is exact and needs no formula of its own, which is the
point: the location derivative is already written, often as a compiled
kernel, and the response derivative is the same number with a sign. It
also fixes the two orders below, where \\\partial\ell/\partial y =
-\partial\ell/ \partial\mu\\ and \\\partial^{2}\ell/\partial y^{2} =
\partial^{2}\ell/\partial\mu^{2}\\, and the tests check the new orders
against exactly that.
