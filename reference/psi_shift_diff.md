# Differences of a Polygamma at an Integer Shift

Returns \\\psi^{(n)}(x+k) - \psi^{(n)}(x)\\ without the cancellation the
direct difference carries. It is the R twin of the `psi_A_rest()` and
`psi_T_rest()` of `src/psi_diff.h`, written for any order rather than
for the two the compiled kernels need.

## Usage

``` r
psi_shift_diff(n, k, x)
```

## Arguments

- n:

  The order of the polygamma: `0` for `digamma`, `1` for `trigamma`, and
  so on.

- k:

  The shift, which for every caller here is a count or a size and
  therefore a non-negative integer. Recycled against `x`.

- x:

  The argument, strictly positive. Recycled against `k`.

## Value

A numeric vector of the same length as the recycled arguments. Exactly
zero wherever `k` is zero, by either branch.

## Details

Several families here carry a boundary limit at which a shape runs to
infinity and the family tends to a simpler one. Every derivative in that
parameter vanishes there, so it is written as a difference of terms that
agree to leading order: both \\\psi^{(n)}(x+k)\\ and \\\psi^{(n)}(x)\\
are of size \\x^{-n}\\ while their difference is of size \\k x^{-n-1}\\,
and the direct subtraction therefore loses a digit for every factor of
ten in \\x/k\\. Consumers that divide the result by a power of the
parameter amplify what is left.

Above \\x = 100\\ the value comes from the asymptotic expansion
\$\$\psi^{(n)}(z) \sim (-1)^{n-1}\Big\[\frac{(n-1)!}{z^{n}} +
\frac{n!}{2 z^{n+1}} + \frac{(n+1)!}{12 z^{n+2}} - \frac{(n+3)!}{720
z^{n+4}} + \frac{(n+5)!}{30240 z^{n+6}}\Big\],\$\$ whose \\n = 0\\ case
is \\\log z - 1/(2z) - 1/(12 z^{2}) + \ldots\\ with the logarithm
differenced as \\\log(1 + k/x)\\. Each power is differenced as
\$\$\frac{1}{(x+k)^{p}} - \frac{1}{x^{p}} =
\frac{1}{x^{p}}\\\big(e^{-p\log(1+k/x)} - 1\big),\$\$ which `expm1` and
`log1p` evaluate to the last bit however small \\k/x\\ is. Below the
crossover the direct difference has all its digits and is used as it
stands.
