# The Fourth Derivative of the Expected Information Through a Reparametrization

Computes the fourth derivative of the expected information in the new
coordinates \\\phi\\ from the parent's expected information and its
first four derivatives, and the map's partial derivatives up to the
fifth order.

## Usage

``` r
dexpected_chain4(Pp, Pn, E, d1, d2, d3, d4, maps, n)
```

## Arguments

- Pp, Pn:

  The parent's parameter names and the new ones.

- E, d1, d2, d3, d4:

  Named lists keyed as
  [`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
  [`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md),
  [`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md),
  [`d3expected_names()`](https://statmodels7.github.io/distributions7/reference/d3expected_names.md)
  and
  [`d4expected_names()`](https://statmodels7.github.io/distributions7/reference/d4expected_names.md)
  on the parent's parameters.

- maps:

  One keyed table per parent parameter to the order 5.

- n:

  The number of observations, to which every component is recycled.

## Value

A named list keyed as
[`d4expected_names()`](https://statmodels7.github.io/distributions7/reference/d4expected_names.md)
on the new parameters.

## Details

The construction is that of
[`dexpected_chain3()`](https://statmodels7.github.io/distributions7/reference/dexpected_chain3.md)
with four positions. The assignments of the positions to \\G\_{ij}\\,
\\J\_{ia}\\ and \\J\_{jb}\\ number \\3^4\\, and the derivative of
\\G\_{ij}\\ along \\S = \\c, d, e, f\\\\ is the sum over the fifteen set
partitions of \\S\\, written by block sizes: \$\$\partial\_{S} G =
\sum_k E\_{,k}\theta\_{k,S} + \sum\_{kl} E\_{,kl}\Big(\sum\_{x \in S}
\theta\_{k,S\setminus x}\theta\_{l,x} + \theta\_{k,cd}\theta\_{l,ef} +
\theta\_{k,ce}\theta\_{l,df} + \theta\_{k,cf}\theta\_{l,de}\Big)\$\$
\$\$\quad + \sum\_{klm} E\_{,klm} \sum\_{\\x,y\\ \subset S}
\theta\_{k,xy}\\\theta\_{l,z}\\\theta\_{m,w} + \sum\_{klmn}
E\_{,klmn}\\\theta\_{k,c}\theta\_{l,d}\theta\_{m,e}\theta\_{n,f},\$\$
with \\z, w\\ the two positions outside the pair \\\\x, y\\\\. The
factor \\J\_{ia}\\ is differentiated up to four times, so the map is
read to the fifth order.
