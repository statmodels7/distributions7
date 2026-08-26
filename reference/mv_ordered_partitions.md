# Ordered Partitions of a Set of Positions

Enumerates all ordered partitions of a set of positions into nonempty
blocks: every set partition, in every ordering of its blocks. This is
how the derivative of an inverse distributes its differentiations,
\$\$P_t = \sum (-1)^q\\ P A\_{B_1} P \cdots A\_{B_q} P,\$\$ the blocks
appearing in a product where order matters even though the partition
itself is unordered. The count is the ordered Bell number: 1, 3, 13, 75
at sizes 1 to 4.

## Usage

``` r
mv_ordered_partitions(pos)
```

## Arguments

- pos:

  An integer vector of positions, of length 1 to 4 in practice. Its
  values are carried through unchanged, so it may hold any labels the
  caller indexes with.

## Value

A list of ordered partitions. Each is a list of integer vectors, the
blocks in the order the product takes them, and the blocks of one
partition together hold every element of `pos` exactly once.

## Details

The set partitions are built by inserting the last element into each
block of each partition of the rest, and into a block of its own; the
orderings by inserting the last index at each position of each
permutation of the rest. Both recursions are written out here, the sizes
involved being at most 4.
[`numericals7::set_partitions()`](https://statmodels7.github.io/numericals7/reference/set_partitions.html)
supplies the unordered enumeration the rest of the toolkit uses, and
would still leave the block orderings to do.

## Notation

\\P\\ is a precision matrix, \\A_B\\ the derivative of its inverse in
the free values \\B\\, and \\q\\ the number of blocks of a partition.

## See also

[`mvg_ptensors()`](https://statmodels7.github.io/distributions7/reference/mvg_ptensors.md),
the only consumer, and
[`numericals7::set_partitions()`](https://statmodels7.github.io/numericals7/reference/set_partitions.html)
for the unordered enumeration the toolkit uses elsewhere.

## Examples

``` r
# Two positions give three ordered partitions: one block, then the two
# orderings of two singleton blocks.
distributions7:::mv_ordered_partitions(1:2)
#> [[1]]
#> [[1]][[1]]
#> [1] 1 2
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#> [1] 2
#> 
#> [[2]][[2]]
#> [1] 1
#> 
#> 
#> [[3]]
#> [[3]][[1]]
#> [1] 1
#> 
#> [[3]][[2]]
#> [1] 2
#> 
#> 

# The counts are the ordered Bell numbers.
vapply(1:4, function(n)
  length(distributions7:::mv_ordered_partitions(seq_len(n))), integer(1))
#> [1]  1  3 13 75

# Every partition covers the positions exactly once.
all(vapply(distributions7:::mv_ordered_partitions(1:3),
           function(p) identical(sort(unlist(p)), 1:3), TRUE))
#> [1] TRUE
```
