# Part 2 of 3 of the continuous families truncated (helper-ccallable.R).

test_that("every continuous family truncated has entries equal to its methods (part 2 of 3)", {
  cont <- test_part(Filter(function(cls) S7::S7_inherits(ccallable_families[[cls]](),
                                                         continuous_distrib),
                           names(ccallable_families)), 2, 3)
  for (cls in cont) {
    d <- ccallable_families[[cls]]()
    pts <- ccallable_cut(d)
    for (tr in list(truncated(d, lower = pts[1]), truncated(d, upper = pts[2]))) {
      expect_false(is.null(distrib_scalar_route(tr)), label = cls)
      ccallable_twin(NULL, c(-0.5, 0.5), 7, d = tr)
    }
  }
})
