#' @include student_t2_distrib.R d34expected_analytic.R
NULL

# The third and fourth derivatives of student_t2's expected information,
# compiled in src/d34e_student_t2.cpp.
register_d34expected(StudentT2Distrib, function(d, y, th, k, t)
  (if (k == 3L) student_t2_dexpected3_cpp(y, th[[1]], th[[2]], th[[3]], t)
    else student_t2_dexpected4_cpp(y, th[[1]], th[[2]], th[[3]], t)))
