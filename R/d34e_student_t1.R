#' @include student_t1_distrib.R dexpected_reparam.R d34expected_analytic.R
NULL

# The third and fourth derivatives of student_t1's expected information,
# compiled in src/d34e_student_t1.cpp.
register_d34expected(StudentT1Distrib, function(d, y, th, k, t)
  (if (k == 3L) student_t1_dexpected3_cpp(y, th[[1]], th[[2]], th[[3]], t)
    else student_t1_dexpected4_cpp(y, th[[1]], th[[2]], th[[3]], t)))
