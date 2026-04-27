#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
double series_cpp(Function f, double start, double end, 
                  int step = 1000, double tol = 1e-10, 
                  int maxit = 1000000, bool reltol = true) {
    
    double s = 0.0;
    double start_internal = start;
    double end_internal = end;
    
    int it = 0;
    bool climbing = true;
    double prev_max_abs = R_PosInf;
    int divergence_counter = 0;
    int flat_counter = 0; 
    
    while (climbing && it < maxit) {
        it++;
        
        // Upper limit of the current block (chunk)
        double upper_limit = std::min(start_internal + step - 1.0, end_internal);
        int current_step_size = static_cast<int>(upper_limit - start_internal + 1.0);
        
        NumericVector x(current_step_size);
        for (int i = 0; i < current_step_size; i++) {
            x[i] = start_internal + i;
        }
        
        // Vectorized evaluation of the R function (called from C++)
        NumericVector vals = f(x);
        
        // Internal analysis of the block
        double chunk_sum = 0.0;
        double chunk_abs_sum = 0.0;
        double max_abs_val = 0.0;
        for (int i = 0; i < current_step_size; i++) {
            chunk_sum += vals[i];
            double abs_val = std::abs(vals[i]);
            chunk_abs_sum += abs_val;
            if (abs_val > max_abs_val) {
                max_abs_val = abs_val;
            }
        }
        
        s += chunk_sum;
        
        // 1. Explicit divergence check (NaN or Infinity)
        if (std::isinf(s) || std::isnan(s)) {
            warning("The series reached Inf or NaN. Stopping.");
            return s;
        }
        
        // 2. Divergence Early-Exit (terms are growing in absolute value)
        if (max_abs_val > prev_max_abs && max_abs_val > tol) {
            divergence_counter++;
            if (divergence_counter >= 3) {
                warning("The series seems to be divergent (terms are growing). Stopping early.");
                return s;
            }
        } else {
            divergence_counter = 0;
        }
        prev_max_abs = max_abs_val;
        
        // 3. Convergence or Underflow Protection
        double scaled_tol = reltol ? (tol * std::max(std::abs(s), 1.0)) : tol;
        
        if (chunk_abs_sum < scaled_tol) {
            flat_counter++;
            // Stop if we have accumulated enough mass (|s| > tol) 
            // or if we have scanned up to 50 empty blocks (50 * step = 50,000 zeros in a row)
            if (std::abs(s) > tol || flat_counter >= 50) {
                climbing = false;
            }
        } else {
            flat_counter = 0; 
        }
        
        // 4. Advancement
        if (climbing) {
            start_internal = upper_limit + 1.0;
            if (start_internal > end_internal) {
                climbing = false;
            }
        }
    }
    
    if (it >= maxit) {
        warning("Maximum number of iterations reached.");
    }
    
    return s;
}
