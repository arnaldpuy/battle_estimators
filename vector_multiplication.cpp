#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericMatrix mmult( NumericMatrix m , NumericVector v , bool byrow = true ){
  if(byrow);
  if(! m.nrow() == v.size()) stop("Non-conformable arrays") ;
  if( ! byrow );
  if(! m.ncol() == v.size()) stop("Non-conformable arrays") ;
  
  NumericMatrix out(m) ;
  
  if(byrow){
    for (int j = 0; j < m.ncol(); j++) {
      for (int i = 0; i < m.nrow(); i++) {
        out(i,j) = m(i,j) * v[j];
      }
    }
  }
  if(! byrow){
    for (int i = 0; i < m.nrow(); i++) {
      for (int j = 0; j < m.ncol(); j++) {
        out(i,j) = m(i,j) * v[i];
      }
    }
  }
  return out ;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

