# include <cppad/cppad.hpp>
# include <cppad/core/graph/cpp_graph.hpp>
#include "scorecompdir_types.h"
#include <Rcpp.h>

// [[Rcpp::export]]
void printgraph(Rcpp::XPtr< CppAD::ADFun<double> > pfun){
  CppAD::cpp_graph graph_obj;
  pfun->to_graph(graph_obj);
  graph_obj.print(Rcpp::Rcout);
  return;
};

/* this might fail because 
Currently ``cppad_lib`` library is only needed if one uses
:ref:`colpack<colpack_prefix-name>` ,
:ref:`json_ad_graph-name` ,
:ref:`cpp_ad_graph-name` , or
:ref:`code_gen_fun-name` .. 
AND I'VE REMOVED CPPAD_LIB! */


