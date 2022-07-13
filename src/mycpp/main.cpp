//R CMD LINK g++ -I/usr/share/R/include -I../inst/include -I/home/kassel/R/x86_64-pc-linux-gnu-library/4.0/Rcpp/include -Wl,--export-dynamic -fopenmp -Wl,-Bsymbolic-functions -Wl,-z,relro -L/usr/lib/R/lib -lR -lpcre2-8 -llzma -lbz2 -lz -lrt -ldl -lm -licuuc -licui18n  mycpp/main.cpp -o klhexec

# include "scm.cpp"

bool print_graph(void)
{   bool ok = true;
    using std::string;
    //
    // AD graph example
    // node_1 : p[0]
    // node_2 : p[1]
    // node_3 : x[0]
    // node_4 : p[0] + p[1]
    // node_5 : x[0] + ( p[0] + p[1] )
    // y[0]   = x[0] + ( p[0] + p[1] )
    //
    // C++ graph object
    CppAD::cpp_graph graph_obj;
    //
    // operator being used
    CppAD::graph::graph_op_enum op_enum;
    //
    // set scalars
    graph_obj.function_name_set("print_graph example");
    size_t n_dynamic_ind = 2;
    graph_obj.n_dynamic_ind_set(n_dynamic_ind);
    size_t n_variable_ind = 1;
    graph_obj.n_variable_ind_set(n_variable_ind);
    //
    // node_4 : p[0] + p[1]
    op_enum = CppAD::graph::add_graph_op;
    graph_obj.operator_vec_push_back(op_enum);
    graph_obj.operator_arg_push_back(1);
    graph_obj.operator_arg_push_back(2);
    //
    // node_5 : x[0] + ( p[0] + p[1] )
    graph_obj.operator_vec_push_back(op_enum);
    graph_obj.operator_arg_push_back(3);
    graph_obj.operator_arg_push_back(4);
    //
    // y[0]   = x[0] + ( p[0] + p[1] )
    graph_obj.dependent_vec_push_back(5);
    //
    // get output of print command
    std::stringstream os;
    graph_obj.print(os);
    //
    std::string check =
        "print_graph example\n"
        "          1      p[0]\n"
        "          2      p[1]\n"
        "          3      x[0]\n"
        "          4       add    1    2\n"
        "          5       add    3    4\n"
        "y nodes = 5\n"
    ;
    std::string str = os.str();
    ok &= str == check;
    //
    return ok;
}


int main(int argc, char** argv)
 {   using CppAD::AD;   // use AD as abbreviation for CppAD::AD
       std::cout << "Reading inputs" << std::endl;
       //read inputs
       size_t n = 3;
       veca1 z_ad(n);
       for(int i=0; i < n; i++){
           z_ad[i] = 0.1;
       }
       veca1 theta_ad(8);
       for(int i=0; i < theta_ad.size(); i++){
           theta_ad[i] = 1. + i/2.;
       }

       Eigen::Matrix<int, Eigen::Dynamic, 1> fixedtheta_e(theta_ad.size());
       for(int i=0; i < fixedtheta_e.size() - 2; i++){
           fixedtheta_e[i] = 0;
       }
       fixedtheta_e[fixedtheta_e.size() - 1] = 1;
       fixedtheta_e[fixedtheta_e.size() - 2] = 1;

       std::cout << "Creating manifold object" << std::endl;
       manifold<a1type> * man;
       man = new Spos<a1type>();

   std::cout << "Preparing to tape" << std::endl;
   CppAD::ADFun<double> out; //returning a pointer
   out = tapell(z_ad,
                 theta_ad,
                 ll_ppi,
                 man,
                 fixedtheta_e,
                 true);

   CppAD::cpp_graph graph_obj;
   std::cout <<"Printing..." << std::endl;
   std::cout <<"Function name is currently: " << graph_obj.function_name_get() << std::endl;
   print_graph();

   out.to_graph(graph_obj);
   std::cout << std::endl << " Graph Created " << std::endl;

   std::stringstream os;
   graph_obj.print(os);

   std::cout << std::endl << " Finished printing graph " << std::endl;

   delete man;

   return 0;
 }
