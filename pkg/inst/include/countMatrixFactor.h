#include <RcppEigen.h>
#include <math.h>
#include <iostream>
#include <string>
//#include <typeinfo>
//#include <functional>
#include <boost/math/special_functions/digamma.hpp>

using namespace Rcpp;
using namespace std;

using boost::math::digamma;

using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;                  // variable size matrix, integer
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::VectorXi;                  // variable size vector, double precision
using Eigen::PartialPivLU;              // for fast matrix inversion
//using Eigen::JacobiSVD;                 // svd decomposition
//using Eigen::ComputeThinU;              // svd decomposition: compute U
//using Eigen::ComputeThinV;              // svd decomposition: compute V
