#ifndef _counMatrixFactor_COUNTMATRIXFACTOR_H
#define _counMatrixFactor_COUNTMATRIXFACTOR_H

#include <RcppEigen.h>
#include <cstdio>
#include <string>
#include <vector>
#include <math.h>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>

using namespace Rcpp;

using boost::math::digamma;
using boost::math::trigamma;
using boost::math::binomial_coefficient;
using boost::math::factorial;

using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::MatrixXi;                  // variable size matrix, integer
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::VectorXi;                  // variable size vector, double precision
using Eigen::PartialPivLU;              // for fast matrix inversion

#endif