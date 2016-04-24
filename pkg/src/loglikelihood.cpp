// Copyright 2016-04 Ghislain Durif
//
// This file is part of the `hClustering' library for R and related languages.
// It is made available under the terms of the GNU General Public
// License, version 2, or at your option, any later version,
// incorporated herein by reference.
//
// This program is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public
// License along with this program; if not, write to the Free
// Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
// MA 02111-1307, USA

/*!
* \file loglikelihood.cpp
* \brief class definition for log-likelihood
* \author Ghislain Durif
* \version 0.1
* \date 22/04/2016
*/

#include <Rcpp.h>
#include <RcppEigen.h>
#include "loglikelihood.h"

// [[Rcpp::depends(RcppEigen)]]
using Eigen::VectorXd;              // variable size matrix, double precision

namespace countMatrixFactor {

    //------------------------------------------------------------------------//
    // CONSTRUCTOR
    loglikelihood::loglikelihood(int iterMax) {
        m_iterMax = iterMax;
        m_margLogLike = VectorXd::Zero(iterMax);
        m_condLogLike = VectorXd::Zero(iterMax);
        m_priorLogLike = VectorXd::Zero(iterMax);
        m_postLogLike = VectorXd::Zero(iterMax);
        m_compLogLike = VectorXd::Zero(iterMax);
    }

    // DESTRUCTOR
    loglikelihood::~loglikelihood() {}

    // getter

    void loglikelihood::getMarginal(VectorXd &res) {
        res = m_margLogLike;
    }
    void loglikelihood::getConditional(VectorXd &res) {
        res = m_condLogLike;
    }
    void loglikelihood::getPrior(VectorXd &res) {
        res = m_priorLogLike;
    }
    void loglikelihood::getPosterior(VectorXd &res) {
        res = m_postLogLike;
    }
    void loglikelihood::getComplete(VectorXd &res) {
        res = m_compLogLike;
    }

    // MEMBER FUNCTIONS
    void computeLogLike() {}

}










