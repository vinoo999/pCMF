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

#ifndef explainedVariance_H
#define explainedVariance_H

/*!
* \file explainedVariance.h
* \brief class definition for explained variance
* \author Ghislain Durif
* \version 0.1
* \date 22/04/2016
*/

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;              // variable size matrix, double precision
using Eigen::MatrixXi;              // variable size matrix, integer
using Eigen::VectorXd;              // variable size matrix, double precision

namespace countMatrixFactor {
    /*!
    * \class explained variance
    * \brief class to define explained variance
    */
    class explainedVariance {
    protected:
        int m_iterMax;

        VectorXd m_expVar0;           /*!< proportion of explained variance as residuals sum of squares */
        VectorXd m_expVarU;           /*!< proportion of variance explained by columns of U (as the ratio between variance of the projection over total variance) */
        VectorXd m_expVarV;           /*!< proportion of variance explained by columns of V (as the ratio between variance of the projection over total variance) */

    public:
        /*!
        * \brief Constructor
        *
        * Constructor of the class explained variance
        */
        explainedVariance(int iterMax);

        /*!
        * \brief Destructor
        *
        * Destructor of the class explained variance
        */
        ~explainedVariance();

    public:
        // getter
        void getExpVar0(VectorXd &res);
        void getExpVarU(VectorXd &res);
        void getExpVarV(VectorXd &res);

        // member functions: doc in src
        double expVar0();
        double expVarU();
        double expVarV();
    };
}

#endif
