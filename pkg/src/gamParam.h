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

#ifndef gamParam_H
#define gamParam_H

/*!
* \file gamParam.h
* \brief class definition for Gamma parameters
* \author Ghislain Durif
* \version 0.1
* \date 22/04/2016
*/

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;              // variable size matrix, double precision

namespace countMatrixFactor {

    /*!
    * \class gamParam
    * \brief class to define Gamma distribution parameters
    */
    class gamParam {
    protected:
        // dimensions
        int m_rows;         /*!< number of rows (dimension
                            of observation or variable space) */
        int m_cols;         /*!< number of factors */

        // hyper-parameters
        MatrixXd m_param1;  /*!< rows x cols, shape parameter (called alpha) */
        MatrixXd m_param2;  /*!< rows x cols, rate parameter (called beta) */

    public:
        /*!
        * \brief Constructor
        *
        * Constructor of the class gamParam without initialization
        */
        gamParam(int rows, int cols);

        /*!
        * \brief Constructor
        *
        * Constructor of the class gamParam with initialization
        */
        gamParam(int rows, int cols,
                 const MatrixXd &param1, const MatrixXd &param2);

        /*!
        * \brief Destructor
        *
        * Destructor of the class gamDistrib
        */
        ~gamParam();

    public:
        // getter
        void getParam1(MatrixXd &param1);
        void getParam2(MatrixXd &param2);

        // setter
        void setParam1(MatrixXd &param1);
        void setParam2(MatrixXd &param2);

        // member functions: documented in src

        // expectation
        void E(MatrixXd &res);

        // log-expectation
        void Elog(MatrixXd &res);

        // entropy
        void entropy(MatrixXd &res);

        // parameter norm
        double parameterNorm2();
    };

    // functions: documented in src

    // expectation
    void Egam(const MatrixXd &param1, const MatrixXd &param2, MatrixXd &res);

    // log-expectation
    void Elgam(const MatrixXd &param1, const MatrixXd &param2, MatrixXd &res);

    // entropy
    void entropyGam(const MatrixXd &param1, const MatrixXd &param2, MatrixXd &res);

    // parameter norm
    double parameterNorm2(const MatrixXd &param1, const MatrixXd &param2);

}

#endif
