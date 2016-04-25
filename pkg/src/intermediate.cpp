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
 * \file intermediate.cpp
 * \brief intermediate functions
 * \author Ghislain Durif
 * \version 0.1
 * \date 25/04/2016
 */

#include "intermediate.h"

namespace intermediate {

    /*!
     * \fn Dirac function in zero
     *
     * Indicate if the considered value is null
     *
     * @param x value tested
     * @return 0 or 1 indicating if x==0 or not
     *
     * @TODO verify the 0-1 return
     */
    double dirac(double x) {

        if(x==0) {
            return 1;
        } else {
            return 0;
        }
    }
}
