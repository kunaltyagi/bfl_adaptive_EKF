// $Id: adaptivekalmanfilter.h  2015-06-05 10:59:21Z kunaltyagi $
// Copyright (C) 2015 Kunal Tyagi <last dot first at live dot com>
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//

// Using ADAPTIVE KALMAN FILTER FOR NOISE IDENTIFICATION by M. Oussalah and J. De Schutter
// available at https://www.isma-isaac.be/publications/PMA_MOD_publications/ISMA25/p1225p1232.pdf

#ifndef __ADAPTIVE_EXTENDED_KALMAN_FILTER__
#define __ADAPTIVE_EXTENDED_KALMAN_FILTER__

#include "extendedkalmanfilter.h"
#include "../pdf/conditionalpdf.h"
#include "../pdf/gaussian.h"
# include <map>

namespace BFL
{
class AdaptiveExtendedKalmanFilter : public ExtendedKalmanFilter  // maybe Filter<MatrixWrapper::ColumnVector,MatrixWrapper::ColumnVector> because of private members in KalmanFilter.h
{
public:
  AdaptiveExtendedKalmanFilter(Gaussian* prior);
  
  virtual ~AdaptiveExtendedKalmanFilter();
  
  //  For realtime use, this function should be called before calling measUpdate
  /*  @param vector containing the dimension of the measurement models which are
      going to be used
  */
  void AllocateMeasModelExt( const vector<unsigned int>& meas_dimensions);
  
  //  For realtime use, this function should be called before calling measUpdate
  /*  @param dimension of the measurement models which is
      going to be used
  */
  void AllocateMeasModelExt( const unsigned int& meas_dimensions);
  
private:
  struct MeasUpdateVariablesAdapt
  {
    SymmetricMatrix _R;
    SymmetricMatrix _Q;
    SymmetricMatrix _W;
    Matrix _Optimal_S;
    ColumnVector _Autocorrelation_C;
    // C[0] = 
    MeasUpdateVariablesExt() {};
    MeasUpdateVariablesExt(unsigned int meas_dimension, unsigned int state_dimension):
      _R(meas_dimension)
    , _H(meas_dimension,state_dimension)
    , _Z(meas_dimension)
{};
  }; //struct
  
  
};  // class
}  // End namespace BFL
#endif  // __ADAPTIVE_EXTENDED_KALMAN_FILTER__
