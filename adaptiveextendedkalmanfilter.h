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
  // @TODO: replace all these by a friend class AdaptiveExtendedKalmanFilter in KalmanFilter and ExtendedKalmanFilter
  struct MeasUpdateVariables
  {
    Matrix _S_Matrix;
    Matrix _K;
    ColumnVector _innov;
    Matrix _postHT;
    MeasUpdateVariables() {};
    MeasUpdateVariables(unsigned int meas_dimension, unsigned int state_dimension):
      _S_Matrix(meas_dimension,meas_dimension)
    , _K(state_dimension,meas_dimension)
    , _innov(meas_dimension)
    , _postHT(state_dimension,meas_dimension)
{};
  }; //struct
  struct MeasUpdateVariablesExt
  {
    SymmetricMatrix _R;
    Matrix _H;
    ColumnVector _Z;
    MeasUpdateVariablesExt() {};
    MeasUpdateVariablesExt(unsigned int meas_dimension, unsigned int state_dimension):
      _R(meas_dimension)
    , _H(meas_dimension,state_dimension)
    , _Z(meas_dimension)
{};
  }; //struct
  struct MeasUpdateVariablesAdapt
  {
    ColumnVector _OptError;
    ColumnVector _Omega;
    ColumnVector _OptOmega;
    SymmetricMatrix _Q;
    Matrix _DeltaR;
    Matric _DeltaQ;
  }
  
protected:
  // variables to avoid allocation during update calls
  ColumnVector  _Mu_new;
  SymmetricMatrix _Sigma_new;
  Matrix _Sigma_temp;
  Matrix _Sigma_temp_par;
  Matrix _SMatrix;
  Matrix _K;
  std::map<unsigned int, MeasUpdateVariables> _mapMeasUpdateVariables;
  std::map<unsigned int, MeasUpdateVariables>::iterator _mapMeasUpdateVariables_it;
  
  void PostSigmaSet( const MatrixWrapper::SymmetricMatrix& s);
  void PostMuSet( const MatrixWrapper::ColumnVector& c);
  void CalculateSysUpdate(const MatrixWrapper::ColumnVector& J, const MatrixWrapper::Matrix& F, const MatrixWrapper::SymmetricMatrix& Q);
  void CalculateMeasUpdate(const MatrixWrapper::ColumnVector& z, const MatrixWrapper::ColumnVector& Z, const MatrixWrapper::Matrix& H, const MatrixWrapper::SymmetricMatrix& R);
  // virtual void SysUpdate(SystemModel<MatrixWrapper::ColumnVector>* const sysmodel,
  //                        const MatrixWrapper::ColumnVector& u) = 0;
  // virtual void MeasUpdate(MeasurementModel<MatrixWrapper::ColumnVector,MatrixWrapper::ColumnVector>* const measmodel,
  //                         const MatrixWrapper::ColumnVector& z,
  //                         const MatrixWrapper::ColumnVector& s) = 0;
  virtual bool UpdateInternal(SystemModel<MatrixWrapper::ColumnVector>* const sysmodel,
                              const MatrixWrapper::ColumnVector& u,
                              MeasurementModel<MatrixWrapper::ColumnVector,MatrixWrapper::ColumnVector>* const measmodel,
                              const MatrixWrapper::ColumnVector& z,
                              const MatrixWrapper::ColumnVector& s);
  
  virtual void SysUpdate(SystemModel<MatrixWrapper::ColumnVector>* const sysmodel,
                         const MatrixWrapper::ColumnVector& u);
  virtual void MeasUpdate(MeasurementModel<MatrixWrapper::ColumnVector,MatrixWrapper::ColumnVector>* const measmodel,
                          const MatrixWrapper::ColumnVector& z,
                          const MatrixWrapper::ColumnVector& s);
  // variables to avoid allocation on the heap
  ColumnVector _x;
  ColumnVector _J;
  Matrix    _F;
  SymmetricMatrix _Q;
  std::map<unsigned int, MeasUpdateVariablesExt> _mapMeasUpdateVariablesExt;
  std::map<unsigned int, MeasUpdateVariablesExt>::iterator _mapMeasUpdateVariablesExt_it;
  
  std::map<unsigned int, MeasUpdateVariablesAdapt> _mapMeasUpdateVariablesAdapt;
  std::map<unsigned int, MeasUpdateVariablesAdapt>::iterator _mapMeasUpdateVariablesAdapt_it;
  
};  // class
}  // End namespace BFL
#endif  // __ADAPTIVE_EXTENDED_KALMAN_FILTER__
