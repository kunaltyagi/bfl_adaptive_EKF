// $Id: adaptiveextendedkalmanfilter.cpp 30022 2015-09-06 16:10:30Z kunaltyagi $
// Copyright (C) 2003 Kunal Tyagi <last dot first at live dot com>
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

#include "adaptiveextendedkalmanfilter.h"
#include <cmath>

namespace BFL
{
using namespace MatrixWrapper;
  
#define AnalyticSys    AnalyticSystemModelGaussianUncertainty
#define AnalyticMeas   AnalyticMeasurementModelGaussianUncertainty
  
AdaptiveExtendedKalmanFilter::AdaptiveExtendedKalmanFilter(Gaussian* prior, double Nr, double Nq)
  : // KalmanFilter(prior)
    Filter<ColumnVector,ColumnVector>(prior)
  , _Mu_new(prior->DimensionGet())
  , _Sigma_new(prior->DimensionGet())
  , _Sigma_temp(prior->DimensionGet(),prior->DimensionGet())
  , _Sigma_temp_par(prior->DimensionGet(),prior->DimensionGet())  // kalmanfilter end, ekf start
  , _x(prior->DimensionGet())
  , _J(prior->DimensionGet())
  , _F(prior->DimensionGet(),prior->DimensionGet())
  , _Q(prior->DimensionGet())
  // changes for adaptive
  , _Nr(Nr)
  , _Nq(Nq)
{
  // create posterior dencity
  _post = new Gaussian(*prior);
  _alpha1 = (_Nr - 1)/_Nr;  // reduce loss in signigicant digits
  _alpha1 = (_Nq - 1)/_Nq;
  _firstData = true;
}

AdaptiveExtendedKalmanFilter::~AdaptiveExtendedKalmanFilter()
{
  delete _post;
}

// kf start
void
AdaptiveExtendedKalmanFilter::AllocateMeasModel(const vector<unsigned int>& meas_dimensions)
{
  unsigned int meas_dimension;
  for(int i = 0 ; i< meas_dimensions.size(); i++)
  {
    // find if variables with size meas_sizes[i] are already allocated
    meas_dimension = meas_dimensions[i];
    _mapMeasUpdateVariables_it =  _mapMeasUpdateVariables.find(meas_dimension);
    if( _mapMeasUpdateVariables_it == _mapMeasUpdateVariables.end())
    {
      //variables with size z.rows() not allocated yet
      _mapMeasUpdateVariables_it = (_mapMeasUpdateVariables.insert
          (std::pair<unsigned int, MeasUpdateVariables>( meas_dimension,MeasUpdateVariables(meas_dimension,_Mu_new.rows()) ))).first;
    }
  }
}

void
AdaptiveExtendedKalmanFilter::AllocateMeasModel(const unsigned int& meas_dimension)
{
  // find if variables with size meas_sizes[i] are already allocated
  _mapMeasUpdateVariables_it =  _mapMeasUpdateVariables.find(meas_dimension);
  if( _mapMeasUpdateVariables_it == _mapMeasUpdateVariables.end())
  {
    //variables with size z.rows() not allocated yet
    _mapMeasUpdateVariables_it = (_mapMeasUpdateVariables.insert
        (std::pair<unsigned int, MeasUpdateVariables>( meas_dimension,MeasUpdateVariables(meas_dimension,_Mu_new.rows()) ))).first;
  }
}

void
AdaptiveExtendedKalmanFilter::CalculateSysUpdate(const ColumnVector& J, const Matrix& F, const SymmetricMatrix& Q)
{
  _Sigma_temp = F * ( (Matrix)_post->CovarianceGet() * F.transpose());
  _Sigma_temp += (Matrix)Q;
  _Sigma_temp.convertToSymmetricMatrix(_Sigma_new);
  
  // set new state gaussian
  PostMuSet   ( J );
  PostSigmaSet( _Sigma_new );
}

void
AdaptiveExtendedKalmanFilter::CalculateMeasUpdate(const ColumnVector& z, const ColumnVector& Z, const Matrix& H, const SymmetricMatrix& R)
{
  // allocate measurement for z.rows() if needed
  AllocateMeasModel(z.rows());
  
  (_mapMeasUpdateVariables_it->second)._postHT =   (Matrix)(_post->CovarianceGet()) * H.transpose() ;
  (_mapMeasUpdateVariables_it->second)._S_Matrix =  H * (_mapMeasUpdateVariables_it->second)._postHT;
  
  // calcutate new state gaussian
  // Mu = expectedValue + K*(z-Z)
  (_mapMeasUpdateVariables_it->second)._innov = z-Z;
  
  // changes for adaptive begin (part 1)
  (_mapMeasUpdateVaraiblesAdapt_it->second)._OptError *= _alpha1;
  (_mapMeasUpdateVaraiblesAdapt_it->second)._OptError += (_mapMeasUpdateVariables_it->second)._innov/_Nr;
  
  (_mapMeasUpdateVaraiblesAdapt_it->second)._DeltaR = ((_mapMeasUpdateVariables_it->second)._innov - (_mapMeasUpdateVaraiblesAdapt_it->second)._OptError);
  (_mapMeasUpdateVaraiblesAdapt_it->second)._DeltaR *= (_mapMeasUpdateVaraiblesAdapt_it->second)._DeltaR.transpose()/(Nr - 1);
  (_mapMeasUpdateVaraiblesAdapt_it->second)._DeltaR -= (_mapMeasUpdateVariables_it->second)._S_Matrix/_Nr;
  
  ((Matrix)(R*_alpha1 + (_mapMeasUpdateVaraiblesAdapt_it->second)._DeltaR)).convertToSymmetricMatrix(R);
  for (int i = 0; i < R.rows(); ++i)
  {
    for (int j = 0; j < R.columns(); ++j)
    {
      if (R(i,j) < 0)
      {
        R(i,j) = -R(i,j);
      }
    }
  }
  // changes for adaptive end (part 1)
  
  (_mapMeasUpdateVariables_it->second)._S_Matrix += (Matrix)R;
  
  // _K = covariance * H' * S(-1)
  (_mapMeasUpdateVariables_it->second)._K =  (_mapMeasUpdateVariables_it->second)._postHT * ( (_mapMeasUpdateVariables_it->second)._S_Matrix.inverse());
  
  _Mu_new  =  (_mapMeasUpdateVariables_it->second)._K * (_mapMeasUpdateVariables_it->second)._innov  ;
  
  // changes for adaptive begin (part 2)
  _alpha2 = (_Nq - 1)/_Nq;
  (_mapMeasUpdateVaraiblesAdapt_it->second).OptOmega *= _alpha2;
  (_mapMeasUpdateVaraiblesAdapt_it->second).OptOmega += _Mu_new/_Nq;
  
  (_mapMeasUpdateVaraiblesAdapt_it->second)._DeltaQ = (_Mu_new._innov - (_mapMeasUpdateVaraiblesAdapt_it->second)._OptOmega);
  (_mapMeasUpdateVaraiblesAdapt_it->second)._DeltaQ *= (_mapMeasUpdateVaraiblesAdapt_it->second)._DeltaQ.transpose()/(Nq - 1);
  // changes for adaptive end (part 2)
  
  _Mu_new  +=  _post->ExpectedValueGet() ;
  
  // Sigma = post - K*H*post
  _Sigma_temp = (_post->CovarianceGet());
  _Sigma_temp_par = (_mapMeasUpdateVariables_it->second)._K * H ;
  _Sigma_temp -=  _Sigma_temp_par * (Matrix)(_post->CovarianceGet());
  // convert to symmetric matrix
  _Sigma_temp.convertToSymmetricMatrix(_Sigma_new);
  
  // changes for adaptive begin (part 3)
  ((Matrix)(_mapMeasUpdateVaraiblesAdapt_it->second)._Q*_alpha2 + (_mapMeasUpdateVaraiblesAdapt_it->second)._DeltaQ)).convertToSymmetricMatrix((_mapMeasUpdateVaraiblesAdapt_it->second)._Q);  // if required, take the diagonal elements only
  for (int i = 0; i < _mapMeasUpdateVaraiblesAdapt_it->second)._Q.rows(); ++i)
  {
    for (int j = 0; j < _mapMeasUpdateVaraiblesAdapt_it->second)._Q.columns(); ++j)
    {
      if (_mapMeasUpdateVaraiblesAdapt_it->second)._Q(i,j) < 0)
      {
        _mapMeasUpdateVaraiblesAdapt_it->second)._Q(i,j) = -_mapMeasUpdateVaraiblesAdapt_it->second)._Q(i,j);
      }
    }
  }
  // changes for adaptive end (part 3)
  
  // set new state gaussian
  PostMuSet   ( _Mu_new );
  PostSigmaSet( _Sigma_new );
  
  /*
  cout << "H:\n" << H << endl;
  cout << "R:\n" << R << endl;
  cout << "Z:\n" << Z << endl;
  cout << "inov:\n" << z-Z << endl;
  cout << "S:\n" << S << endl;
  cout << "S.inverse:\n" << S.inverse() << endl;
  cout << "K:\n" << K << endl;
  cout << "Mu_new:\n" << Mu_new << endl;
  cout << "sigma_new\n" << Sigma_new << endl;
  */
}
 
bool
AdaptiveExtendedKalmanFilter::UpdateInternal(SystemModel<ColumnVector>* const sysmodel,
                             const ColumnVector& u,
                             MeasurementModel<ColumnVector,ColumnVector>* const measmodel,
                             const ColumnVector& z, const ColumnVector& s)
{
  if (sysmodel != NULL)
  {
    SysUpdate(sysmodel,u);
  }
  if (measmodel != NULL)
  {
    MeasUpdate(measmodel,z,s);
  }
  return true;
}

void
AdaptiveExtendedKalmanFilter::PostSigmaSet( const SymmetricMatrix& s)
{
  dynamic_cast<Gaussian *>(_post)->CovarianceSet(s);
}
 
void
AdaptiveExtendedKalmanFilter::PostMuSet( const ColumnVector& c)
{
  dynamic_cast<Gaussian *>(_post)->ExpectedValueSet(c);
}

Gaussian*
AdaptiveExtendedKalmanFilter::PostGet()
{
  return (Gaussian*)Filter<ColumnVector,ColumnVector>::PostGet();
}
// kf ends

// ekf starts

void
AdaptiveExtendedKalmanFilter::AllocateMeasModelAdapt(const vector<unsigned int>& meas_dimensions)
{
  unsigned int meas_dimension;
  for(int i = 0 ; i< meas_dimensions.size(); i++)
  {
    // find if variables with size meas_sizes[i] are already allocated
    meas_dimension = meas_dimensions[i];
    _mapMeasUpdateVariablesAdapt_it =  _mapMeasUpdateVariablesAdapt.find(meas_dimension);
    if( _mapMeasUpdateVariablesAdapt_it == _mapMeasUpdateVariablesAdapt.end())
    {
      //variables with size z.rows() not allocated yet
      _mapMeasUpdateVariablesAdapt_it = (_mapMeasUpdateVariablesAdapt.insert
          (std::pair<unsigned int, MeasUpdateVariablesAdapt>( meas_dimension,MeasUpdateVariablesAdapt(meas_dimension,_x.rows()) ))).first;
    }
  }
}

void
AdaptiveExtendedKalmanFilter::AllocateMeasModelAdapt(const unsigned int& meas_dimension)
{
  // find if variables with size meas_sizes[i] are already allocated
  _mapMeasUpdateVariablesAdapt_it =  _mapMeasUpdateVariablesAdapt.find(meas_dimension);
  if( _mapMeasUpdateVariablesAdapt_it == _mapMeasUpdateVariablesAdapt.end())
  {
    //variables with size z.rows() not allocated yet
    _mapMeasUpdateVariablesAdapt_it = (_mapMeasUpdateVariablesAdapt.insert
        (std::pair<unsigned int, MeasUpdateVariablesAdapt>( meas_dimension,MeasUpdateVariablesAdapt(meas_dimension,_x.rows()) ))).first;
  }
}

void
AdaptiveExtendedKalmanFilter::AllocateMeasModelExt(const vector<unsigned int>& meas_dimensions)
{
  unsigned int meas_dimension;
  for(int i = 0 ; i< meas_dimensions.size(); i++)
  {
    // find if variables with size meas_sizes[i] are already allocated
    meas_dimension = meas_dimensions[i];
    _mapMeasUpdateVariablesExt_it =  _mapMeasUpdateVariablesExt.find(meas_dimension);
    if( _mapMeasUpdateVariablesExt_it == _mapMeasUpdateVariablesExt.end())
    {
      //variables with size z.rows() not allocated yet
      _mapMeasUpdateVariablesExt_it = (_mapMeasUpdateVariablesExt.insert
          (std::pair<unsigned int, MeasUpdateVariablesExt>( meas_dimension,MeasUpdateVariablesExt(meas_dimension,_x.rows()) ))).first;
    }
  }
}

void
AdaptiveExtendedKalmanFilter::AllocateMeasModelExt(const unsigned int& meas_dimension)
{
  // find if variables with size meas_sizes[i] are already allocated
  _mapMeasUpdateVariablesExt_it =  _mapMeasUpdateVariablesExt.find(meas_dimension);
  if( _mapMeasUpdateVariablesExt_it == _mapMeasUpdateVariablesExt.end())
  {
    //variables with size z.rows() not allocated yet
    _mapMeasUpdateVariablesExt_it = (_mapMeasUpdateVariablesExt.insert
        (std::pair<unsigned int, MeasUpdateVariablesExt>( meas_dimension,MeasUpdateVariablesExt(meas_dimension,_x.rows()) ))).first;
  }
}
 
void
AdaptiveExtendedKalmanFilter::SysUpdate(SystemModel<ColumnVector>* const sysmodel,
                                const ColumnVector& u)
{
  _x = _post->ExpectedValueGet();
  _J = ((AnalyticSys*)sysmodel)->PredictionGet(u,_x);
  _F = ((AnalyticSys*)sysmodel)->df_dxGet(u,_x);
  _Q = ((AnalyticSys*)sysmodel)->CovarianceGet(u,_x);
  
  // @TODO: confirm this
  CalculateSysUpdate(_J, _F, (_mapMeasUpdateVaraiblesAdapt_it->second)._Q);
}
 
void
AdaptiveExtendedKalmanFilter::MeasUpdate(MeasurementModel<ColumnVector,ColumnVector>* const measmodel,
                                 const ColumnVector& z,
                                 const ColumnVector& s)
{
  // allocate measurement for z.rows() if needed
  AllocateMeasModelExt(z.rows());

  _x = _post->ExpectedValueGet();
  (_mapMeasUpdateVariablesExt_it->second)._Z = ((AnalyticMeas*)measmodel)->PredictionGet(s,_x);
  (_mapMeasUpdateVariablesExt_it->second)._H = ((AnalyticMeas*)measmodel)->df_dxGet(s,_x);
  // remove this
  if (true == _firstData)
  {
    (_mapMeasUpdateVariablesExt_it->second)._R = ((AnalyticMeas*)measmodel)->CovarianceGet(s,_x);
    (_mapMeasUpdateVariablesAdapt_it->second)._Q = _Q;
    _firstData = false;
  }
 
  CalculateMeasUpdate(z, (_mapMeasUpdateVariablesExt_it->second)._Z, (_mapMeasUpdateVariablesExt_it->second)._H, (_mapMeasUpdateVariablesExt_it->second)._R);
}
// ekf ends

} // end namespace BFL
