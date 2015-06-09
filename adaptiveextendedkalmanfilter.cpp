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
  
AdaptiveExtendedKalmanFilter::AdaptiveExtendedKalmanFilter(Gaussian* prior)
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
{
  // create posterior dencity
  _post = new Gaussian(*prior);
}

AdaptiveExtendedKalmanFilter::~AdaptiveExtendedKalmanFilter()
{
  delete _post;
}

