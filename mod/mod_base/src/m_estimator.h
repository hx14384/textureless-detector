/*
  Computer Vision Group, Department of Computer Science. University of      
  Bristol. This code can not be used or copied without express permission.   
*/

#ifndef M_ESTIMATOR_H
#define M_ESTIMATOR_H
#include <vector>
#include <algorithm>
#include <assert.h>
#include <cmath>
#include "DefineTypes.h"
/**
 @file m_estimator.h
 @brief There are weighting functions for M-Robust Estimator.
 */
struct EdgeTrackerWeight
{
    inline static REAL_TYPE weight( REAL_TYPE rError, const REAL_TYPE rSigma );
    inline static REAL_TYPE estimateSigma( const std::vector<REAL_TYPE> &vrError );
};

inline REAL_TYPE EdgeTrackerWeight::weight( REAL_TYPE rError, const REAL_TYPE rSigma )
{
  return 1.f / ( rSigma + fabs( rError ) );
}
;

inline REAL_TYPE EdgeTrackerWeight::estimateSigma( const std::vector<REAL_TYPE> &vrError )
{
  std::vector<REAL_TYPE> vrErrorDummy = vrError;
  assert( vrErrorDummy.size() > 0 );
  std::sort( vrErrorDummy.begin(), vrErrorDummy.end() );
  REAL_TYPE rMedian = vrErrorDummy[vrErrorDummy.size() / 2];
  std::vector<REAL_TYPE> vrAbsDelta;
  vrAbsDelta.reserve( vrErrorDummy.size() );
  for( std::vector<REAL_TYPE>::iterator it = vrErrorDummy.begin(); it != vrErrorDummy.end(); ++it )
    vrAbsDelta.push_back( std::fabs( *it - rMedian ) );
  std::sort( vrAbsDelta.begin(), vrAbsDelta.end() );
  REAL_TYPE rMAD = vrAbsDelta[vrAbsDelta.size() / 2];
  return ( 1.4826 * rMAD );
}

struct Tukey
{
    inline static REAL_TYPE weight( REAL_TYPE rErrorSquared, const REAL_TYPE rSigmaSquared );
    inline static REAL_TYPE estimateSigmaSquared( std::vector<REAL_TYPE> & vrError );
    inline static REAL_TYPE objectiveScore( REAL_TYPE rErrorSquared, const REAL_TYPE rSigmaSquared );
    inline static REAL_TYPE estimateSigmaSquaredR( const std::vector<REAL_TYPE> &vrError, const std::vector<REAL_TYPE> &vrErrorSquared );
};

inline REAL_TYPE Tukey::weight( REAL_TYPE rErrorSquared, const REAL_TYPE rSigmaSquared )
{
  if( rErrorSquared > rSigmaSquared ) return 0;
  else
  {
    REAL_TYPE rSqrtWeight = 1.0 - (rErrorSquared / rSigmaSquared);
    return rSqrtWeight * rSqrtWeight;
  }
}
;

inline REAL_TYPE Tukey::estimateSigmaSquaredR( const std::vector<REAL_TYPE> &vrError, const std::vector<REAL_TYPE> &vrErrorSquared )
{
  int n = vrError.size();
  REAL_TYPE sum=0.0, sumSquared=0.0;
  for( int i = 0; i < n; ++i )
  {
    sum += vrError[i];
    sumSquared += vrErrorSquared[i];
  }
  return 4.6851*4.6851*(sumSquared - sum*sum/n )/n;
}
;

inline REAL_TYPE Tukey::estimateSigmaSquared( std::vector<REAL_TYPE> & vrError )
{
  assert( vrError.size() > 0 );
  std::sort( vrError.begin(), vrError.end() );
  REAL_TYPE rMedian = vrError[ vrError.size() / 2 ];
  for( std::vector<REAL_TYPE>::iterator it = vrError.begin(); it != vrError.end(); ++it )
    (*it) = fabs( (*it) - rMedian );
  std::sort( vrError.begin(), vrError.end() );
  REAL_TYPE rMAD = vrError[ vrError.size() / 2 ];
  return (1.4826*4.6851) * (1.4826*4.6851) * rMAD * rMAD;
}
;

inline REAL_TYPE Tukey::objectiveScore( REAL_TYPE rErrorSquared, const REAL_TYPE rSigmaSquared )
{
  if( rErrorSquared > rSigmaSquared )
    return 1.0;
  REAL_TYPE rBraces = 1.0 - rErrorSquared/rSigmaSquared;
  return (1.0 - rBraces*rBraces*rBraces );
}
;
#endif // M_ESTIMATOR_H
