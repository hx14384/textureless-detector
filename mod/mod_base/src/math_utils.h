/*
  Computer Vision Group, Department of Computer Science. University of      
  Bristol. This code can not be used or copied without express permission.   
*/

#ifndef MATH_UTILS_H_
#define MATH_UTILS_H_
#include <TooN/TooN.h>
#include <TooN/SVD.h>
#include <assert.h>
#include <vector>
#include <cvd/image_ref.h>
/**
 @file math_utils.h
 @brief Thera are some math helping functions.
 */
template< class T > inline T sum ( const std::vector<T> & v );
template<int S, class Accessor> inline double sum( const TooN::Vector<S, double, Accessor>& v );
template<int S, class Accessor> inline double max( const TooN::Vector<S, double, Accessor>& v );
template<int S, class Accessor> inline double norm2( const TooN::Vector<S, double, Accessor>& v );
template<int M> TooN::Matrix<M, M> inv( const TooN::Matrix<M, M>& A );
inline double norm2square( const TooN::Vector<3>& A );
inline TooN::Matrix<3, 3> V2MforCrossProduct( TooN::Vector<3> & t );

struct Matrix3
{
    inline static TooN::Matrix<3> inv( const TooN::Matrix<3> & mat );
    inline static TooN::Matrix<3> invSym( const TooN::Matrix<3> & mat );
    inline static double det( const TooN::Matrix<3> & mat );
};

struct Matrix2
{
    inline static TooN::Matrix<2> inv( const TooN::Matrix<2> & mat );
    inline static TooN::Matrix<2> invSym( const TooN::Matrix<2> & mat );
    inline static double det( const TooN::Matrix<2> & mat );
};

struct SymMatrix
{
    inline static void MakeSymMatrix6( TooN::Matrix<6>& m6, const TooN::Vector<6>& v6 );
    inline static void SumSymMatrix6( TooN::Matrix<6>& m6, const TooN::Vector<6>& v6, double weight = 1.0 );
    inline static void LLtoFullMatrix6( TooN::Matrix<6>& m6 );
};

template< class T > inline T sum ( const std::vector<T> & v )
{
  T result = 0;
  for( int i = 0; i < v.size(); ++i )
    result += v[i];
  return result;
};

template<int S, class Accessor> inline double sum( const TooN::Vector<S, double, Accessor>& v )
{
  double result = 0.0;
  for( int i = 0; i < v.size(); ++i )
    result += v[i];
  return result;
};

template<int S, class Accessor> inline double max( const TooN::Vector<S, double, Accessor>& v )
{
  double maximum = -1.7E308;
  for( int i = 0; i < v.size(); ++i )
    maximum = ( v[i] > maximum ) ? v[i] : maximum;
  return maximum;
};

template<int S, class Accessor> inline double norm2( const TooN::Vector<S, double, Accessor>& v )
{
  double norm = 0.0;
  for( int i = 0; i < v.size(); ++i )
    norm += v[i] * v[i];
  return sqrt( norm );
};

template<int M> TooN::Matrix<M, M> inv( const TooN::Matrix<M, M>& A )
{
  TooN::SVD<M, M> svdA( A );
  return svdA.get_pinv( 1e6 );
};

inline double norm2square( const TooN::Vector<3>&A )
{
  return A[0] * A[0] + A[1] * A[1] + A[2] * A[2];
};

inline TooN::Matrix<3, 3> V2MforCrossProduct( TooN::Vector<3> & t )
{
  //	double data[9] = {0.0, -t[2], t[1], t[2], 0.0, -t[0], -t[1], t[0], 0.0};
  TooN::Matrix<3, 3> tx;
  tx[0] = TooN::makeVector( 0.0, -t[2], t[1] );
  tx[1] = TooN::makeVector( t[2], 0.0, -t[0] );
  tx[2] = TooN::makeVector( -t[1], t[0], 0.0 );
  return tx;
};

inline double Matrix3::det( const TooN::Matrix<3> & mat )
{
  return mat[0][0] * mat[1][1] * mat[2][2] - mat[0][0] * mat[1][2] * mat[2][1] - mat[1][0] * mat[0][1]
      * mat[2][2] + mat[1][0] * mat[0][2] * mat[2][1] + mat[2][0] * mat[0][1] * mat[1][2] - mat[2][0] * mat[0][2]
      * mat[1][1]; 
};

inline TooN::Matrix<3> Matrix3::inv( const TooN::Matrix<3> & mat )
{
  TooN::Matrix<3> result;
  double detM = det(mat);
  assert( detM != 0 );
  double oneOverDetM = 1.f / detM;
  result[0] = TooN::makeVector( mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1], mat[0][2] * mat[2][1] - mat[0][1]
      * mat[2][2], mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1] );
  result[1] = TooN::makeVector( mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2], mat[0][0] * mat[2][2] - mat[0][2]
      * mat[2][0], mat[0][2] * mat[1][0] - mat[0][0] * mat[1][2] );
  result[2] = TooN::makeVector( mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0], mat[0][1] * mat[2][0] - mat[0][0]
      * mat[2][1], mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0] );
  result = oneOverDetM * result;
  return result;
};

inline TooN::Matrix<3> Matrix3::invSym( const TooN::Matrix<3> & mat )
{
  TooN::Matrix<3> result;
  double detM = det(mat);
  assert( detM != 0 );
  double oneOverDetM = 1.f / detM;
  result[2] = TooN::makeVector( mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0], mat[1][0] * mat[2][0] - mat[0][0]
      * mat[2][1], mat[0][0] * mat[1][1] - mat[1][0] * mat[1][0] );
  result[1] = TooN::makeVector( mat[2][1] * mat[2][0] - mat[1][0] * mat[2][2], mat[0][0] * mat[2][2] - mat[2][0]
      * mat[2][0], result[2][1] );
  result[0] = TooN::makeVector( mat[1][1] * mat[2][2] - mat[2][1] * mat[2][1], result[1][0], result[2][0] );
  result = oneOverDetM * result;
  return result;
};

inline double Matrix2::det( const TooN::Matrix<2> & mat )
{
  return  mat[0][1] * mat[1][0] - mat[0][0] * mat[1][1];
};

inline TooN::Matrix<2> Matrix2::inv( const TooN::Matrix<2> & mat )
{
  TooN::Matrix<2> result;
  double detM = det(mat);
  assert( detM != 0 );
  double oneOverDetM = 1.f / detM;
  result[0] = TooN::makeVector( -mat[1][1], mat[0][1] );
  result[1] = TooN::makeVector( mat[1][0], -mat[0][0] );
  result = oneOverDetM * result;
  return result;
}
;

inline TooN::Matrix<2> Matrix2::invSym( const TooN::Matrix<2> & mat )
{
  TooN::Matrix<2> result;
  double detM = det(mat);
  assert( detM != 0 );
  double oneOverDetM = 1.f / detM;
  result[0] = TooN::makeVector( -mat[1][1], mat[1][0] );
  result[1] = TooN::makeVector( mat[1][0], -mat[0][0] );
  result = oneOverDetM * result;
  return result;
}
;

inline void SymMatrix::MakeSymMatrix6( TooN::Matrix<6>& m6, const TooN::Vector<6>& v6 )
{
  m6[5] = v6[5] * v6;
  m6[4].slice( 0, 5 ) = v6[4] * v6.slice( 0, 5 );
  m6[3].slice( 0, 4 ) = v6[3] * v6.slice( 0, 4 );
  m6[2].slice( 0, 3 ) = v6[2] * v6.slice( 0, 3 );
  m6[1].slice( 0, 2 ) = v6[1] * v6.slice( 0, 2 );
  m6[0][0] = v6[0] * v6[0];
}
;

inline void SymMatrix::SumSymMatrix6( TooN::Matrix<6>& m6, const TooN::Vector<6>& v6, double weight )
{
  m6[5] += weight * v6[5] * v6;
  m6[4].slice( 0, 5 ) += weight * v6[4] * v6.slice( 0, 5 );
  m6[3].slice( 0, 4 ) += weight * v6[3] * v6.slice( 0, 4 );
  m6[2].slice( 0, 3 ) += weight * v6[2] * v6.slice( 0, 3 );
  m6[1].slice( 0, 2 ) += weight * v6[1] * v6.slice( 0, 2 );
  m6[0][0] += weight * v6[0] * v6[0];
}
;

inline void SymMatrix::LLtoFullMatrix6( TooN::Matrix<6>& m6 )
{
  m6[4][5] = m6[5][4];
  m6[3].slice( 4, 2 ) = m6.T()[3].slice( 4, 2 );
  m6[2].slice( 3, 3 ) = m6.T()[2].slice( 3, 3 );
  m6[1].slice( 2, 4 ) = m6.T()[1].slice( 2, 4 );
  m6[0].slice( 1, 5 ) = m6.T()[0].slice( 1, 5 );
}
;

inline bool isColinear( const TooN::Vector<3> & v3One, const TooN::Vector<3> & v3Two, const TooN::Vector<3> & v3Three )
{
  return fabs( (v3One^v3Two)*v3Three ) < 0.0000001 ;
}
;

inline bool isColinear( const TooN::Vector<2> & v2One, const TooN::Vector<2> & v2Two, const TooN::Vector<2> & v2Three )
{
  return fabs( (v2Three[1]-v2One[1])*(v2Two[0]-v2One[0]) - (v2Two[1]-v2One[1])*(v2Three[0]-v2One[0]) ) < 0.0000001;
}
;

inline bool isColinear( const CVD::ImageRef & irOne, const CVD::ImageRef & irTwo, const CVD::ImageRef & irThree )
{
  return fabs( (irThree.y-irOne.y)*(irTwo.x-irOne.x) - (irTwo.y-irOne.y)*(irThree.x-irOne.x) ) < 0.0000001;
}
;
#endif                           /*MATH_UTILS_H_*/
