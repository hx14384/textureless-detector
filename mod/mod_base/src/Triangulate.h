/*
  Computer Vision Group, Department of Computer Science. University of      
  Bristol. This code can not be used or copied without express permission.   
*/

#ifndef TRIANGULATE_H
#define TRIANGULATE_H

#include <TooN/TooN.h>
#include <TooN/SVD.h>
#include <TooN/se3.h>
/**
 @file Triangulate.h
 @brief Getting a 3D point by doing triangulation from a transformation matrix and its position in from 2 views.
 */
struct Triangulate
{
  /**
   @param[in] se3Xa2Xb is a rigid body transformation.
   @param[in] v2XaInCam is a homogenous position in a camera plane from the 1st view.
   @param[in] v2XbInCam is a homogenous position in a camera plane from the 2nd view.
   @return a 3D position.
   */
  inline static TooN::Vector<3> compute( const TooN::SE3<> & se3Xa2Xb, const TooN::Vector<2> & v2XaInCam, const TooN::Vector<2> & v2XbInCam );
};

inline TooN::Vector<3> Triangulate::compute( const TooN::SE3<> & se3Xa2Xb, const TooN::Vector<2> & v2XaInCam, const TooN::Vector<2> & v2XbInCam )
{
// X = [X; Y; Z; 1];
// v2XxInCam = [Ux; Vx] = normalize( v3XxInCam );
// v3XaInCam = [I|0] X; --> [-1 0 Ua 0;0 -1 Va 0] X = 0;
// v3XbInCam = [R|t] X = P; --> [ Ub P[2] - P[0]; Vb P[2] - P[1] ] X = 0;
// AX = 0; Find X from SVD.
  TooN::Matrix<3,4> P;
  P.slice<0,0,3,3>() = se3Xa2Xb.get_rotation().get_matrix();
  P.slice<0,3,3,1>() = se3Xa2Xb.get_translation().as_col();

  TooN::Matrix<4> A;
  A[0] = TooN::makeVector(-1.0, 0.0, v2XaInCam[0], 0.0);
  A[1] = TooN::makeVector(0.0, -1.0, v2XaInCam[1], 0.0);
  A[2] = v2XbInCam[0]*P[2] - P[0];
  A[3] = v2XbInCam[1]*P[2] - P[1];
  
  TooN::SVD<4,4> svd(A);
  TooN::Vector<4> v4 = svd.get_VT()[3]; // related to smallest singular value.
  if( v4[3] == 0.0 ) // it is an infinity point.
    v4[3] = 0.000001;
  return project(v4);
};
#endif
