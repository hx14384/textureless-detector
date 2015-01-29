/********************************************************************************/
/**  Computer Vision Group, Department of Computer Science. University of      **/
/**  Bristol. This code can not be used or copied without express permission.  **/ 
/********************************************************************************/

#include "Homography.h"
#include <TooN/SVD.h>
#include <gvars3/instances.h>
#include <assert.h>
#include <algorithm>
#include "math_utils.h"
#include "Triangulate.h"
#define SQRT2 1.41421356
Homography::Homography() : mvse3First2Second(2)
{
  mv5NormalizeMatchedPairs.resize(5);
}

Homography::~Homography()
{

}

bool Homography::compute( const std::vector<MatchedPair> & vMatchedPairs)
{
  mvMatchedPairs = vMatchedPairs;
  TooN::Matrix<3> H;
  estimateHomography( H );
  if( estimateSE3( H ) )
  {
    triangulate( mvse3First2Second[0] );
    cal_plane();
    return true;
  }
  else
  {
    return false;
  }
  return true;
}

bool Homography::computeForRectangle( const std::vector<MatchedPair> & vMatchedPairs )
{
  mvMatchedPairs = vMatchedPairs;
  TooN::Matrix<3> H;
  estimateHomographyFromRectangle( H );
  if( estimateSE3FromRectangle( H ) )
  {
    triangulate( mvse3First2Second[0] );
    cal_plane();
    return true; 
  }
  else
  {
    return false;
  }
  return true;
}

int Homography::max (int a, int b)
{
   if (a > b) return a;
   return b;
}

int Homography::min (int a, int b)
{
   if (a < b) return a;
   return b;
}

double Homography::estimateRotationAngle(double & stdDev, std::vector<MatchedPair> & vMatchedPairs)
{
  double angle;
  // estimate the mean
  MatchedPair meanPair, minPair, maxPair;
  meanPair.v2FstView[0] = 0;
  meanPair.v2FstView[1] = 0;
  meanPair.v2SndView[0] = 0;
  meanPair.v2SndView[1] = 0;
  minPair.v2FstView[0] = 1000;
  minPair.v2FstView[1] = 1000;
  minPair.v2SndView[0] = 1000;
  minPair.v2SndView[1] = 1000;
  maxPair.v2FstView[0] = 0;
  maxPair.v2FstView[1] = 0;
  maxPair.v2SndView[0] = 0;
  maxPair.v2SndView[1] = 0;
  const int n = vMatchedPairs.size();
  for (int i = 0; i < n; i++)
  {
     meanPair.v2FstView[0] += vMatchedPairs[i].v2FstView[0];
     meanPair.v2FstView[1] += vMatchedPairs[i].v2FstView[1];
     meanPair.v2SndView[0] += vMatchedPairs[i].v2SndView[0];
     meanPair.v2SndView[1] += vMatchedPairs[i].v2SndView[1];
   //  std::cout << vMatchedPairs[i].v2FstView[0] << " " << vMatchedPairs[i].v2FstView[1] << " " << vMatchedPairs[i].v2SndView[0] << " " << vMatchedPairs[i].v2SndView[1] << std::endl;
  }
  meanPair.v2FstView[0] /= n;
  meanPair.v2FstView[1] /= n;
  meanPair.v2SndView[0] /= n;
  meanPair.v2SndView[1] /= n;

  for (int i = 0; i < n; i++)
  {
     vMatchedPairs[i].v2FstView[0] -= meanPair.v2FstView[0];
     minPair.v2FstView[0] = min (minPair.v2FstView[0], vMatchedPairs[i].v2FstView[0]);
     maxPair.v2FstView[0] = max (maxPair.v2FstView[0], vMatchedPairs[i].v2FstView[0]);

     vMatchedPairs[i].v2FstView[1] -= meanPair.v2FstView[1];
     minPair.v2FstView[1] = min (minPair.v2FstView[1], vMatchedPairs[i].v2FstView[1]);
     maxPair.v2FstView[1] = max (maxPair.v2FstView[1], vMatchedPairs[i].v2FstView[1]);

     vMatchedPairs[i].v2SndView[0] -= meanPair.v2SndView[0];
     minPair.v2SndView[0] = min (minPair.v2SndView[0], vMatchedPairs[i].v2SndView[0]);
     maxPair.v2SndView[0] = max (maxPair.v2SndView[0], vMatchedPairs[i].v2SndView[0]);

     vMatchedPairs[i].v2SndView[1] -= meanPair.v2SndView[1];
     minPair.v2SndView[1] = min (minPair.v2SndView[1], vMatchedPairs[i].v2SndView[1]);
     maxPair.v2SndView[1] = max (maxPair.v2SndView[1], vMatchedPairs[i].v2SndView[1]);
  }

  // scaled values
 /* double scale = 0;
  scale += ((double)(minPair.v2FstView[0]))/minPair.v2SndView[0];
  scale += ((double)(minPair.v2FstView[1]))/minPair.v2SndView[1];
  scale += ((double)(maxPair.v2FstView[0]))/maxPair.v2SndView[0];
  scale += ((double)(maxPair.v2FstView[1]))/maxPair.v2SndView[1];
  scale /= 4;*/
  double scale = 1;

  TooN::Matrix<> B (n,2);
  for (int i = 0; i < n; i++)
  {
      B[i][0] = vMatchedPairs[i].v2SndView[0] * scale;
      B[i][1] = vMatchedPairs[i].v2SndView[1] * scale;
  }
  // svd to get the pinv
//  TooN::Matrix<2,n> M(B);
/*  TooN::SVD<> svdM(B);
  TooN::Matrix<> x = svdM.get_pinv();  
  angle = -1*acos(x(1,1));*/

  double anSum = 0;
  std::vector<double> angles;
  angles.reserve (n);
  for (int i = 0; i < n; i++)
  {
     double an1 = atan2 (vMatchedPairs[i].v2FstView[1],vMatchedPairs[i].v2FstView[0]);
     double an2 = atan2 (B[i][1],B[i][0]);
   //  std::cout << vMatchedPairs[i].v2FstView[0] << " " << vMatchedPairs[i].v2FstView[1] << " " << vMatchedPairs[i].v2SndView[0] << " " << vMatchedPairs[i].v2SndView[1] << " " << B[i][0] << " " << B[i][1] << " " << an1 << " " << an2 << std::endl;
     double diff = an2-an1;
     if (diff > 3.14)
        diff = 2*3.14 - diff;
     if (diff < -3.14)
        diff = 2*3.14 + diff;
     anSum += diff;
     angles.push_back (diff);
  }
  angle = anSum / n;

  // calculate standard deviation
  double anSD = 0;
  for (int i = 0; i < n; i++)
  {
     anSD = (angles[i]-angle)*(angles[i]-angle);
  }  
  stdDev = sqrt (anSD/n);
 // std::cout << x(1,1) << " " << x(1,2) << " " << x(2,1) << " " << x(2,2) << std::endl;
  //std::cout << angle << " with standard deviation " << stdDev << std::endl;
  //exit(-1);
  return angle;
}

bool Homography::estimateHomography_group_of_5(TooN::Matrix<3> & H, const std::vector<MatchedPair> & vMatchedPairs)
{
  if( isDegenerate( vMatchedPairs ) )
  {
/*    for( int i = 0; i < vMatchedPairs.size(); ++i )
      std::cout << vMatchedPairs[i].v2FstView << " " << vMatchedPairs[i].v2SndView << std::endl; */
   // std::cout <<"Degenerate case ... " << std::endl;
    return false;
  }
  TooN::Matrix<3> T1, T2;
  normalise2dpts( mv5NormalizeMatchedPairs, T1, T2, vMatchedPairs );
  TooN::Matrix<10,9> mx9;
  for( int i = 0; i < 5; ++i )
  {
    const double & u0 = mv5NormalizeMatchedPairs[i].v2FstView[0];
    const double & v0 = mv5NormalizeMatchedPairs[i].v2FstView[1];
    const double & u1 = mv5NormalizeMatchedPairs[i].v2SndView[0];
    const double & v1 = mv5NormalizeMatchedPairs[i].v2SndView[1];
//    mx9[2*i] = TooN::makeVector(-u0, -v0, -1.f, 0, 0, 0, u1*u0, u1*v0, u1);
//    mx9[2*i+1] = TooN::makeVector(0, 0, 0, -u0, -v0, -1, v1*u0, v1*v0, v1);
    int idx = 2*i;
    mx9[idx][0] = -u0;    mx9[idx][1] = -v0;    mx9[idx][2] = -1.f;
    mx9[idx][3] = 0.f;    mx9[idx][4] = 0.f;    mx9[idx][5] = 0.f;
    mx9[idx][6] = u1*u0;  mx9[idx][7] = u1*v0;  mx9[idx][8] = u1;
    idx += 1;
    mx9[idx][0] = 0.f;    mx9[idx][1] = 0.f;    mx9[idx][2] = 0.f;
    mx9[idx][3] = -u0;    mx9[idx][4] = -v0;    mx9[idx][5] = -1.f;
    mx9[idx][6] = v1*u0;  mx9[idx][7] = v1*v0;  mx9[idx][8] = v1;
  }
//  m9[8] = TooN::Zeros;
  TooN::SVD<> svd(mx9);
  TooN::Vector<9> h = svd.get_VT()[8];
  H[0] = h.slice<0,3>(); H[1] = h.slice<3,3>(); H[2] = h.slice<6,3>();
  H = Matrix3::inv(T2)*H*T1;
//  std::cout << "Homography : " << H << std::endl;
//  return checkNiceHomography( H );
// Test
  if( checkNiceHomography( H ) )
  {
//    std::vector<int> vInts; vInts.push_back(0);
/*    for( int i = 0; i < 5; ++i )
    {
      std::cout << vMatchedPairs[i].v2FstView << " " << vMatchedPairs[i].v2SndView << " " << mv5NormalizeMatchedPairs[i].v2FstView << " " << mv5NormalizeMatchedPairs[i].v2SndView << std::endl;
    } */
//    estimateSE3( H, vInts, vMatchedPairs );
    return true;
  }
  else
    return false;
}

void Homography::estimateHomographyFromRectangle(TooN::Matrix<3> & H)
{
/*
  TooN::Matrix<9> m9;
  for( int i = 0; i < 4; i++ )
  {
    double & u0 = mvMatchedPairs[i].v2FstView[0];
    double & v0 = mvMatchedPairs[i].v2FstView[1];
    double & u1 = mvMatchedPairs[i].v2SndView[0];
    double & v1 = mvMatchedPairs[i].v2SndView[1];
    m9[2*i] = TooN::makeVector(-u0, -v0, -1.f, 0, 0, 0, u1*u0, u1*v0, u1);
    m9[2*i+1] = TooN::makeVector(0, 0, 0, -u0, -v0, -1, v1*u0, v1*v0, v1);
  }
  m9[8] = TooN::Zeros;

  TooN::SVD<> svd(m9);
  TooN::Vector<9> h = svd.get_VT()[8];
  H[0] = h.slice<0,3>(); H[1] = h.slice<3,3>(); H[2] = h.slice<6,3>();
*/ 
  TooN::Matrix<9> m9;
  TooN::Matrix<3> T1, T2;
  std::vector<MatchedPair> vMatchedPairs ( mvMatchedPairs.size() );
  normalise2dpts( vMatchedPairs, T1, T2, mvMatchedPairs );
  for( int i = 0; i < 4; i++ )
  {
    double & u0 = vMatchedPairs[i].v2FstView[0];
    double & v0 = vMatchedPairs[i].v2FstView[1];
    double & u1 = vMatchedPairs[i].v2SndView[0];
    double & v1 = vMatchedPairs[i].v2SndView[1];
    m9[2*i] = TooN::makeVector(-u0, -v0, -1.f, 0, 0, 0, u1*u0, u1*v0, u1);
    m9[2*i+1] = TooN::makeVector(0, 0, 0, -u0, -v0, -1, v1*u0, v1*v0, v1);
  }
  m9[8] = TooN::Zeros;

  TooN::SVD<> svd(m9);
  TooN::Vector<9> h = svd.get_VT()[8];
  H[0] = h.slice<0,3>(); H[1] = h.slice<3,3>(); H[2] = h.slice<6,3>();
  H=Matrix3::inv(T2)*H*T1; 
//  std::cout << "Homography : " << H << std::endl;
}

void Homography::estimateHomography(TooN::Matrix<3> & H)
{
  assert( mvMatchedPairs.size() >= 4 );//assume there are more than 4 pairs.
  TooN::Matrix<9> m9;
  TooN::Matrix<3> T1, T2;
  std::vector<MatchedPair> vMatchedPairs( mvMatchedPairs.size() );
  normalise2dpts( vMatchedPairs, T1, T2, mvMatchedPairs );
  for( unsigned int i = 0; i < mvMatchedPairs.size(); i++ )
  {
    double & u0 = vMatchedPairs[i].v2FstView[0];
    double & v0 = vMatchedPairs[i].v2FstView[1];
    double & u1 = vMatchedPairs[i].v2SndView[0];
    double & v1 = vMatchedPairs[i].v2SndView[1];
    m9[2*i] = TooN::makeVector(-u0, -v0, -1.f, 0, 0, 0, u1*u0, u1*v0, u1);
    m9[2*i+1] = TooN::makeVector(0, 0, 0, -u0, -v0, -1, v1*u0, v1*v0, v1);
  }
  m9[8] = TooN::Zeros;

  TooN::SVD<> svd(m9);
  TooN::Vector<9> h = svd.get_VT()[8];
  H[0] = h.slice<0,3>(); H[1] = h.slice<3,3>(); H[2] = h.slice<6,3>();
  H = Matrix3::inv(T2)*H*T1;
//  std::cout << "Homography : " << H << std::endl;
}

bool Homography::estimateHomography(TooN::Matrix<3> & H, const std::vector<MatchedPair> & vMatchedPairs, int size )
{
  int matchedSize = ( size == -1 ) ? vMatchedPairs.size() : size;
  assert( matchedSize > 4 );//assume there are more than 4 pairs.
  TooN::Matrix<> mx9( 2*matchedSize, 9 );
  TooN::Matrix<3> T1, T2;
  std::vector<MatchedPair> vCopyMatchedPairs = vMatchedPairs;
  vCopyMatchedPairs.resize(matchedSize);
  std::vector<MatchedPair> vNormalizeMatchedPairs = vCopyMatchedPairs;
  normalise2dpts( vNormalizeMatchedPairs, T1, T2, vCopyMatchedPairs );
  for( int i = 0; i < matchedSize; i++ )
  {
    const double & u0 = vNormalizeMatchedPairs[i].v2FstView[0];
    const double & v0 = vNormalizeMatchedPairs[i].v2FstView[1];
    const double & u1 = vNormalizeMatchedPairs[i].v2SndView[0];
    const double & v1 = vNormalizeMatchedPairs[i].v2SndView[1];
    mx9[2*i] = TooN::makeVector(-u0, -v0, -1.f, 0, 0, 0, u1*u0, u1*v0, u1);
    mx9[2*i+1] = TooN::makeVector(0, 0, 0, -u0, -v0, -1, v1*u0, v1*v0, v1);
  }
//  m9[8] = TooN::Zeros;
  TooN::SVD<> svd(mx9);
  TooN::Vector<9> h = svd.get_VT()[8];
  H[0] = h.slice<0,3>(); H[1] = h.slice<3,3>(); H[2] = h.slice<6,3>();
  H = Matrix3::inv(T2)*H*T1;
//  std::cout << "Homography : " << H << std::endl;
  return checkNiceHomography( H );
}

void Homography::normalise2dpts( std::vector<MatchedPair> & vNormalizeMatchedPairs, TooN::Matrix<3> & T1, TooN::Matrix<3> & T2, const std::vector<MatchedPair> & vMatchedPairs )
{
  //assume vNormalizeMatchedPairs and vMatchedPairs have the same size;
  TooN::Vector<2> v2Centroid1 = TooN::Zeros;
  TooN::Vector<2> v2Centroid2 = TooN::Zeros;
  int no_matchedPairs = vMatchedPairs.size();
  for( int i = 0; i < no_matchedPairs; ++i )
  {
    v2Centroid1 += vMatchedPairs[i].v2FstView;
    v2Centroid2 += vMatchedPairs[i].v2SndView;
  }
  v2Centroid1 /= no_matchedPairs;
  v2Centroid2 /= no_matchedPairs;
  REAL_TYPE dist1 = 0.f, dist2 = 0.f;
  TooN::Vector<2> v2;
  for( int i = 0; i < no_matchedPairs; ++i )
  {
    v2 = vMatchedPairs[i].v2FstView - v2Centroid1;
    dist1 += norm(v2);
    v2 = vMatchedPairs[i].v2SndView - v2Centroid2;
    dist2 += norm(v2);
  }
  dist1 /= no_matchedPairs; dist2 /= no_matchedPairs;
  REAL_TYPE scale1 = SQRT2/dist1, scale2 = SQRT2/dist2;
  
  T1[0][0]=scale1 ; T1[0][1]=0.f    ; T1[0][2]= -scale1*v2Centroid1[0];
  T1[1][0]=0.f    ; T1[1][1]=scale1 ; T1[1][2]= -scale1*v2Centroid1[1];
  T1[2][0]=0.f    ; T1[2][1]=0.f    ; T1[2][2]= 1.f;

  T2[0][0]=scale2 ; T2[0][1]=0.f    ; T2[0][2]= -scale2*v2Centroid2[0];
  T2[1][0]=0.f    ; T2[1][1]=scale2 ; T2[1][2]= -scale2*v2Centroid2[1];
  T2[2][0]=0.f    ; T2[2][1]=0.f    ; T2[2][2]= 1.f;
  for( int i = 0; i < no_matchedPairs; ++i )
  {
    MatchedPair & normal = vNormalizeMatchedPairs[i];
    const MatchedPair & org = vMatchedPairs[i];
    normal.v2FstView[0] = scale1*org.v2FstView[0] + T1[0][2];
    normal.v2FstView[1] = scale1*org.v2FstView[1] + T1[1][2];
    normal.v2SndView[0] = scale2*org.v2SndView[0] + T2[0][2];
    normal.v2SndView[1] = scale2*org.v2SndView[1] + T2[1][2];
  }
}

int Homography::estimateSE3(const TooN::Matrix<3> & H, std::vector<int> & vInlierIdx, const std::vector<MatchedPair> & vMatchedPairs)
{
  mvMatchedPairsInCam.resize(5*vInlierIdx.size());
  for( unsigned int idx = 0; idx < vInlierIdx.size(); ++idx )
  {
    for( int i = vInlierIdx[idx], j = idx; i < vInlierIdx[idx]+5; ++ i, ++ j )
    {
//      mvMatchedPairsInCam[j].v2FstView = gUndistortedCameraModel.unproject( vMatchedPairs[i].v2FstView );
//      mvMatchedPairsInCam[j].v2SndView = gUndistortedCameraModel.unproject( vMatchedPairs[i].v2SndView );
    }
  }

  TooN::Matrix<3> Hcam = mInvK*H*mK;
  std::vector<HomographySolution> vHomographySolution;
  TooN::SVD<3> svd(Hcam);
  TooN::Vector<3> v3Diagonal = svd.get_diagonal();
  REAL_TYPE d1 = fabs( v3Diagonal[0] );
  REAL_TYPE d2 = fabs( v3Diagonal[1] );
  REAL_TYPE d3 = fabs( v3Diagonal[2] );

  if( d1 == d2 || d2 == d3 || d1 == d3 )
  {
    std::cout <<"Not implement this case!" << std::endl;
    return false;
  }

  TooN::Matrix<3> U = svd.get_U();
  TooN::Matrix<3> V = svd.get_VT().T();

  REAL_TYPE s = Matrix3::det(U)*Matrix3::det(V);
  REAL_TYPE dPrime = d2;
  REAL_TYPE x1, x2, x3; //eq(12)
  x1 = sqrt( (d1*d1 - d2*d2) / (d1*d1 - d3*d3) );
  x2 = 0;
  x3 = sqrt( (d2*d2 - d3*d3) / (d1*d1 - d3*d3) );
  REAL_TYPE e1[4] = {1.f, -1.f, 1.f, -1.f};
  REAL_TYPE e3[4] = {1.f, 1.f, -1.f, -1.f};

  TooN::Vector<3> v3np;
  HomographySolution solution;

  solution.d = s * dPrime;
  REAL_TYPE invD2 = 1.f/d2;
  // d' > 0
  for( int signs = 0; signs < 4; ++signs ) // eq(13,14)
  {
    solution.Rp = TooN::Identity;
    REAL_TYPE dSinTheta = (d1-d3)*x1*x3*e1[signs]*e3[signs]*invD2;
    REAL_TYPE dCosTheta = (d1*x3*x3+d3*x1*x1)*invD2;
    solution.Rp[0][0] = dCosTheta; solution.Rp[0][2] = -dSinTheta;
    solution.Rp[2][0] = dSinTheta; solution.Rp[2][2] = dCosTheta;

    solution.tp = TooN::makeVector( (d1-d3)*x1*e1[signs], 0.0, -(d1-d3)*x3*e3[signs] );
    v3np = TooN::makeVector( x1*e1[signs], x2, x3*e3[signs] );
    solution.n = V*v3np; // eq(8)
    vHomographySolution.push_back(solution);
  }
  // d' < 0
  solution.d = -s*dPrime;
  for( int signs = 0; signs < 4; ++signs )  // eq(15,16)
  {
    solution.Rp = -1*TooN::Identity;
    REAL_TYPE dSinPhi = (d1+d3)*x1*x3*e1[signs]*e3[signs]*invD2;
    REAL_TYPE dCosPhi = (d3*x1*x1-d1*x3*x3)*invD2;
    solution.Rp[0][0] = dCosPhi; solution.Rp[0][2] = dSinPhi;
    solution.Rp[2][0] = dSinPhi; solution.Rp[2][2] = -dCosPhi;

    solution.tp = TooN::makeVector( (d1+d3)*x1*e1[signs], 0.0, (d1+d3)*x3*e3[signs] );
    v3np = TooN::makeVector( x1*e1[signs], x2, x3*e3[signs] );
    solution.n = V*v3np; // eq(8)
    vHomographySolution.push_back(solution);
  }
  for( unsigned int i = 0; i < vHomographySolution.size(); ++i ) // eq(8)
  {
    vHomographySolution[i].se3First2Second.get_rotation() = s*U*vHomographySolution[i].Rp*V.T();
    vHomographySolution[i].se3First2Second.get_translation() = U*vHomographySolution[i].tp;
  }

  return selectSE3(vHomographySolution);
}

int Homography::estimateSE3(const TooN::Matrix<3> & H)
{
  mvMatchedPairsInCam.resize(mvMatchedPairs.size());
  for( unsigned int i = 0; i < mvMatchedPairsInCam.size(); ++i )
  {
//    mvMatchedPairsInCam[i].v2FstView = gUndistortedCameraModel.unproject( mvMatchedPairs[i].v2FstView );
//    mvMatchedPairsInCam[i].v2SndView = gUndistortedCameraModel.unproject( mvMatchedPairs[i].v2SndView );
  }

  TooN::Matrix<3> Hcam = mInvK*H*mK;
  std::vector<HomographySolution> vHomographySolution;
  TooN::SVD<3> svd(Hcam);
  TooN::Vector<3> v3Diagonal = svd.get_diagonal();
  REAL_TYPE d1 = fabs( v3Diagonal[0] );
  REAL_TYPE d2 = fabs( v3Diagonal[1] );
  REAL_TYPE d3 = fabs( v3Diagonal[2] );

  if( d1 == d2 || d2 == d3 || d1 == d3 )
  {
    std::cout <<"Not implement this case!" << std::endl;
    return false;
  }

  TooN::Matrix<3> U = svd.get_U();
  TooN::Matrix<3> V = svd.get_VT().T();

  REAL_TYPE s = Matrix3::det(U)*Matrix3::det(V);
  REAL_TYPE dPrime = d2;
  REAL_TYPE x1, x2, x3; //eq(12)
  x1 = sqrt( (d1*d1 - d2*d2) / (d1*d1 - d3*d3) );
  x2 = 0;
  x3 = sqrt( (d2*d2 - d3*d3) / (d1*d1 - d3*d3) );
  REAL_TYPE e1[4] = {1.f, -1.f, 1.f, -1.f};
  REAL_TYPE e3[4] = {1.f, 1.f, -1.f, -1.f};

  TooN::Vector<3> v3np;
  HomographySolution solution;

  solution.d = s * dPrime;
  REAL_TYPE invD2 = 1.f/d2;
  // d' > 0
  for( int signs = 0; signs < 4; ++signs ) // eq(13,14)
  {
    solution.Rp = TooN::Identity;
    REAL_TYPE dSinTheta = (d1-d3)*x1*x3*e1[signs]*e3[signs]*invD2;
    REAL_TYPE dCosTheta = (d1*x3*x3+d3*x1*x1)*invD2;
    solution.Rp[0][0] = dCosTheta; solution.Rp[0][2] = -dSinTheta;
    solution.Rp[2][0] = dSinTheta; solution.Rp[2][2] = dCosTheta;

    solution.tp = TooN::makeVector( (d1-d3)*x1*e1[signs], 0.0, -(d1-d3)*x3*e3[signs] );
    v3np = TooN::makeVector( x1*e1[signs], x2, x3*e3[signs] );
    solution.n = V*v3np; // eq(8)
    vHomographySolution.push_back(solution);
  }
  // d' < 0
  solution.d = -s*dPrime;
  for( int signs = 0; signs < 4; ++signs )  // eq(15,16)
  {
    solution.Rp = -1*TooN::Identity;
    REAL_TYPE dSinPhi = (d1+d3)*x1*x3*e1[signs]*e3[signs]*invD2;
    REAL_TYPE dCosPhi = (d3*x1*x1-d1*x3*x3)*invD2;
    solution.Rp[0][0] = dCosPhi; solution.Rp[0][2] = dSinPhi;
    solution.Rp[2][0] = dSinPhi; solution.Rp[2][2] = -dCosPhi;

    solution.tp = TooN::makeVector( (d1+d3)*x1*e1[signs], 0.0, (d1+d3)*x3*e3[signs] );
    v3np = TooN::makeVector( x1*e1[signs], x2, x3*e3[signs] );
    solution.n = V*v3np; // eq(8)
    vHomographySolution.push_back(solution);
  }
  for( unsigned int i = 0; i < vHomographySolution.size(); ++i ) // eq(8)
  {
    vHomographySolution[i].se3First2Second.get_rotation() = s*U*vHomographySolution[i].Rp*V.T();
    vHomographySolution[i].se3First2Second.get_translation() = U*vHomographySolution[i].tp;
  }

  return selectSE3(vHomographySolution);
}

int Homography::estimateSE3FromRectangle(const TooN::Matrix<3> & H)
{
  mvMatchedPairsInCam.resize(mvMatchedPairs.size());
  for( unsigned int i = 0; i < mvMatchedPairsInCam.size(); ++i )
  {
//    mvMatchedPairsInCam[i].v2FstView = gUndistortedCameraModel.unproject( mvMatchedPairs[i].v2FstView );
//    mvMatchedPairsInCam[i].v2SndView = gUndistortedCameraModel.unproject( mvMatchedPairs[i].v2SndView );
  }
  TooN::Matrix<3> Hcam = mInvK*H*mK;
  std::vector<HomographySolution> vHomographySolution;
  TooN::SVD<3> svd(Hcam);
//  TooN::SVD<3> svd(H);
  TooN::Vector<3> v3Diagonal = svd.get_diagonal();
  REAL_TYPE d1 = fabs( v3Diagonal[0] );
  REAL_TYPE d2 = fabs( v3Diagonal[1] );
  REAL_TYPE d3 = fabs( v3Diagonal[2] );

  if( d1 == d2 || d2 == d3 || d1 == d3 )
  {
    std::cout <<"Not implement this case!" << std::endl;
    return false;
  }

  TooN::Matrix<3> U = svd.get_U();
  TooN::Matrix<3> V = svd.get_VT().T();

  REAL_TYPE s = Matrix3::det(U)*Matrix3::det(V);
  REAL_TYPE dPrime = d2;
  REAL_TYPE x1, x2, x3; //eq(12)
  x1 = sqrt( (d1*d1 - d2*d2) / (d1*d1 - d3*d3) );
  x2 = 0;
  x3 = sqrt( (d2*d2 - d3*d3) / (d1*d1 - d3*d3) );
  REAL_TYPE e1[4] = {1.f, -1.f, 1.f, -1.f};
  REAL_TYPE e3[4] = {1.f, 1.f, -1.f, -1.f};

  TooN::Vector<3> v3np;
  HomographySolution solution;

  solution.d = s * dPrime;
  REAL_TYPE invD2 = 1.f/d2;
  // d' > 0
  for( int signs = 0; signs < 4; ++signs ) // eq(13,14)
  {
    solution.Rp = TooN::Identity;
    REAL_TYPE dSinTheta = (d1-d3)*x1*x3*e1[signs]*e3[signs]*invD2;
    REAL_TYPE dCosTheta = (d1*x3*x3+d3*x1*x1)*invD2;
    solution.Rp[0][0] = dCosTheta; solution.Rp[0][2] = -dSinTheta;
    solution.Rp[2][0] = dSinTheta; solution.Rp[2][2] = dCosTheta;

    solution.tp = TooN::makeVector( (d1-d3)*x1*e1[signs], 0.0, -(d1-d3)*x3*e3[signs] );
    v3np = TooN::makeVector( x1*e1[signs], x2, x3*e3[signs] );
    solution.n = V*v3np; // eq(8)
    vHomographySolution.push_back(solution);
  }
  // d' < 0
  solution.d = -s*dPrime;
  for( int signs = 0; signs < 4; ++signs ) // eq(15,16)
  {
    solution.Rp = -1*TooN::Identity;
    REAL_TYPE dSinPhi = (d1+d3)*x1*x3*e1[signs]*e3[signs]*invD2;
    REAL_TYPE dCosPhi = (d3*x1*x1-d1*x3*x3)*invD2;
    solution.Rp[0][0] = dCosPhi; solution.Rp[0][2] = dSinPhi;
    solution.Rp[2][0] = dSinPhi; solution.Rp[2][2] = -dCosPhi;

    solution.tp = TooN::makeVector( (d1+d3)*x1*e1[signs], 0.0, (d1+d3)*x3*e3[signs] );
    v3np = TooN::makeVector( x1*e1[signs], x2, x3*e3[signs] );
    solution.n = V*v3np; // eq(8)
    vHomographySolution.push_back(solution);
  }
  for( unsigned int i = 0; i < vHomographySolution.size(); ++i ) // eq(8)
  {
    vHomographySolution[i].se3First2Second.get_rotation() = s*U*vHomographySolution[i].Rp*V.T();
    vHomographySolution[i].se3First2Second.get_translation() = U*vHomographySolution[i].tp;
  }

  return selectSE3FromRectangle(vHomographySolution);
}

bool operator<(const HomographySolution lhs, const HomographySolution rhs)
{
  return lhs.rScore < rhs.rScore;
}

int Homography::selectSE3FromRectangle( std::vector<HomographySolution> & vHomographySolution )
{
  std::vector<HomographySolution> vVisibleHomographySolution;
// z component must be more than 0 because this point is infront of a camera.
  for( unsigned int i = 0; i < vHomographySolution.size(); ++i )
  {
    unsigned int j; 
    for( j = 0; j < mvMatchedPairsInCam.size(); ++j )
    {
      TooN::Vector<3> v3 = unproject( mvMatchedPairsInCam[j].v2FstView );
      HomographySolution & sol = vHomographySolution[i];
      REAL_TYPE z = v3*sol.n / sol.d;
      if( z < 0.f )
        break;
    }
    if( j == mvMatchedPairsInCam.size() )
    {
      HomographySolution & sol = vHomographySolution[i];
      triangulate( sol.se3First2Second );
      TooN::Vector<3> v3a = mvv3WorldCoordinate[1] - mvv3WorldCoordinate[0];
      TooN::Vector<3> v3b = mvv3WorldCoordinate[2] - mvv3WorldCoordinate[1];
      REAL_TYPE dSumError = fabs(v3a*v3b/(TooN::norm(v3a)*TooN::norm(v3b)));
      v3a = mvv3WorldCoordinate[3] - mvv3WorldCoordinate[2];
      v3b = mvv3WorldCoordinate[0] - mvv3WorldCoordinate[3];
      dSumError += fabs(v3a*v3b/(TooN::norm(v3a)*TooN::norm(v3b)));
      vVisibleHomographySolution.push_back( vHomographySolution[i] );
      vVisibleHomographySolution.back().rScore = dSumError;
    }
  }
  if( vVisibleHomographySolution.size() > 0 )
  {
    if( vVisibleHomographySolution.size() == 1 )
    {
      mvse3First2Second[0] = vVisibleHomographySolution[0].se3First2Second;
      return 1;
    }
    std::sort( vVisibleHomographySolution.begin(), vVisibleHomographySolution.end() );
    if( fabs(vVisibleHomographySolution[0].rScore - vVisibleHomographySolution[1].rScore) < 0.0000001 )
    {
      std::cout << "Hello !!!" << std::endl;
      if( TooN::norm( vVisibleHomographySolution[0].se3First2Second.ln() ) < TooN::norm( vVisibleHomographySolution[1].se3First2Second.ln() ) )
      {
        mvse3First2Second[0] = vVisibleHomographySolution[0].se3First2Second;
        mvse3First2Second[1] = vVisibleHomographySolution[1].se3First2Second;
      }
      else
      {
        mvse3First2Second[0] = vVisibleHomographySolution[1].se3First2Second;
        mvse3First2Second[1] = vVisibleHomographySolution[0].se3First2Second;
      }
      return 2;
    }
    else
    {
//      for( int i = 0; i < vVisibleHomographySolution.size(); ++i )
//        std::cout << vVisibleHomographySolution[i].rScore << std::endl;
      mvse3First2Second[0] = vVisibleHomographySolution[0].se3First2Second;
      mvse3First2Second[1] = vVisibleHomographySolution[1].se3First2Second;
      return 2;
    }
  }
  else
  {
    std::cout <<"Invisible Homography!" << std::endl;
    return 0;
  }
  return 0;
}

int Homography::selectSE3( std::vector<HomographySolution> & vHomographySolution )
{
  std::vector<HomographySolution> vVisibleHomographySolution;
// z component must be more than 0 because this point is infront of a camera.
  for( unsigned int i = 0; i < vHomographySolution.size(); ++i )
  {
    unsigned int j; 
    for( j = 0; j < mvMatchedPairsInCam.size(); ++j )
    {
      TooN::Vector<3> v3 = unproject( mvMatchedPairsInCam[j].v2FstView );
      HomographySolution & sol = vHomographySolution[i];
      REAL_TYPE z = v3*sol.n / sol.d;
      if( z < 0.f )
        break;
    }
    if( j == mvMatchedPairsInCam.size() )
      vVisibleHomographySolution.push_back( vHomographySolution[i] );
  }
  if( vVisibleHomographySolution.size() > 0 )
  {
    if( vVisibleHomographySolution.size() == 1 )
    {
      mvse3First2Second[0] = vVisibleHomographySolution[0].se3First2Second;
      return 1;
    }
    for( unsigned int i = 0; i < vVisibleHomographySolution.size(); ++i )
    {
      TooN::SE3<> & se3 = vVisibleHomographySolution[i].se3First2Second;
      TooN::Matrix<3> essentialMatrix;
      for( int j = 0; j < 3; ++j )
        essentialMatrix.T()[j] = se3.get_translation() ^ se3.get_rotation().get_matrix().T()[j];
      REAL_TYPE dSumError = 0.f;
      for( unsigned int k = 0; k < mvMatchedPairsInCam.size(); ++k )
        dSumError += SampsonDistanceError( mvMatchedPairsInCam[k].v2SndView, essentialMatrix, mvMatchedPairsInCam[k].v2FstView );
      vVisibleHomographySolution[i].rScore = dSumError;
    }
    std::sort( vVisibleHomographySolution.begin(), vVisibleHomographySolution.end() );
    if( fabs(vVisibleHomographySolution[0].rScore - vVisibleHomographySolution[1].rScore) < 0.0000001 )
    {
      std::cout <<"Hello " << std::endl;
      if( TooN::norm(vVisibleHomographySolution[0].se3First2Second.ln()) < TooN::norm( vVisibleHomographySolution[1].se3First2Second.ln() ) )
      {
        mvse3First2Second[0] = vVisibleHomographySolution[0].se3First2Second;
        mvse3First2Second[1] = vVisibleHomographySolution[1].se3First2Second;
      }
      else
      {
        mvse3First2Second[0] = vVisibleHomographySolution[1].se3First2Second;
        mvse3First2Second[1] = vVisibleHomographySolution[0].se3First2Second;
      }
      return 2;
    }
    else
    {
//      for( int i = 0; i < vVisibleHomographySolution.size(); ++i )
//        std::cout << "Score ... " << vVisibleHomographySolution[i].rScore << std::endl;
      mvse3First2Second[0] = vVisibleHomographySolution[0].se3First2Second;
      mvse3First2Second[1] = vVisibleHomographySolution[1].se3First2Second;
      return 2;
    }
  }
  else
  {
    std::cout <<"Invisible Homography!" << std::endl;
    return 0;
  }
  return 0;
}

REAL_TYPE Homography::SampsonDistanceError( const TooN::Vector<2> & v2Second, const TooN::Matrix<3> & essentialMatrix, const TooN::Vector<2> & v2First )
{
  TooN::Vector<3> v3Second = unproject(v2Second);
  TooN::Vector<3> v3First = unproject(v2First);
  REAL_TYPE dError = v3Second*essentialMatrix*v3First;

  TooN::Vector<3> ex1 = essentialMatrix*v3First;
  TooN::Vector<3> ex2 = essentialMatrix.T()*v3Second;

  return ( dError*dError / (ex1*ex1+ex2*ex2) );
}

void Homography::triangulate( const TooN::SE3<> & se3First2Second )
{
  mvv3WorldCoordinate.clear();
  mvv3WorldCoordinate.resize( mvMatchedPairsInCam.size() );
  for( unsigned int i = 0; i < mvMatchedPairsInCam.size(); ++i )
  {
    TooN::Vector<3> v3;
    v3 = Triangulate::compute( se3First2Second, mvMatchedPairsInCam[i].v2FstView, mvMatchedPairsInCam[i].v2SndView );
    mvv3WorldCoordinate[i] = v3;
  }
}

void Homography::cal_plane()
{
// generate normal vector
  TooN::Vector<3> v3Plane = ( mvv3WorldCoordinate[1] - mvv3WorldCoordinate[0] ) ^ ( mvv3WorldCoordinate[2] - mvv3WorldCoordinate[1] );
  assert( norm(v3Plane) > 0.000000001 );
  normalize( v3Plane );
  if( v3Plane[2] < 0 )
    v3Plane *= -1;
//  std::cout <<"Homo Plane : " << v3Plane << std::endl;
  mv4Plane.slice(0,3) = v3Plane;
// set mvv3WorldCoordinate[3] as a point on the plane
  mv4Plane[3] = -1.f * ( v3Plane * mvv3WorldCoordinate[3] );
}

REAL_TYPE Homography::homogdist2d_group_of_5( std::vector<int> & vInlierIdx, TooN::Matrix<3> & H, int iNoGroupIdx, std::vector<MatchedPair> & vMatchedPairs, int no_of_groups, REAL_TYPE threshold )
{
  REAL_TYPE rSumSquaredError = 0;
  TooN::Matrix<3> invH = Matrix3::inv(H);
  TooN::Vector<2> vv2Hx1, vv2invHx2;
//  int matchedSize = vMatchedPairs.size();
  for( int idx = 0; idx < 5*no_of_groups; idx += 5 )
  {
//    if( idx == iNoGroupIdx )
//      continue;
    int count = 0;
    REAL_TYPE rSquaredError = 0;
    while( count != 5 )
    {
      int currentIdx = idx + count;
      TooN::Vector<2> & v2FstView = vMatchedPairs[ currentIdx ].v2FstView;
      TooN::Vector<2> & v2SndView = vMatchedPairs[ currentIdx ].v2SndView;
      vv2Hx1 = v2SndView - project( H*unproject( v2FstView ) );
      vv2invHx2 = v2FstView - project( invH*unproject( v2SndView ) );
      REAL_TYPE d2 = vv2Hx1*vv2Hx1 + vv2invHx2*vv2invHx2;
      if( fabs(d2) > threshold )
        break;
      rSquaredError += d2;
      ++count;
    }
    if( count == 5 )
    {
      rSumSquaredError += rSquaredError;
      vInlierIdx.push_back( idx );
    }
  }
  if( vInlierIdx.size() != 0 )
  {
    return rSumSquaredError/vInlierIdx.size();
  }
  else
    return 1000000.0f;
}

void Homography::ransacHomography_group_of_5(TooN::Matrix<3> & H, std::vector<int> & vInlierIdx, std::vector<TooN::Matrix<3> > & vPossibleHomographys, std::vector<MatchedPair> & vMatchedPairs)
{
  assert( vPossibleHomographys.size() <= vMatchedPairs.size()/5 );
#if 1
  int no_of_groups = vPossibleHomographys.size();
  std::vector<TooN::Matrix<3> >::iterator beginIter = vPossibleHomographys.begin();
  std::vector<TooN::Matrix<3> >::iterator endIter = vPossibleHomographys.end();
  std::vector<TooN::Matrix<3> >::iterator bestIter = beginIter;
  unsigned int max_no_inlier = 0;
  int min_error = 1000000.0f;
  for( ; beginIter != endIter; ++beginIter )
  {
    std::vector<int> vInlierIdxDummy;
    REAL_TYPE err = homogdist2d_group_of_5( vInlierIdxDummy, *beginIter, 0, vMatchedPairs, no_of_groups, 1000.f);
    // need to choose H based on whether no. of inlier or sum squred error.
    // choose no of inlier at this time.
    if( max_no_inlier < vInlierIdxDummy.size() )
    {
      max_no_inlier = vInlierIdxDummy.size();
      min_error = err;
      bestIter = beginIter;
      vInlierIdx = vInlierIdxDummy;
    }
    else if( max_no_inlier == vInlierIdxDummy.size() )
    {
      if( min_error > err )
      {
        min_error = err;
        bestIter = beginIter;
        vInlierIdx = vInlierIdxDummy;
      }
    }
  }
// estimate homography from all inliers pairs
  std::vector<MatchedPair> vInlierMatchedPairs(5*vInlierIdx.size());
  for( unsigned int idx = 0; idx < vInlierIdx.size(); ++idx )
  {
    for( int i = vInlierIdx[idx], j = idx; i < vInlierIdx[idx]+5; ++ i, ++ j )
    {
      vInlierMatchedPairs[j] = vMatchedPairs[i];
    }
  }
  estimateHomography( H, vInlierMatchedPairs, -1 );
//  std::cout <<" Support : " << max_no_inlier << std::endl;
//  H = *bestIter; 
  std::cout << " Check HOMO " << std::endl << H << " " <<  *bestIter << std::endl;
#else
  int no_groups = vPossibleHomographys.size();
  int max_no_inlier = 0;
  int max_group = 0;
  int min_error = 1000000.0f;
  for( int i = 0; i < no_groups; ++ i )
  {
    std::vector<int> vInlierIdxDummy;
    REAL_TYPE err = homogdist2d_group_of_5( vInlierIdxDummy, vPossibleHomographys[i], 0, vMatchedPairs, no_of_groups, 1000.f);
    // need to choose H based on whether no. of inlier or sum squred error.
    // choose no of inlier at this time.
    if( max_no_inlier < vInlierIdxDummy.size() )
    {
      max_no_inlier = vInlierIdxDummy.size();
      max_group = i;
      min_error = err;
      vInlierIdx = vInlierIdxDummy;
    }
    else if( max_no_inlier == vInlierIdxDummy.size() )
    {
      if( min_error > err )
      {
        max_group = i;
        min_error = err;
        vInlierIdx = vInlierIdxDummy;
      }
    }
  }
//  std::cout <<"Group : " << max_group << " Support : " << max_no_inlier << std::endl;
  H = vPossibleHomographys[max_group];
#endif  
}

bool Homography::isDegenerate( const std::vector<MatchedPair> & vMatchedPairs )
{// assume the 1st view is already check from the training phase.
  const TooN::Vector<2> & v2_First = vMatchedPairs[0].v2SndView;
  const TooN::Vector<2> & v2_Second = vMatchedPairs[1].v2SndView;
  const TooN::Vector<2> & v2_Third = vMatchedPairs[2].v2SndView;
  const TooN::Vector<2> & v2_Fourth = vMatchedPairs[3].v2SndView;
  return isColinear( v2_First, v2_Second, v2_Third ) ||
          isColinear( v2_First, v2_Second, v2_Fourth ) ||
          isColinear( v2_First, v2_Third, v2_Fourth ) ||
          isColinear( v2_Second, v2_Third, v2_Fourth );
}
