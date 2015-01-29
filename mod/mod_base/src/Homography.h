/********************************************************************************/
/**  Computer Vision Group, Department of Computer Science. University of      **/
/**  Bristol. This code can not be used or copied without express permission.  **/ 
/********************************************************************************/

#ifndef HOMOGRAPHYRECONSTRUCTION_H
#define HOMOGRAPHYRECONSTRUCTION_H
#include <TooN/TooN.h>
#include <TooN/se3.h>
#include "DefineTypes.h"
#include "math_utils.h"

struct MatchedPair
{
  TooN::Vector<2> v2FstView;
  TooN::Vector<2> v2SndView;
};

struct HomographySolution
{
  TooN::Matrix<3> Rp;
  TooN::Vector<3> tp;
  TooN::Vector<3> n;
  TooN::SE3<> se3First2Second;
  REAL_TYPE d;
  REAL_TYPE rScore;
};

class Homography
{
  public:
    Homography(); 
    ~Homography();
    bool compute( const std::vector<MatchedPair> & vMatchedPairs );
    bool computeForRectangle( const std::vector<MatchedPair> & vMatchedPairs );
    const std::vector<TooN::SE3<> > & get_se3() { return mvse3First2Second; };
    const TooN::Vector<4> & get_plane() { return mv4Plane; };
    const std::vector<TooN::Vector<3> > & get_world_pose() { return mvv3WorldCoordinate; };
    void estimateHomographyFromRectangle(TooN::Matrix<3> & H);
    void estimateHomography(TooN::Matrix<3> & H);
    int estimateSE3(const TooN::Matrix<3> & H);
    int estimateSE3(const TooN::Matrix<3> & H, std::vector<int> & vInlierIdx, const std::vector<MatchedPair> & vMatchedPairs );
    int estimateSE3FromRectangle( const TooN::Matrix<3> & H );

    bool estimateHomography(TooN::Matrix<3> & H, const std::vector<MatchedPair> & vMatchedPairs, int size );
    bool estimateHomography_group_of_5(TooN::Matrix<3> & H, const std::vector<MatchedPair> & vMatchedPairs );
    double estimateRotationAngle(double & stdDev, std::vector<MatchedPair> & vMatchedPairs);
    void ransacHomography_group_of_5(TooN::Matrix<3> & H, std::vector<int> & vInlierIdx, std::vector<TooN::Matrix<3> > & vPossibleHomographys, std::vector<MatchedPair> & vMatchedPairs);
    inline bool checkNiceHomography(const TooN::Matrix<3> & H)
    {
#if 1
      REAL_TYPE D = H(0,0)*H(1,1)-H(1,0)*H(0,1);
      if( D < 0 ) return false;
      D = sqrt(D);
      if( D > 5 ) return false;
      if( D < 0.05 ) return false;
      D = sqrt( H(0,1)*H(0,1)+H(1,1)*H(1,1) );
      if( D > 5 ) return false;
      if( D < 0.05 ) return false;
      D = sqrt( H(2,0)*H(2,0)+H(2,1)*H(2,1));
      if( D > 0.005 ) return false;
      return true;
#else 
      REAL_TYPE N1 = H(0,0)*H(1,1) - H(1,0)*H(0,1);
      if( N1 > 25 || N1 < 0.0025 )
        return false;
      REAL_TYPE N2 = H(0,1)*H(0,1) + H(1,1)*H(1,1);
      if( N2 > 25 || N2 < 0.0025 )
        return false;
      if( H(2,0)*H(2,0) + H(2,1)*H(2,1) > 0.000025 )
        return false;
      return true;
#endif   
    };

    REAL_TYPE homogdist2d_group_of_5( std::vector<int> & inlierIdx, TooN::Matrix<3> & H, int iNoGroupIdx, std::vector<MatchedPair> & vMatchedPairs, int no_of_groups, REAL_TYPE threshold = 50.f );
    void normalise2dpts( std::vector<MatchedPair> & vNormalizeMatchedPairs, TooN::Matrix<3> & T1, TooN::Matrix<3> & T2, const std::vector<MatchedPair> & vMatchedPairs );
    bool isDegenerate( const std::vector<MatchedPair> & vMatchedPairs );
    
    void set_camera_parameters( TooN::Matrix<3> & K )
    {
      mK = K;
      mInvK = Matrix3::inv(K);
    }

  private:
    REAL_TYPE SampsonDistanceError(const TooN::Vector<2> & v2Second, const TooN::Matrix<3> & essentialMatrix, const TooN::Vector<2> & v2First );// 11.4.3 in Multiple View Geometry in Compuer Vision
    int max (int a, int b);
    int min (int a, int b);
    void triangulate( const TooN::SE3<> & se3First2Second );
    void cal_plane();
    int selectSE3( std::vector<HomographySolution> & vHomographySolution );
    int selectSE3FromRectangle( std::vector<HomographySolution> & vHomographySolution );
    std::vector<TooN::SE3<> >mvse3First2Second;
    TooN::Vector<4> mv4Plane;
    std::vector<MatchedPair> mvMatchedPairs;
    std::vector<MatchedPair> mvMatchedPairsInCam;
    std::vector<TooN::Vector<3> > mvv3WorldCoordinate;
    std::vector<MatchedPair> mv5NormalizeMatchedPairs;
    TooN::Matrix<3> mK;
    TooN::Matrix<3> mInvK;
};

#endif
