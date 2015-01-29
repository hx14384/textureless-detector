/*
  Computer Vision Group, Department of Computer Science. University of      
  Bristol. This code can not be used or copied without express permission.   
*/

#ifndef EDGEDETECTOR_H
#define EDGEDETECTOR_H
#include <TooN/TooN.h>
#include <cvd/vector_image_ref.h>
#include <cvd/image.h>
#include <cvd/byte.h>
#include <vector>
#include "DefineTypes.h"
#include "WireFrame.h"

#define MAX_SEARCH_RANGE 15 
//TODO: MAX_SEARCH_RANGE should be varied depended on uncertainty of a camera movement.
/**
 @file EdgeDetector.h
 @brief EdgeDetector is a special edge detector for the tracker class. It is less complexity compared to Canny edge detector. 
*/
class EdgeDetector
{
  public:
    struct ClosestEdge
    {
      bool bNotFound;
      REAL_TYPE rDistance;
      int iEdgeIndex;
    };
    struct DetectedEdge
    {
      bool bNotFound;
      REAL_TYPE rDistance;
      REAL_TYPE x;
      REAL_TYPE y;
    };
    struct EdgeData
    {
      TooN::Vector<6> v6fi;
      unsigned char no_candidate;
      TooN::Vector<2> v2ExpectedImagePoint;
      TooN::Vector<2> v2Normal;
      std::vector<REAL_TYPE> vrDistance;
      std::vector<REAL_TYPE> vrDistanceSquared;
      EdgeData() : no_candidate(0) {};
    };

    EdgeDetector( int iMaxSearchRange = MAX_SEARCH_RANGE, REAL_TYPE fEdgeThreshold = 100.f );
    ~EdgeDetector();
    void mark_edge_point_along_search_path( std::vector<int>& viEdgeIndex,
        std::vector<TooN::Vector<3> > & vv3EdgeStrs, const CVD::Image<REAL_TYPE>& image,
        const TooN::Vector<2>& v2ImagePoint, const TooN::Vector<2>& v2Normal);
    int CompPartialMax( const std::vector<TooN::Vector<3> > & vv3EdgeStrs, int from, int to );
    unsigned int select_edge_point( DetectedEdge& result, unsigned int no_edge,
        const std::vector<int>& viEdgeIndex, const std::vector<TooN::Vector<3> > & vv3EdgeStrs, const TooN::Vector<2> & v2Normal );

    void edgeStrength( TooN::Vector<3>& edgeStr, const CVD::Image<REAL_TYPE>& image, int x, int y );
    void gen_Normal( TooN::Vector<2>& v2Normal, TooN::Vector<2>& v2Direction );
    void search_for_closest_edge_point( ClosestEdge& result, std::vector<TooN::Vector<3> > & vv3EdgeStrs, 
        const CVD::Image<REAL_TYPE> & image, const TooN::Vector<2> & v2ImagePoint, const TooN::Vector<2>& v2Normal );
    void search_for_closest_edge_point_with_pole( ClosestEdge& result, std::vector<TooN::Vector<3> > & vv3EdgeStrs, const EdgeFeature & edgeFeature,
        const CVD::Image<REAL_TYPE> & image, const TooN::Vector<2> & v2ImagePoint, const TooN::Vector<2>& v2Normal );
    void search_for_edge_point( EdgeData& ed, std::vector<TooN::Vector<3> > & vv3EdgeStrs, const CVD::Image<REAL_TYPE> & image,
        const TooN::Vector<2> & v2ImagePoint, const TooN::Vector<2>& v2Normal ); 
    void setEdgeThreshold( REAL_TYPE fEdgeThreshold ) 
    { 
      mfEdgeThreshold = fEdgeThreshold;
    }
    ;
    void setMaxSearchRange( int iMaxSearchRange )
    {
      miMaxSearchRange = iMaxSearchRange;
      mvrPMAX.resize( 2*miMaxSearchRange + 1 );
    }
    ;
  private:
    int miMaxSearchRange;
    REAL_TYPE mfEdgeThreshold;
    std::vector<REAL_TYPE> mvrPMAX;
};
#endif
