/*
  Computer Vision Group, Department of Computer Science. University of      
  Bristol. This code can not be used or copied without express permission.   
*/

#ifndef CANNYEDGEDETECTOR_H
#define CANNYEDGEDETECTOR_H

#include <TooN/TooN.h>
#include <assert.h>
#include "DefineTypes.h"
#include "LowLevelImageData.h"
#include <string.h>

struct EdgeProperties
{
  CVD::ImageRef irPossibleEdge;
  CVD::Image<REAL_TYPE>::const_iterator pGradMagnitude;
};
/**
 @file CannyEdgeDetector.h
 @brief It contains functions used for finding edge image using Canny Edge Detection.
*/
class CannyEdgeDetector
{
  public:
    CannyEdgeDetector() {};
    ~CannyEdgeDetector() {};
    /**
      @param[out] imEdge is an result edge image.
      @param[in] llimSource is holding all image information such as raw, smoothed and image.
      @param[in] maxThreshold is a threshold for upper bound.
      @param[in] minThreshold is a threshold for lower bound.
    */
    void compute(CVD::Image<CVD::byte> & imEdge, Edges & edges, LowLevelImageData<CVD::byte> & llimSource, unsigned int uiShortEdge = 3, REAL_TYPE maxThreshold = 0.7f, REAL_TYPE minThreshold = 0.4f );
    void non_max_suppression( CVD::Image<CVD::byte> & imEdge, std::vector<EdgeProperties> & vPossibleEdges, const LowLevelImageData<CVD::byte> & llimSource, REAL_TYPE lowThreshold = 0.01f );
    void hysteresis( CVD::Image<CVD::byte> & imEdge, Edges & edges, std::vector<EdgeProperties> & vPossibleEdges, unsigned int uiShortEdge = 3, REAL_TYPE highThreshold = 0.05f );
    void estimateThreshold( TooN::Vector<2> & threshold, const CVD::Image<REAL_TYPE> & magImage, REAL_TYPE maxThreshold = 0.7f, REAL_TYPE minThreshold = 0.4f );
    void follow_edges( CVD::Image<CVD::byte> & imEdge, std::vector<CVD::ImageRef> & virEdges, CVD::ImageRef & irPossibleEdge );
    void thinning( CVD::Image<CVD::byte> & imEdge, std::vector<CVD::ImageRef> & virEdges );
    void link_edges( CVD::Image<CVD::byte> & imEdge, std::vector<CVD::ImageRef> & virLinkedEdges, std::vector<int> & viLinkedEdgesIndexes, std::vector<CVD::ImageRef> & virEdges, unsigned uiShortEdge = 3 );
  private:
    inline void get_eightNeighbours( unsigned char eightNeighbours[], CVD::Image<CVD::byte> & imEdge, CVD::ImageRef & irCentre );
    inline void no_crossing0to1( int & no_cross, int & no_edges, unsigned char eightNeighbours[] );
    inline bool is_endpoint( unsigned char eightNeighbours[] );
    inline bool is_junction( unsigned char eightNeighbours[] );
};

void CannyEdgeDetector::get_eightNeighbours( unsigned char eightNeighbours[], CVD::Image<CVD::byte> & imEdge, CVD::ImageRef & irCentre )
{
  eightNeighbours[0] = ( imEdge[irCentre+CVD::ImageRef(1,0)] == 255 ) ? 1 : 0;
  eightNeighbours[1] = ( imEdge[irCentre+CVD::ImageRef(1,1)] == 255 ) ? 1 : 0;
  eightNeighbours[2] = ( imEdge[irCentre+CVD::ImageRef(0,1)] == 255 ) ? 1 : 0;
  eightNeighbours[3] = ( imEdge[irCentre+CVD::ImageRef(-1,1)] == 255 ) ? 1 : 0;
  eightNeighbours[4] = ( imEdge[irCentre+CVD::ImageRef(-1,0)] == 255 ) ? 1 : 0;
  eightNeighbours[5] = ( imEdge[irCentre+CVD::ImageRef(-1,-1)] == 255 ) ? 1 : 0;
  eightNeighbours[6] = ( imEdge[irCentre+CVD::ImageRef(0,-1)] == 255 ) ? 1 : 0;
  eightNeighbours[7] = ( imEdge[irCentre+CVD::ImageRef(1,-1)] == 255 ) ? 1 : 0;
};

void CannyEdgeDetector::no_crossing0to1( int & no_cross, int & no_edges, unsigned char eightNeighbours[] )
{
  no_cross = no_edges = 0;
  for( int i = 0; i < 7; ++i )
  {
    if( eightNeighbours[i] == 1 )
      ++no_edges;
    else if( eightNeighbours[i+1] == 1 )
      ++no_cross;
  }
  if( eightNeighbours[7] == 0 && eightNeighbours[0] == 1 )
    ++no_cross;
  if( eightNeighbours[7] == 1 )
    ++no_edges;
};

bool CannyEdgeDetector::is_endpoint( unsigned char eightNeighbours[] )
{
  int sr = 0;
  for( int i = 0; i < 7; ++i )
    if( eightNeighbours[i] == 0 && eightNeighbours[i+1] == 1 )
      ++sr;
  if( eightNeighbours[7] == 0 && eightNeighbours[0] == 1 )
    ++sr;
  return (sr == 1);
};

bool CannyEdgeDetector::is_junction( unsigned char eightNeighbours[] )
{
  int sr = 0;
  for( int i = 0; i < 7; ++i )
    if( eightNeighbours[i] == 0 && eightNeighbours[i+1] == 1 )
      ++sr;
  if( eightNeighbours[7] == 0 && eightNeighbours[0] == 1 )
    ++sr;
  return (sr >= 3);
};

#endif
