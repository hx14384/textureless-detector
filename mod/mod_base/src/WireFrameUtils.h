/*
  Computer Vision Group, Department of Computer Science. University of      
  Bristol. This code can not be used or copied without express permission.   
*/

#ifndef WIREFRAMEUTILS_H
#define WIREFRAMEUTILS_H

#include <TooN/TooN.h>
#include <TooN/se3.h>
#include <vector>
#include <list>
#include <fstream>
#include "WireFrame.h"
/**
 @file WireFrameUtils.h
 @brief There are helping functions for selecting a plane from models and save/load models into/from file.
 */
struct WireFrameUtils
{
  public:
    inline static bool selectPlane( ListPolygonOnPlaneIterator & selectedPlaneIterator, WireFrame * pWireFrame, const TooN::SE3<> & se3W2C, const TooN::Vector<2> & v2Pointer );
    inline static bool selectPlane( ListPolygonOnPlaneIterator & selectedPlaneIterator, ListObjectIterator & listObjectIterator, const TooN::SE3<> & se3W2C, const TooN::Vector<2> & v2Pointer );
    inline static bool selectPlane( ListPolygonOnPlaneIterator & selectedPlaneIterator, ListObjectIterator & listObjectIterator, const TooN::SE3<> & se3W2C, const CVD::ImageRef & irPointer );
    inline static bool save_to_file( std::string filename, WireFrame * pWireFrame, const TooN::SE3<> & se3W2C );
    inline static bool load_from_file( std::string filename, WireFrame * pWireFrame, TooN::SE3<> & se3W2C );
  private:
    inline static int isLeft( const TooN::Vector<2> & v2P0, const TooN::Vector<2> & v2P1, const TooN::Vector<2> & v2P2 );
    inline static int wn_PnPoly( const TooN::Vector<2> & v2P, const std::vector<TooN::Vector<2> > & vv2V );
}
;

bool WireFrameUtils::selectPlane( ListPolygonOnPlaneIterator & selectedPlaneIterator, WireFrame * pWireFrame, const TooN::SE3<> & se3W2C, const TooN::Vector<2> & v2Pointer )
{
  TooN::Vector<2> v2ImagePose;
  std::vector<TooN::Vector<2> > vv2ImagePoses;
  bool result = false;
  ListObjectIterator listObjectIterator = pWireFrame->mListObject.begin();
  ListObjectIterator listObjectEndIterator = pWireFrame->mListObject.end();
  ListPolygonOnPlaneIterator listPolygonOnPlaneIterator;
  for( ; listObjectIterator != listObjectEndIterator; ++listObjectIterator ) // check for every object
  {
    listPolygonOnPlaneIterator = (*listObjectIterator).listPolygonOnPlane.begin();
    ListPolygonOnPlaneIterator listPolygonOnPlaneEndIterator = (*listObjectIterator).listPolygonOnPlane.end();
    for( ; listPolygonOnPlaneIterator != listPolygonOnPlaneEndIterator; ++ listPolygonOnPlaneIterator ) // check for every planes
    {
      vv2ImagePoses.clear();
      TooN::Vector<3> v3PlaneNormal = (*listPolygonOnPlaneIterator).v4Plane.slice(0,3);
      TooN::Vector<3> v3RayPlaneToCam = se3W2C.inverse().get_translation() - (*listPolygonOnPlaneIterator).v3Centroid;
      normalize( v3RayPlaneToCam );
      if( v3RayPlaneToCam * v3PlaneNormal <= 0.02 )
        continue;
      ListListPolygonOnPlaneVertexIteratorIterator listListPolygonOnPlaneVertexIteratorIterator = (*listPolygonOnPlaneIterator).controlListListPolygonOnPlaneVertexIterator.begin();
      ListListPolygonOnPlaneVertexIteratorIterator listListPolygonOnPlaneVertexEndIteratorIterator = (*listPolygonOnPlaneIterator).controlListListPolygonOnPlaneVertexIterator.end();
      for( ; listListPolygonOnPlaneVertexIteratorIterator != listListPolygonOnPlaneVertexEndIteratorIterator; ++ listListPolygonOnPlaneVertexIteratorIterator )
      {
        v2ImagePose = gUndistortedCameraModel.project( project( se3W2C*(*listListPolygonOnPlaneVertexIteratorIterator)->v3WorldCoordinate ) );
        vv2ImagePoses.push_back( v2ImagePose );
      }
      if( vv2ImagePoses.size() <= 2 )
        continue;
      if( vv2ImagePoses.front() != vv2ImagePoses.back() )
        vv2ImagePoses.push_back( vv2ImagePoses.front() );
      if( wn_PnPoly( v2Pointer, vv2ImagePoses ) != 0 )
      {
        result = true;
        selectedPlaneIterator = listPolygonOnPlaneIterator;
        break;
      }
    }
    if( result )
      break;
  }
  return result;
}
;

bool WireFrameUtils::selectPlane( ListPolygonOnPlaneIterator & selectedPlaneIterator, ListObjectIterator & listObjectIterator, const TooN::SE3<> & se3W2C, const TooN::Vector<2> & v2Pointer )
{
  TooN::Vector<2> v2ImagePose;
  std::vector<TooN::Vector<2> > vv2ImagePoses;
  bool result = false;
  ListPolygonOnPlaneIterator listPolygonOnPlaneIterator = (*listObjectIterator).listPolygonOnPlane.begin();
  ListPolygonOnPlaneIterator listPolygonOnPlaneEndIterator = (*listObjectIterator).listPolygonOnPlane.end();
  for( ; listPolygonOnPlaneIterator != listPolygonOnPlaneEndIterator; ++ listPolygonOnPlaneIterator ) // check for every planes
  {
    vv2ImagePoses.clear();
    TooN::Vector<3> v3PlaneNormal = (*listPolygonOnPlaneIterator).v4Plane.slice(0,3);
    TooN::Vector<3> v3RayPlaneToCam = se3W2C.inverse().get_translation() - (*listPolygonOnPlaneIterator).v3Centroid;
    normalize( v3RayPlaneToCam );
    if( v3RayPlaneToCam * v3PlaneNormal <= 0.02 )
      continue;
    ListListPolygonOnPlaneVertexIteratorIterator listListPolygonOnPlaneVertexIteratorIterator = (*listPolygonOnPlaneIterator).controlListListPolygonOnPlaneVertexIterator.begin();
    ListListPolygonOnPlaneVertexIteratorIterator listListPolygonOnPlaneVertexEndIteratorIterator = (*listPolygonOnPlaneIterator).controlListListPolygonOnPlaneVertexIterator.end();
    for( ; listListPolygonOnPlaneVertexIteratorIterator != listListPolygonOnPlaneVertexEndIteratorIterator; ++ listListPolygonOnPlaneVertexIteratorIterator )
    {
      v2ImagePose = gUndistortedCameraModel.project( project( se3W2C*(*listListPolygonOnPlaneVertexIteratorIterator)->v3WorldCoordinate ) );
      vv2ImagePoses.push_back( v2ImagePose );
    }
    if( vv2ImagePoses.size() <= 2 )
      continue;
    if( vv2ImagePoses.front() != vv2ImagePoses.back() )
      vv2ImagePoses.push_back( vv2ImagePoses.front() );
    if( wn_PnPoly( v2Pointer, vv2ImagePoses ) != 0 )
    {
      result = true;
      selectedPlaneIterator = listPolygonOnPlaneIterator;
      break;
    }
  }
  return result;
}
;

bool WireFrameUtils::selectPlane( ListPolygonOnPlaneIterator & selectedPlaneIterator, ListObjectIterator & listObjectIterator, const TooN::SE3<> & se3W2C, const CVD::ImageRef & irPointer )
{
  TooN::Vector<2> v2 = CVD::vec(irPointer);
  return selectPlane( selectedPlaneIterator, listObjectIterator, se3W2C, v2 ); 
}
;

// Copyright 2001, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.
// wn_PnPoly(): winding number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  wn = the winding number (=0 only if P is outside V[])
inline int WireFrameUtils::wn_PnPoly( const TooN::Vector<2> & v2P, const std::vector<TooN::Vector<2> > & vv2V )
{
  int wn = 0; // the winding number counter

  // loop through all edges of the polygon
  // edge from V[i] to V[i+1]
  for( int i = 0; i < vv2V.size() - 1; ++i )
  {
    if( vv2V[i][1] <= v2P[1] ) // start y <= P.y
    {
      if( vv2V[i + 1][1] > v2P[1] ) // an upward crossing
      // P left of edge
        if( isLeft( vv2V[i], vv2V[i + 1], v2P ) > 0 ) ++wn; // have a valid up intersect
    }
    else // start y > P.y (no test needed)
    {
      // a downward crossing
      if( vv2V[i + 1][1] <= v2P[1] )
      // P right of edge
        if( isLeft( vv2V[i], vv2V[i + 1], v2P ) < 0 ) --wn; // have a valid down intersect
    }
  }
  return wn;
}
;

// Copyright 2001, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.

// isLeft(): tests if a point is Left|On|Right of an infinite line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 and P1
//            =0 for P2 on the line
//            <0 for P2 right of the line
//    See: the January 2001 Algorithm "Area of 2D and 3D Triangles and Polygons"
//===================================================================
inline int WireFrameUtils::isLeft( const TooN::Vector<2>& v2P0, const TooN::Vector<2>& v2P1, const TooN::Vector<2>& v2P2 )
{
  return ( ( v2P1[0] - v2P0[0] ) * ( v2P2[1] - v2P0[1] ) - ( v2P2[0] - v2P0[0] ) * ( v2P1[1] - v2P0[1] ) );
}
;

inline bool WireFrameUtils::save_to_file( std::string filename, WireFrame * pWireFrame, const TooN::SE3<> & se3W2C )
{
  std::ofstream ofModels(filename.c_str());
  ofModels << se3W2C << std::endl;
  ofModels << pWireFrame->miNoObject << std::endl;
  ListObject & mListObject = pWireFrame->mListObject;
  for( ReverseListObjectIterator listObjectIterator = mListObject.rbegin(); listObjectIterator != mListObject.rend(); ++listObjectIterator )
  {
    Object & object = (*listObjectIterator);
// Free Vertex
    ofModels << object.listFreeVertex.size() << std::endl;
    for( ReverseListFreeVertexIterator itr = object.listFreeVertex.rbegin(); itr != object.listFreeVertex.rend(); ++itr )
    {
      ofModels << (*itr).bTracked << " " << (*itr).iID << std::endl; 
      ofModels << (*itr).v3CameraCoordinate << (*itr).v3WorldCoordinate << (*itr).v2Image_position << (*itr).v3Normal << std::endl;
    }
// Polygon Vertex
    ofModels << object.listPolygonVertex.size() << std::endl;
    for( ReverseListPolygonVertexIterator itr = object.listPolygonVertex.rbegin(); itr != object.listPolygonVertex.rend(); ++itr )
    {
      ofModels << (*itr).bTracked << " " << (*itr).iID << std::endl;
      ofModels << (*itr).v3CameraCoordinate << (*itr).v3WorldCoordinate << (*itr).v2Image_position << (*itr).v3Normal << (*itr).v3WorldCoordinate_plus_direction << std::endl;
      ofModels << (*itr).bAuto << " " << (*itr).edgeFeature.ucBits << std::endl;
    }
// PolygonOnPlane Vertex
    ofModels << object.listPolygonOnPlaneVertex.size() << std::endl;
    for( ReverseListPolygonOnPlaneVertexIterator itr = object.listPolygonOnPlaneVertex.rbegin(); itr != object.listPolygonOnPlaneVertex.rend(); ++itr )
    {
      ofModels << (*itr).bTracked << " " << (*itr).iID << std::endl;
      ofModels << (*itr).v3CameraCoordinate << (*itr).v3WorldCoordinate << (*itr).v2Image_position << (*itr).v3Normal << (*itr).v3WorldCoordinate_plus_direction << std::endl;
      ofModels << (*itr).bAuto << " " << (*itr).edgeFeature.ucBits << std::endl;
    }
// Boundary Vertex
    ofModels << object.listBoundaryVertex.size() << std::endl;
    for( ReverseListPolygonOnPlaneVertexIterator itr = object.listBoundaryVertex.rbegin(); itr != object.listBoundaryVertex.rend(); ++itr )
    {
      ofModels << (*itr).bTracked << " " << (*itr).iID << std::endl;
      ofModels << (*itr).v3CameraCoordinate << (*itr).v3WorldCoordinate << (*itr).v2Image_position << (*itr).v3Normal << (*itr).v3WorldCoordinate_plus_direction << std::endl;
      ofModels << (*itr).bAuto << " " << (*itr).edgeFeature.ucBits << std::endl;
    }
// Polygon
    ofModels << object.listPolygon.size() << std::endl;
    for( ReverseListPolygonIterator itr = object.listPolygon.rbegin(); itr != object.listPolygon.rend(); ++itr )
    {
      ofModels << (*itr).bTracked << " " << (*itr).iID << std::endl; 
      ofModels << (*itr).iNoVertex << std::endl;
      for( ReverseListListPolygonVertexIteratorIterator itrVert = (*itr).listListPolygonVertexIterator.rbegin(); itrVert
        != (*itr).listListPolygonVertexIterator.rend(); ++itrVert )
      {
        ofModels << (*(*itrVert)).iID << " ";
      }
      ofModels << std::endl;
      ofModels << (*itr).iNoControlVertex << std::endl;
      for( ReverseListListPolygonVertexIteratorIterator itrVert = (*itr).controlListListPolygonVertexIterator.rbegin(); itrVert
        != (*itr).controlListListPolygonVertexIterator.rend(); ++itrVert )
      {
        ofModels << (*(*itrVert)).iID << " ";  
      }
      ofModels << std::endl;
    }
// PolygonOnPlane
    ofModels << object.listPolygonOnPlane.size() << std::endl;
    for( ReverseListPolygonOnPlaneIterator itr = object.listPolygonOnPlane.rbegin(); itr != object.listPolygonOnPlane.rend(); ++itr )
    {
      ofModels << (*itr).bTracked << " " << (*itr).iID << std::endl; 
      ofModels << (*itr).type << " " << (*itr).v4Plane << (*itr).v3Centroid << std::endl;
      ofModels << (*itr).iNoVertex << std::endl;
      for( ReverseListListPolygonOnPlaneVertexIteratorIterator itrVert = (*itr).listListPolygonOnPlaneVertexIterator.rbegin(); itrVert
        != (*itr).listListPolygonOnPlaneVertexIterator.rend(); ++itrVert )
      {
        ofModels << (*(*itrVert)).iID << " ";
      }
      ofModels << std::endl;
      ofModels << (*itr).iNoControlVertex << std::endl;
      for( ReverseListListPolygonOnPlaneVertexIteratorIterator itrVert = (*itr).controlListListPolygonOnPlaneVertexIterator.rbegin(); itrVert
        != (*itr).controlListListPolygonOnPlaneVertexIterator.rend(); ++itrVert )
      {
        ofModels << (*(*itrVert)).iID << " ";
      }
      ofModels << std::endl;
      ofModels << (*itr).iNoBoundaryVertex << std::endl;
      for( ReverseListListPolygonOnPlaneVertexIteratorIterator itrVert = (*itr).listListBoundaryVertexIterator.rbegin(); itrVert
        != (*itr).listListBoundaryVertexIterator.rend(); ++itrVert )
      {
        ofModels << (*(*itrVert)).iID << " ";
      }
      ofModels << std::endl;
    }
  }
  ofModels.close();
  return true;
}
;

inline bool WireFrameUtils::load_from_file( std::string filename, WireFrame * pWireFrame, TooN::SE3<> & se3W2C )
{
  std::ifstream ifModels(filename.c_str());
  pWireFrame->empty();
  int iNoOfObject;
  ifModels >> se3W2C;
  ifModels >> iNoOfObject;
  std::vector<ListPolygonVertexIterator> vlPolygonVertexIterator;
  std::vector<ListPolygonOnPlaneVertexIterator> vlPolygonOnPlaneVertexIterator;
  std::vector<ListPolygonOnPlaneVertexIterator> vlBoundaryVertexIterator;
  for( ; iNoOfObject > 0; --iNoOfObject )
  {
    pWireFrame->createObject();
    ListObjectIterator & mActiveObject = pWireFrame->mActiveObject;
// Free Vertex    
    int iNoOfFreeVertex;
    ifModels >> iNoOfFreeVertex;
    pWireFrame->miNoFreeVertex += iNoOfFreeVertex;
    for(; iNoOfFreeVertex > 0; --iNoOfFreeVertex )
    {
      FreeVertex fv;
      ifModels >> fv.bTracked; ifModels >> fv.iID; 
      ifModels >> fv.v3CameraCoordinate; ifModels >> fv.v3WorldCoordinate; ifModels >> fv.v2Image_position; ifModels >> fv.v3Normal;
      mActiveObject->listFreeVertex.push_front(fv);
    }
// Polygon Vertex
    int iNoOfPolygonVertex;
    ifModels >> iNoOfPolygonVertex;
    pWireFrame->miNoPolygonVertex += iNoOfPolygonVertex;
    vlPolygonVertexIterator.resize(pWireFrame->miNoPolygonVertex);
    for(; iNoOfPolygonVertex > 0; --iNoOfPolygonVertex )
    {
      PolygonVertex pgv;
      ifModels >> pgv.bTracked; ifModels >> pgv.iID;
      ifModels >> pgv.v3CameraCoordinate; ifModels >> pgv.v3WorldCoordinate; ifModels >> pgv.v2Image_position; ifModels >> pgv.v3Normal; ifModels >> pgv.v3WorldCoordinate_plus_direction;
      ifModels >> pgv.bAuto; ifModels >> pgv.edgeFeature.ucBits;
      mActiveObject->listPolygonVertex.push_front(pgv);
      vlPolygonVertexIterator[pgv.iID] = mActiveObject->listPolygonVertex.begin();
    }
// Polygon on Plane Vertex
    int iNoOfPolygonOnPlaneVertex;
    ifModels >> iNoOfPolygonOnPlaneVertex;
    pWireFrame->miNoPolygonOnPlaneVertex += iNoOfPolygonOnPlaneVertex;
    vlPolygonOnPlaneVertexIterator.resize(pWireFrame->miNoPolygonOnPlaneVertex);
    for( ; iNoOfPolygonOnPlaneVertex > 0; --iNoOfPolygonOnPlaneVertex )
    {
      PolygonOnPlaneVertex pgv;
      ifModels >> pgv.bTracked; ifModels >> pgv.iID;
      ifModels >> pgv.v3CameraCoordinate; ifModels >> pgv.v3WorldCoordinate; ifModels >> pgv.v2Image_position; ifModels >> pgv.v3Normal; ifModels >> pgv.v3WorldCoordinate_plus_direction;
      ifModels >> pgv.bAuto; ifModels >> pgv.edgeFeature.ucBits;
      mActiveObject->listPolygonOnPlaneVertex.push_front(pgv);
      vlPolygonOnPlaneVertexIterator[pgv.iID] = mActiveObject->listPolygonOnPlaneVertex.begin();
    }
// Boundary Vertex
    int iNoOfBoundaryVertex;
    ifModels >> iNoOfBoundaryVertex;
    pWireFrame->miNoBoundaryVertex += iNoOfBoundaryVertex;
    vlBoundaryVertexIterator.resize(pWireFrame->miNoBoundaryVertex);
    for( ; iNoOfBoundaryVertex > 0; --iNoOfBoundaryVertex )
    {
      PolygonOnPlaneVertex pgv;
      ifModels >> pgv.bTracked; ifModels >> pgv.iID;
      ifModels >> pgv.v3CameraCoordinate; ifModels >> pgv.v3WorldCoordinate; ifModels >> pgv.v2Image_position; ifModels >> pgv.v3Normal; ifModels >> pgv.v3WorldCoordinate_plus_direction;
      ifModels >> pgv.bAuto; ifModels >> pgv.edgeFeature.ucBits;
      mActiveObject->listBoundaryVertex.push_front(pgv);
      vlBoundaryVertexIterator[pgv.iID] = mActiveObject->listBoundaryVertex.begin();
    }
// Polygon
    int iNoOfPolygon;
    ifModels >> iNoOfPolygon;
    for( ; iNoOfPolygon > 0; --iNoOfPolygon )
    {
      pWireFrame->createPolygon();
      ifModels >> pWireFrame->mActivePolygon->bTracked; ifModels >> pWireFrame->mActivePolygon->iID;
      ifModels >> pWireFrame->mActivePolygon->iNoVertex;
      for( int iNoVertex = pWireFrame->mActivePolygon->iNoVertex; iNoVertex > 0; --iNoVertex )
      {
        int iVertexID;
        ifModels >> iVertexID;
        pWireFrame->mActivePolygon->listListPolygonVertexIterator.push_front( vlPolygonVertexIterator[iVertexID] );
        vlPolygonVertexIterator[iVertexID]->listListPolygonIterator.push_front( pWireFrame->mActivePolygon );
      }
      ifModels >> pWireFrame->mActivePolygon->iNoControlVertex;
      for( int iNoControlVertex = pWireFrame->mActivePolygon->iNoControlVertex; iNoControlVertex > 0; --iNoControlVertex )
      {
        int iVertexID;
        ifModels >> iVertexID;
        pWireFrame->mActivePolygon->controlListListPolygonVertexIterator.push_front( vlPolygonVertexIterator[iVertexID] );
      }
    }
// Polygon On Plane
    int iNoOfPolygonOnPlane;
    ifModels >> iNoOfPolygonOnPlane;
    for(; iNoOfPolygonOnPlane > 0; --iNoOfPolygonOnPlane )
    {
      pWireFrame->createPolygonOnPlane();
      ifModels >> pWireFrame->mActivePolygonOnPlane->bTracked; ifModels >> pWireFrame->mActivePolygonOnPlane->iID;
      int type;
      ifModels >> type;
      switch( type )
      {
        case 0:
          pWireFrame->mActivePolygonOnPlane->type = PolygonOnPlaneType::NONE; 
          break;
        case 1:
          pWireFrame->mActivePolygonOnPlane->type = PolygonOnPlaneType::RECTANGLE;
          break;
        case 2:
          pWireFrame->mActivePolygonOnPlane->type = PolygonOnPlaneType::CIRCLE;
          break;
        case 3:
          pWireFrame->mActivePolygonOnPlane->type = PolygonOnPlaneType::SQUARE;
          break;
        default:
          break;
      }
      ifModels >> pWireFrame->mActivePolygonOnPlane->v4Plane; ifModels >> pWireFrame->mActivePolygonOnPlane->v3Centroid;
      ifModels >> pWireFrame->mActivePolygonOnPlane->iNoVertex;
      for( int iNoVertex = pWireFrame->mActivePolygonOnPlane->iNoVertex; iNoVertex > 0; -- iNoVertex )
      {
        int iVertexID;
        ifModels >> iVertexID;
        pWireFrame->mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.push_front( vlPolygonOnPlaneVertexIterator[iVertexID] );
        vlPolygonOnPlaneVertexIterator[iVertexID]->listListPolygonOnPlaneIterator.push_front( pWireFrame->mActivePolygonOnPlane );
      }
      ifModels >> pWireFrame->mActivePolygonOnPlane->iNoControlVertex;
      for( int iNoControlVertex = pWireFrame->mActivePolygonOnPlane->iNoControlVertex; iNoControlVertex > 0; -- iNoControlVertex )
      {
        int iVertexID;
        ifModels >> iVertexID;
        pWireFrame->mActivePolygonOnPlane->controlListListPolygonOnPlaneVertexIterator.push_front( vlPolygonOnPlaneVertexIterator[iVertexID] );
      }
      ifModels >> pWireFrame->mActivePolygonOnPlane->iNoBoundaryVertex;
      for( int iNoBoundaryVertex = pWireFrame->mActivePolygonOnPlane->iNoBoundaryVertex; iNoBoundaryVertex > 0; -- iNoBoundaryVertex )
      {
        int iVertexID;
        ifModels >> iVertexID;
        pWireFrame->mActivePolygonOnPlane->listListBoundaryVertexIterator.push_front( vlBoundaryVertexIterator[iVertexID] );
      }
    }
  }
  ifModels.close();
  return true;
}
;
#endif
