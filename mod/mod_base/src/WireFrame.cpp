/*
  Computer Vision Group, Department of Computer Science. University of      
  Bristol. This code can not be used or copied without express permission.   
*/

#include "WireFrame.h"
#include <assert.h>
#include <algorithm>
#include "EdgeDetector.h"

WireFrame::WireFrame() :
  miNoObject(0), miNoFreeVertex(0), miNoPolygonVertex(0), miNoPolygonOnPlaneVertex(0), miNoBoundaryVertex(0), miNoPolygon(0), miNoPolygonOnPlane(0)
{
  mActiveObject = mListObject.end();
}

WireFrame::~WireFrame()
{

}

void WireFrame::empty()
{
  mListObject.clear();
  mActiveObject = mListObject.end();
  miNoObject =0;
  miNoFreeVertex =0;
  miNoPolygonVertex =0;
  miNoPolygonOnPlaneVertex =0;
  miNoBoundaryVertex =0;
  miNoPolygon =0;
  miNoPolygonOnPlane =0;
}

void WireFrame::createObject( )
{
  Object ob( miNoObject++ );
  mListObject.push_front( ob );
  mActiveObject = mListObject.begin();
  mActiveFreeVertex = mActiveObject->listFreeVertex.end();
  mActivePolygonVertex = mActiveObject->listPolygonVertex.end();
  mActivePolygonOnPlaneVertex = mActiveObject->listPolygonOnPlaneVertex.end();
  mActivePolygon = mActiveObject->listPolygon.end();
  mActivePolygonOnPlane = mActiveObject->listPolygonOnPlane.end();
}

void WireFrame::addFreeVertex( const TooN::Vector<3> & v3pose, const TooN::Vector<3> & v3Normal )
{
  assert( mActiveObject != mListObject.end() );
  FreeVertex fv;
  fv.iID = miNoFreeVertex++;
  fv.v3Normal = v3Normal;
  fv.v3WorldCoordinate = v3pose;
  addFreeVertex( fv );
}

void WireFrame::addFreeVertex( const FreeVertex & freeVertex )
{
  assert( mActiveObject != mListObject.end() );
  mActiveObject->listFreeVertex.push_front( freeVertex );
  mActiveFreeVertex = mActiveObject->listFreeVertex.begin();
}

void WireFrame::createPolygon()
{
  assert( mActiveObject != mListObject.end() );
  Polygon pg;
  pg.iID = miNoPolygon++;
  mActiveObject->listPolygon.push_front( pg );
  mActivePolygon = mActiveObject->listPolygon.begin();
}

void WireFrame::addPolygonVertex( const TooN::Vector<3> & v3pose, const TooN::Vector<3> & v3Normal, bool bGenerateSamplePoints )
{
  assert( mActiveObject != mListObject.end() && mActivePolygon != mActiveObject->listPolygon.end() );
  PolygonVertex pgv;
  pgv.v3WorldCoordinate = v3pose;
  pgv.v3Normal = v3Normal;
  if( !mActivePolygon->listListPolygonVertexIterator.empty() )
  {
    TooN::Vector<3> v3direction = v3pose - mActivePolygon->listListPolygonVertexIterator.front()->v3WorldCoordinate;
    REAL_TYPE no_segment = 8.0*sqrt( v3direction * v3direction );
    TooN::Vector<3> v3Step = v3direction / no_segment;
    normalize( v3direction );
    pgv.v3WorldCoordinate_plus_direction = v3pose - 0.01f * v3direction;
    if( mActivePolygon->listListPolygonVertexIterator.size() == 1 ) mActivePolygon->listListPolygonVertexIterator.back()->v3WorldCoordinate_plus_direction
        = mActivePolygon->listListPolygonVertexIterator.back()->v3WorldCoordinate + 0.01f * v3direction;
    if( bGenerateSamplePoints )
    {
      // start
      // generate auto PolygonVertex
      TooN::Vector<3> v3AutoPose = mActivePolygon->listListPolygonVertexIterator.front()->v3WorldCoordinate + v3Step;
      while( no_segment > 1.0f )
      {
        PolygonVertex pgvAuto;
        pgvAuto.v3WorldCoordinate = v3AutoPose;
        pgvAuto.v3WorldCoordinate_plus_direction = v3AutoPose - 0.01f * v3direction;
        pgvAuto.v3Normal = v3Normal;
        pgvAuto.bAuto = true;
        pgvAuto.listListPolygonIterator.push_front( mActivePolygon );
        pgvAuto.iID = miNoPolygonVertex++;
        addPolygonVertex( pgvAuto );
        no_segment -= 1.f;
        v3AutoPose += v3Step;
      }
      // end
    }
  }
  pgv.listListPolygonIterator.push_front( mActivePolygon );
  pgv.iID = miNoPolygonVertex++;
  addPolygonVertex( pgv );
}

void WireFrame::addPolygonVertex( const PolygonVertex & polygonVertex )
{
  mActiveObject->listPolygonVertex.push_front( polygonVertex );
  mActivePolygonVertex = mActiveObject->listPolygonVertex.begin();
  mActivePolygonVertex->listListPolygonIterator.push_front( mActivePolygon );
  mActivePolygon->iNoVertex += 1;
  mActivePolygon->listListPolygonVertexIterator.push_front( mActivePolygonVertex );
  if( !polygonVertex.bAuto ) 
  {
    mActivePolygon->iNoControlVertex += 1;
    mActivePolygon->controlListListPolygonVertexIterator.push_front( mActivePolygonVertex );
  }
}

void WireFrame::closePolygon( bool bGenerateSamplePoints )
{
  assert( mActivePolygon->listListPolygonVertexIterator.size() >= 3 );
  if( bGenerateSamplePoints )
  {
    TooN::Vector<3> v3direction = mActivePolygon->listListPolygonVertexIterator.back()->v3WorldCoordinate
        - mActivePolygon->listListPolygonVertexIterator.front()->v3WorldCoordinate;
    REAL_TYPE no_segment = 8.0*sqrt( v3direction * v3direction );
    TooN::Vector<3> v3Step = v3direction / no_segment;
    normalize( v3direction );
    TooN::Vector<3> v3Normal = mActivePolygon->listListPolygonVertexIterator.front()->v3Normal;
    TooN::Vector<3> v3AutoPose = mActivePolygon->listListPolygonVertexIterator.front()->v3WorldCoordinate + v3Step;
    while( no_segment > 1.0f )
    {
      PolygonVertex pgvAuto;
      pgvAuto.v3WorldCoordinate = v3AutoPose;
      pgvAuto.v3WorldCoordinate_plus_direction = v3AutoPose - 0.01f * v3direction;
      pgvAuto.v3Normal = v3Normal;
      pgvAuto.bAuto = true;
      pgvAuto.listListPolygonIterator.push_front( mActivePolygon );
      pgvAuto.iID = miNoPolygonVertex++;
      addPolygonVertex( pgvAuto );
      no_segment -= 1.f;
      v3AutoPose += v3Step;
    }
  }
  mActivePolygonVertex = *mActivePolygon->listListPolygonVertexIterator.begin();
  mActivePolygonVertex->listListPolygonIterator.push_front( mActivePolygon );
  mActivePolygon->iNoVertex += 1;
  mActivePolygon->listListPolygonVertexIterator.push_front( mActivePolygonVertex );
}

void WireFrame::createPolygonOnPlane( PolygonOnPlaneType::Type type )
{
  assert( mActiveObject != mListObject.end() );
  PolygonOnPlane pg;
  pg.iID = miNoPolygonOnPlane++;
  mActiveObject->listPolygonOnPlane.push_front( pg );
  mActivePolygonOnPlane = mActiveObject->listPolygonOnPlane.begin();
  mActivePolygonOnPlane->type = type;
}

void WireFrame::addPolygonOnPlaneVertex( const TooN::Vector<3> & v3pose, bool bGenerateSamplePoints )
{
  assert( mActiveObject != mListObject.end() && mActivePolygonOnPlane != mActiveObject->listPolygonOnPlane.end() );
  PolygonOnPlaneVertex pgv;
  pgv.v3WorldCoordinate = v3pose;
  pgv.v3Normal = mActivePolygonOnPlane->v4Plane.slice( 0, 3 );
  if( !mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.empty() ) // if it is not the first vertex.
  {
    TooN::Vector<3> v3direction = v3pose
        - mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.front()->v3WorldCoordinate;
    REAL_TYPE no_segment = 8.0*sqrt( v3direction * v3direction );
    TooN::Vector<3> v3Step = v3direction / no_segment;
    normalize( v3direction );
    pgv.v3WorldCoordinate_plus_direction = v3pose - 0.01f * v3direction;
    if( mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.size() == 1 ) mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.back()->v3WorldCoordinate_plus_direction
        = mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.back()->v3WorldCoordinate + 0.01f * v3direction;
    if( bGenerateSamplePoints )
    {
      // start
      // generate auto PolygonVertex
      TooN::Vector<3> v3AutoPose = mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.front()->v3WorldCoordinate
          + v3Step;
      while( no_segment > 1.0f )
      {
        PolygonOnPlaneVertex pgvAuto;
        pgvAuto.v3WorldCoordinate = v3AutoPose;
        pgvAuto.v3WorldCoordinate_plus_direction = v3AutoPose - 0.01f * v3direction;
        pgvAuto.v3Normal = mActivePolygonOnPlane->v4Plane.slice( 0, 3 );
        pgvAuto.bAuto = true;
        pgvAuto.iID = miNoPolygonOnPlaneVertex++;
        addPolygonOnPlaneVertex( pgvAuto );
        no_segment -= 1.f;
        v3AutoPose += v3Step;
      }
      // end
    }
  }
  pgv.iID = miNoPolygonOnPlaneVertex++;
  addPolygonOnPlaneVertex( pgv );
}

void WireFrame::addPolygonOnPlaneVertex( const TooN::Vector<3> & v3pose, const TooN::SE3<> & se3W2C, const CVD::Image<
    CVD::byte> & image, bool bGenerateSamplePoints )
{
  assert( mActiveObject != mListObject.end() && mActivePolygonOnPlane != mActiveObject->listPolygonOnPlane.end() );
  PolygonOnPlaneVertex pgv;
  pgv.v3WorldCoordinate = v3pose;
  pgv.v3CameraCoordinate = se3W2C * v3pose;
  pgv.v2Image_position = gQuinticCameraModel.project( project( pgv.v3CameraCoordinate ) );
  pgv.v3Normal = mActivePolygonOnPlane->v4Plane.slice( 0, 3 );
  std::vector<int> viEdgeIndex;
  std::vector<TooN::Vector<3> > vv3EdgeStrs( 2*MAX_SEARCH_RANGE + 1 );
  if( !mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.empty() ) // if it is not the first vertex.
  {
    TooN::Vector<3> v3direction = v3pose
        - mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.front()->v3WorldCoordinate;
    REAL_TYPE no_segment = 8.0*sqrt( v3direction * v3direction );
    TooN::Vector<3> v3Step = v3direction / no_segment;
    normalize( v3direction );
    pgv.v3WorldCoordinate_plus_direction = v3pose - 0.01f * v3direction;
    // determine edge pole
    //    TooN::Vector<2> v2Direction = pgv.v2Image_position - gQuinticCameraModel.project(project(se3W2C*pgv.v3WorldCoordinate_plus_direction));
    //    TooN::Vector<2> v2Normal;
    //    EdgeDetector::gen_Normal(v2Normal, v2Direction);
    //    EdgeDetector::mark_edge_point_along_search( viEdgeIndex, vv3EdgeStrs, image, pgv.v2Image_position, v2Normal, MAX_SEARCH_RANGE, 10 );
    // end determine edge pole.
    if( mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.size() == 1 ) mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.back()->v3WorldCoordinate_plus_direction
        = mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.back()->v3WorldCoordinate + 0.01f * v3direction;
    if( bGenerateSamplePoints )
    {
      // start
      // generate auto PolygonVertex
      TooN::Vector<3> v3AutoPose = mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.front()->v3WorldCoordinate
          + v3Step;
      while( no_segment > 1.0f )
      {
        PolygonOnPlaneVertex pgvAuto;
        pgvAuto.v3WorldCoordinate = v3AutoPose;
        pgvAuto.v3CameraCoordinate = se3W2C * v3AutoPose;
        pgvAuto.v2Image_position = gQuinticCameraModel.project( project( pgvAuto.v3CameraCoordinate ) );
        pgvAuto.v3WorldCoordinate_plus_direction = v3AutoPose - 0.01f * v3direction;
        pgvAuto.v3Normal = mActivePolygonOnPlane->v4Plane.slice( 0, 3 );
        pgvAuto.bAuto = true;
        pgvAuto.iID = miNoPolygonOnPlaneVertex++;
        addPolygonOnPlaneVertex( pgvAuto );
        no_segment -= 1.f;
        v3AutoPose += v3Step;
      }
      // end
    }
  }
  pgv.iID = miNoPolygonOnPlaneVertex++;
  addPolygonOnPlaneVertex( pgv );
}

void WireFrame::addPolygonOnPlaneVertex( const PolygonOnPlaneVertex & polygonOnPlaneVertex )
{
  mActiveObject->listPolygonOnPlaneVertex.push_front( polygonOnPlaneVertex );
  mActivePolygonOnPlaneVertex = mActiveObject->listPolygonOnPlaneVertex.begin();
  mActivePolygonOnPlaneVertex->listListPolygonOnPlaneIterator.push_front( mActivePolygonOnPlane );
  mActivePolygonOnPlane->iNoVertex += 1;
  mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.push_front( mActivePolygonOnPlaneVertex );
  if( !polygonOnPlaneVertex.bAuto )
  {
    int no_control_points = mActivePolygonOnPlane->controlListListPolygonOnPlaneVertexIterator.size();
    mActivePolygonOnPlane->v3Centroid = ( no_control_points * mActivePolygonOnPlane->v3Centroid
        + mActivePolygonOnPlaneVertex->v3WorldCoordinate ) / ( 1.f + no_control_points );
    mActivePolygonOnPlane->iNoControlVertex += 1;
    mActivePolygonOnPlane->controlListListPolygonOnPlaneVertexIterator.push_front( mActivePolygonOnPlaneVertex );
  }
}

void WireFrame::addPolygonOnPlaneVertex( PolygonOnPlaneVertex & polygonOnPlaneVertex, ListObjectIterator & listObjectIterator, ListPolygonOnPlaneIterator & listPolygonOnPlaneIterator )
{
  polygonOnPlaneVertex.iID = miNoPolygonOnPlaneVertex++;
  listObjectIterator->listPolygonOnPlaneVertex.push_front( polygonOnPlaneVertex );
  listObjectIterator->listPolygonOnPlaneVertex.begin()->listListPolygonOnPlaneIterator.push_front( listPolygonOnPlaneIterator );
  listPolygonOnPlaneIterator->iNoVertex += 1;
  listPolygonOnPlaneIterator->listListPolygonOnPlaneVertexIterator.push_front( listObjectIterator->listPolygonOnPlaneVertex.begin() );
}

void WireFrame::addBoundaryVertex( const TooN::Vector<3> & v3pose )
{
  assert( mActiveObject != mListObject.end() && mActivePolygonOnPlane != mActiveObject->listPolygonOnPlane.end() );
  PolygonOnPlaneVertex pgv;
  pgv.v3WorldCoordinate = v3pose;
  pgv.v3Normal = mActivePolygonOnPlane->v4Plane.slice( 0, 3 );
  if( !mActivePolygonOnPlane->listListBoundaryVertexIterator.empty() ) // if it is not the first vertex.
  {
    TooN::Vector<3> v3direction = v3pose
        - mActivePolygonOnPlane->listListBoundaryVertexIterator.front()->v3WorldCoordinate;
    REAL_TYPE no_segment = 4.0*sqrt( v3direction * v3direction );
    TooN::Vector<3> v3Step = v3direction / no_segment;
    normalize( v3direction );
    pgv.v3WorldCoordinate_plus_direction = v3pose - 0.01f * v3direction;
    if( mActivePolygonOnPlane->listListBoundaryVertexIterator.size() == 1 ) mActivePolygonOnPlane->listListBoundaryVertexIterator.back()->v3WorldCoordinate_plus_direction
        = mActivePolygonOnPlane->listListBoundaryVertexIterator.back()->v3WorldCoordinate + 0.01f * v3direction;
    TooN::Vector<3> v3AutoPose = mActivePolygonOnPlane->listListBoundaryVertexIterator.front()->v3WorldCoordinate
        + v3Step;
    while( no_segment > 1.0f )
    {
      PolygonOnPlaneVertex pgvAuto;
      pgvAuto.v3WorldCoordinate = v3AutoPose;
      pgvAuto.v3WorldCoordinate_plus_direction = v3AutoPose - 0.01f * v3direction;
      pgvAuto.v3Normal = mActivePolygonOnPlane->v4Plane.slice( 0, 3 );
      pgvAuto.iID = miNoBoundaryVertex++;
      addBoundaryVertex( pgvAuto );
      no_segment -= 1.f;
      v3AutoPose += v3Step;
    }
  }
  pgv.iID = miNoBoundaryVertex++;
  addBoundaryVertex( pgv );
}

void WireFrame::addBoundaryVertex( const PolygonOnPlaneVertex & polygonOnPlaneVertex )
{
  mActiveObject->listBoundaryVertex.push_front( polygonOnPlaneVertex );
  mActiveBoundaryVertex = mActiveObject->listBoundaryVertex.begin();
  mActiveBoundaryVertex->listListPolygonOnPlaneIterator.push_front( mActivePolygonOnPlane );
  mActivePolygonOnPlane->iNoBoundaryVertex += 1;
  mActivePolygonOnPlane->listListBoundaryVertexIterator.push_front( mActiveBoundaryVertex );
}

void WireFrame::addPolygonOnPlaneVertex( ListPolygonOnPlaneVertexIterator & listPolygonOnPlaneVertexIterator,
    bool bGenerateSamplePoints )
{
  if( bGenerateSamplePoints )
  {
    if( !mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.empty() ) // if it is not the first vertex.
    {
      TooN::Vector<3> v3direction = listPolygonOnPlaneVertexIterator->v3WorldCoordinate
          - mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.front()->v3WorldCoordinate;
      REAL_TYPE no_segment = 8.0*sqrt( v3direction * v3direction );
      TooN::Vector<3> v3Step = v3direction / no_segment;
      normalize( v3direction );
      if( mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.size() == 1 ) mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.back()->v3WorldCoordinate_plus_direction
          = mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.back()->v3WorldCoordinate + 0.01f * v3direction;
      // start
      // generate auto PolygonVertex
      TooN::Vector<3> v3AutoPose =
          mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.front()->v3WorldCoordinate + v3Step;
      while( no_segment > 1.0f )
      {
        PolygonOnPlaneVertex pgvAuto;
        pgvAuto.v3WorldCoordinate = v3AutoPose;
        pgvAuto.v3WorldCoordinate_plus_direction = v3AutoPose - 0.01f * v3direction;
        pgvAuto.v3Normal = mActivePolygonOnPlane->v4Plane.slice( 0, 3 );
        pgvAuto.bAuto = true;
        pgvAuto.iID = miNoPolygonOnPlaneVertex++;
        addPolygonOnPlaneVertex( pgvAuto );
        no_segment -= 1.f;
        v3AutoPose += v3Step;
      }
      // end
    }
  }
  mActivePolygonOnPlaneVertex = listPolygonOnPlaneVertexIterator;
  mActivePolygonOnPlaneVertex->v3Normal = ( mActivePolygonOnPlaneVertex->v3Normal
      + mActivePolygonOnPlane->v4Plane.slice( 0, 3 ) ) * 0.5f; // may not need this
  mActivePolygonOnPlaneVertex->listListPolygonOnPlaneIterator.push_front( mActivePolygonOnPlane );
  mActivePolygonOnPlane->iNoVertex += 1;
  mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.push_front( mActivePolygonOnPlaneVertex );
  if( !( *listPolygonOnPlaneVertexIterator ).bAuto )
  {
    int no_control_points = mActivePolygonOnPlane->controlListListPolygonOnPlaneVertexIterator.size();
    mActivePolygonOnPlane->v3Centroid = ( no_control_points * mActivePolygonOnPlane->v3Centroid
        + mActivePolygonOnPlaneVertex->v3WorldCoordinate ) / ( 1.f + no_control_points );
    mActivePolygonOnPlane->iNoControlVertex += 1;
    mActivePolygonOnPlane->controlListListPolygonOnPlaneVertexIterator.push_front( mActivePolygonOnPlaneVertex );
  }
}

void WireFrame::addPolygonOnPlaneVertex( ListPolygonOnPlaneVertexIterator & listPolygonOnPlaneVertexIterator,
    const TooN::SE3<> & se3W2C, const CVD::Image<CVD::byte> & image, bool bGenerateSamplePoints )
{
  if( !mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.empty() ) // if it is not the first vertex.
  {
    TooN::Vector<3> v3direction = listPolygonOnPlaneVertexIterator->v3WorldCoordinate
        - mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.front()->v3WorldCoordinate;
    REAL_TYPE no_segment = 8.0*sqrt( v3direction * v3direction );
    TooN::Vector<3> v3Step = v3direction / no_segment;
    normalize( v3direction );
    if( mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.size() == 1 ) mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.back()->v3WorldCoordinate_plus_direction
        = mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.back()->v3WorldCoordinate + 0.01f * v3direction;
    if( bGenerateSamplePoints )
    {
      // start
      // generate auto PolygonVertex
      TooN::Vector<3> v3AutoPose = mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.front()->v3WorldCoordinate
          + v3Step;
      while( no_segment > 1.0f )
      {
        PolygonOnPlaneVertex pgvAuto;
        pgvAuto.v3WorldCoordinate = v3AutoPose;
        pgvAuto.v3CameraCoordinate = se3W2C * v3AutoPose;
        pgvAuto.v2Image_position = gQuinticCameraModel.project( project( pgvAuto.v3CameraCoordinate ) );
        pgvAuto.v3WorldCoordinate_plus_direction = v3AutoPose - 0.01f * v3direction;
        pgvAuto.v3Normal = mActivePolygonOnPlane->v4Plane.slice( 0, 3 );
        pgvAuto.bAuto = true;
        pgvAuto.iID = miNoPolygonOnPlaneVertex++;
        addPolygonOnPlaneVertex( pgvAuto );
        no_segment -= 1.f;
        v3AutoPose += v3Step;
      }
    // end
    }
  }
  mActivePolygonOnPlaneVertex = listPolygonOnPlaneVertexIterator;
  mActivePolygonOnPlaneVertex->v3Normal = ( mActivePolygonOnPlaneVertex->v3Normal
      + mActivePolygonOnPlane->v4Plane.slice( 0, 3 ) ) * 0.5f; // may not need this
  mActivePolygonOnPlaneVertex->listListPolygonOnPlaneIterator.push_front( mActivePolygonOnPlane );
  mActivePolygonOnPlane->iNoVertex += 1;
  mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.push_front( mActivePolygonOnPlaneVertex );
  if( !( *listPolygonOnPlaneVertexIterator ).bAuto )
  {
    int no_control_points = mActivePolygonOnPlane->controlListListPolygonOnPlaneVertexIterator.size();
    mActivePolygonOnPlane->v3Centroid = ( no_control_points * mActivePolygonOnPlane->v3Centroid
        + mActivePolygonOnPlaneVertex->v3WorldCoordinate ) / ( 1.f + no_control_points );
    mActivePolygonOnPlane->iNoControlVertex += 1;
    mActivePolygonOnPlane->controlListListPolygonOnPlaneVertexIterator.push_front( mActivePolygonOnPlaneVertex );
  }
}

void WireFrame::addLineToPolygonOnPlane( ListPolygonOnPlaneVertexIterator & firstVertexIterator,
    ListPolygonOnPlaneVertexIterator & secondVertexIterator, bool bGenerateSamplePoints )
{
  ListListPolygonOnPlaneIteratorIterator pPlaneIterIterator =
      firstVertexIterator->listListPolygonOnPlaneIterator.begin();
  ListListPolygonOnPlaneIteratorIterator pPlaneEndIterIterator =
      firstVertexIterator->listListPolygonOnPlaneIterator.end();
  for( ; pPlaneIterIterator != pPlaneEndIterIterator; ++pPlaneIterIterator ) //for each plane that firstVertexIterator belongs to.
  {
    if( ( *( *pPlaneIterIterator ) ).controlListListPolygonOnPlaneVertexIterator.size() < 2 ) continue;
    ListListPolygonOnPlaneVertexIteratorIterator nextVertexIterIterator = find(
        ( *( *pPlaneIterIterator ) ).controlListListPolygonOnPlaneVertexIterator.begin(),
        ( *( *pPlaneIterIterator ) ).controlListListPolygonOnPlaneVertexIterator.end(), firstVertexIterator );
    ListListPolygonOnPlaneVertexIteratorIterator previousVertexIterIter = nextVertexIterIterator;
    if( nextVertexIterIterator != ( *( *pPlaneIterIterator ) ).controlListListPolygonOnPlaneVertexIterator.end() )
    {
      ++nextVertexIterIterator;
      if( *nextVertexIterIterator == secondVertexIterator ) // found in this plane.
      {
        ListListPolygonOnPlaneVertexIteratorIterator firstListListPolygonOnPlaneVertexIteratorIterator = find(
            ( *( *pPlaneIterIterator ) ).listListPolygonOnPlaneVertexIterator.begin(),
            ( *( *pPlaneIterIterator ) ).listListPolygonOnPlaneVertexIterator.end(), firstVertexIterator );
        ++firstListListPolygonOnPlaneVertexIteratorIterator; // not include the firstListListPolygonOnPlaneVertexIteratorIterator
        ListListPolygonOnPlaneVertexIteratorIterator secondListListPolygonOnPlaneVertexIteratorIterator = find(
            firstListListPolygonOnPlaneVertexIteratorIterator,
            ( *( *pPlaneIterIterator ) ).listListPolygonOnPlaneVertexIterator.end(), secondVertexIterator );
        while( firstListListPolygonOnPlaneVertexIteratorIterator != secondListListPolygonOnPlaneVertexIteratorIterator )
        {
          //          std::cout << (*(*firstListListPolygonOnPlaneVertexIteratorIterator)).v3WorldCoordinate << "\t" << (*(*secondListListPolygonOnPlaneVertexIteratorIterator)).v3WorldCoordinate << std::endl;
          addPolygonOnPlaneVertex( *firstListListPolygonOnPlaneVertexIteratorIterator, false );
          ++firstListListPolygonOnPlaneVertexIteratorIterator;
        }
        addPolygonOnPlaneVertex( *secondListListPolygonOnPlaneVertexIteratorIterator, false );
        break;
      }
    }
    if( previousVertexIterIter != ( *( *pPlaneIterIterator ) ).controlListListPolygonOnPlaneVertexIterator.begin() )
    {
      --previousVertexIterIter;
      if( *previousVertexIterIter == secondVertexIterator ) // found in this plane.
      {
        ListListPolygonOnPlaneVertexIteratorIterator secondListListPolygonOnPlaneVertexIteratorIterator = find(
            ( *( *pPlaneIterIterator ) ).listListPolygonOnPlaneVertexIterator.begin(),
            ( *( *pPlaneIterIterator ) ).listListPolygonOnPlaneVertexIterator.end(), secondVertexIterator );
        ListListPolygonOnPlaneVertexIteratorIterator firstListListPolygonOnPlaneVertexIteratorIterator = find(
            secondListListPolygonOnPlaneVertexIteratorIterator,
            ( *( *pPlaneIterIterator ) ).listListPolygonOnPlaneVertexIterator.end(), firstVertexIterator );
        --firstListListPolygonOnPlaneVertexIteratorIterator; // not include the firstListListPolygonOnPlaneVertexIteratorIterator
        while( firstListListPolygonOnPlaneVertexIteratorIterator != secondListListPolygonOnPlaneVertexIteratorIterator )
        {
          addPolygonOnPlaneVertex( *firstListListPolygonOnPlaneVertexIteratorIterator, false );
          --firstListListPolygonOnPlaneVertexIteratorIterator;
        }
        addPolygonOnPlaneVertex( *secondListListPolygonOnPlaneVertexIteratorIterator, false );
        break;
      }
    }
  }
  if( pPlaneIterIterator == pPlaneEndIterIterator )// not found this pair
  {
    addPolygonOnPlaneVertex( secondVertexIterator, bGenerateSamplePoints );
  }
}

void WireFrame::addLineToPolygonOnPlane( ListPolygonOnPlaneVertexIterator & firstVertexIterator,
    ListPolygonOnPlaneVertexIterator & secondVertexIterator, const TooN::SE3<> & se3W2C,
    const CVD::Image<CVD::byte> & image, bool bGenerateSamplePoints )
{
  ListListPolygonOnPlaneIteratorIterator pPlaneIterIterator =
      firstVertexIterator->listListPolygonOnPlaneIterator.begin();
  ListListPolygonOnPlaneIteratorIterator pPlaneEndIterIterator =
      firstVertexIterator->listListPolygonOnPlaneIterator.end();
  for( ; pPlaneIterIterator != pPlaneEndIterIterator; ++pPlaneIterIterator ) //for each plane that firstVertexIterator belongs to.
  {
    if( ( *( *pPlaneIterIterator ) ).controlListListPolygonOnPlaneVertexIterator.size() < 2 ) continue;
    ListListPolygonOnPlaneVertexIteratorIterator nextVertexIterIterator = find(
        ( *( *pPlaneIterIterator ) ).controlListListPolygonOnPlaneVertexIterator.begin(),
        ( *( *pPlaneIterIterator ) ).controlListListPolygonOnPlaneVertexIterator.end(), firstVertexIterator );
    ListListPolygonOnPlaneVertexIteratorIterator previousVertexIterIter = nextVertexIterIterator;
    if( nextVertexIterIterator != ( *( *pPlaneIterIterator ) ).controlListListPolygonOnPlaneVertexIterator.end() )
    {
      ++nextVertexIterIterator;
      if( *nextVertexIterIterator == secondVertexIterator ) // found in this plane.
      {
        ListListPolygonOnPlaneVertexIteratorIterator firstListListPolygonOnPlaneVertexIteratorIterator = find(
            ( *( *pPlaneIterIterator ) ).listListPolygonOnPlaneVertexIterator.begin(),
            ( *( *pPlaneIterIterator ) ).listListPolygonOnPlaneVertexIterator.end(), firstVertexIterator );
        ++firstListListPolygonOnPlaneVertexIteratorIterator; // not include the firstListListPolygonOnPlaneVertexIteratorIterator
        ListListPolygonOnPlaneVertexIteratorIterator secondListListPolygonOnPlaneVertexIteratorIterator = find(
            firstListListPolygonOnPlaneVertexIteratorIterator,
            ( *( *pPlaneIterIterator ) ).listListPolygonOnPlaneVertexIterator.end(), secondVertexIterator );
        while( firstListListPolygonOnPlaneVertexIteratorIterator != secondListListPolygonOnPlaneVertexIteratorIterator )
        {
          //          std::cout << (*(*firstListListPolygonOnPlaneVertexIteratorIterator)).v3WorldCoordinate << "\t" << (*(*secondListListPolygonOnPlaneVertexIteratorIterator)).v3WorldCoordinate << std::endl;
          addPolygonOnPlaneVertex( *firstListListPolygonOnPlaneVertexIteratorIterator, false );
          ++firstListListPolygonOnPlaneVertexIteratorIterator;
        }
        addPolygonOnPlaneVertex( *secondListListPolygonOnPlaneVertexIteratorIterator, false );
        break;
      }
    }
    if( previousVertexIterIter != ( *( *pPlaneIterIterator ) ).controlListListPolygonOnPlaneVertexIterator.begin() )
    {
      --previousVertexIterIter;
      if( *previousVertexIterIter == secondVertexIterator ) // found in this plane.
      {
        ListListPolygonOnPlaneVertexIteratorIterator secondListListPolygonOnPlaneVertexIteratorIterator = find(
            ( *( *pPlaneIterIterator ) ).listListPolygonOnPlaneVertexIterator.begin(),
            ( *( *pPlaneIterIterator ) ).listListPolygonOnPlaneVertexIterator.end(), secondVertexIterator );
        ListListPolygonOnPlaneVertexIteratorIterator firstListListPolygonOnPlaneVertexIteratorIterator = find(
            secondListListPolygonOnPlaneVertexIteratorIterator,
            ( *( *pPlaneIterIterator ) ).listListPolygonOnPlaneVertexIterator.end(), firstVertexIterator );
        --firstListListPolygonOnPlaneVertexIteratorIterator; // not include the firstListListPolygonOnPlaneVertexIteratorIterator
        while( firstListListPolygonOnPlaneVertexIteratorIterator != secondListListPolygonOnPlaneVertexIteratorIterator )
        {
          addPolygonOnPlaneVertex( *firstListListPolygonOnPlaneVertexIteratorIterator, false );
          --firstListListPolygonOnPlaneVertexIteratorIterator;
        }
        addPolygonOnPlaneVertex( *secondListListPolygonOnPlaneVertexIteratorIterator, false );
        break;
      }
    }
  }
  if( pPlaneIterIterator == pPlaneEndIterIterator )// not found this pair
  {
    addPolygonOnPlaneVertex( secondVertexIterator, se3W2C, image, bGenerateSamplePoints );
  }
}

void WireFrame::addPolygonOnPlane( const std::vector<TooN::Vector<3> > & vv3poses, const TooN::Vector<4> & v4Plane, PolygonOnPlaneType::Type type, bool bGenerateSamplePoints )
{
  createPolygonOnPlane( type );
  setPlane( v4Plane );
  std::vector<TooN::Vector<3> >::const_iterator v3Iter = vv3poses.begin();
  std::vector<TooN::Vector<3> >::const_iterator v3EndIter = vv3poses.end();
  bool bFoundFirstVertex = false;
  bool bFound = false;
  ListPolygonOnPlaneVertexIterator firstVertexIterator;
  for( ; v3Iter != v3EndIter; ++v3Iter )
  {
    ListPolygonOnPlaneIterator pPlaneIterator = mActiveObject->listPolygonOnPlane.begin();
    ListPolygonOnPlaneIterator pPlaneEndIterator = mActiveObject->listPolygonOnPlane.end();
    bFound = false;
    // check new points againt all control points
    for( ; pPlaneIterator != pPlaneEndIterator; ++pPlaneIterator ) // for every plane.
    {
      ListListPolygonOnPlaneVertexIteratorIterator listListVertexIteratorIterator =
          ( *pPlaneIterator ).controlListListPolygonOnPlaneVertexIterator.begin();
      ListListPolygonOnPlaneVertexIteratorIterator listListVertexEndIteratorIterator =
          ( *pPlaneIterator ).controlListListPolygonOnPlaneVertexIterator.end();
      for( ; listListVertexIteratorIterator != listListVertexEndIteratorIterator; ++listListVertexIteratorIterator )
      {
        if( const_cast<TooN::Vector<3> &> ( *v3Iter ) == ( *( *listListVertexIteratorIterator ) ).v3WorldCoordinate )
        {
          bFound = true;
          if( bFoundFirstVertex ) // bFoundFirstVertex == true
          {
            addLineToPolygonOnPlane( firstVertexIterator, *listListVertexIteratorIterator, bGenerateSamplePoints );
          }
          else // bFoundFirstVertex == false
          {
            bFoundFirstVertex = true;
            addPolygonOnPlaneVertex( *listListVertexIteratorIterator, bGenerateSamplePoints );
          }
          firstVertexIterator = *listListVertexIteratorIterator;
          break;
        }
      }
      if( bFound ) break;
    }
    if( !bFound )
    {
      bFoundFirstVertex = false;
      addPolygonOnPlaneVertex( *v3Iter, bGenerateSamplePoints );
    }
  }
  if( type == PolygonOnPlaneType::RECTANGLE && bGenerateSamplePoints == false )
  {
    TooN::Vector<3> v3DistanceFstEdge = vv3poses[1] - vv3poses[0];
    TooN::Vector<3> v3DistanceSndEdge = vv3poses[2] - vv3poses[1];

    TooN::Vector<3> v3Start, v3Stop;
    if( TooN::norm(v3DistanceFstEdge) < 0.1f*TooN::norm(v3DistanceSndEdge) )
    {
      v3Start = 0.5f*( vv3poses[1] + vv3poses[0] );
      v3Stop = 0.5f*( vv3poses[3] + vv3poses[2] );
    }
    else
    {
      v3Start = 0.5f*( vv3poses[2] + vv3poses[1] );
      v3Stop = 0.5f*( vv3poses[4] + vv3poses[3] );
    }
    addBoundaryVertex( v3Start );
    addBoundaryVertex( v3Stop );
  }
}

void WireFrame::addPolygonOnPlane( const std::vector<TooN::Vector<3> > & vv3poses, const TooN::Vector<4> & v4Plane,
    const TooN::SE3<> & se3W2C, const CVD::Image<CVD::byte> & image, PolygonOnPlaneType::Type type, bool bGenerateSamplePoints )
{
  createPolygonOnPlane( type );
  setPlane( v4Plane );
  std::vector<TooN::Vector<3> >::const_iterator v3Iter = vv3poses.begin();
  std::vector<TooN::Vector<3> >::const_iterator v3EndIter = vv3poses.end();
  bool bFoundFirstVertex = false;
  bool bFound = false;
  ListPolygonOnPlaneVertexIterator firstVertexIterator;
  for( ; v3Iter != v3EndIter; ++v3Iter )
  {
    ListPolygonOnPlaneIterator pPlaneIterator = mActiveObject->listPolygonOnPlane.begin();
    ListPolygonOnPlaneIterator pPlaneEndIterator = mActiveObject->listPolygonOnPlane.end();
    bFound = false;
    // check new points againt all control points
    for( ; pPlaneIterator != pPlaneEndIterator; ++pPlaneIterator ) // for every plane.
    {
      ListListPolygonOnPlaneVertexIteratorIterator listListVertexIteratorIterator =
          ( *pPlaneIterator ).controlListListPolygonOnPlaneVertexIterator.begin();
      ListListPolygonOnPlaneVertexIteratorIterator listListVertexEndIteratorIterator =
          ( *pPlaneIterator ).controlListListPolygonOnPlaneVertexIterator.end();
      for( ; listListVertexIteratorIterator != listListVertexEndIteratorIterator; ++listListVertexIteratorIterator )
      {
        if( const_cast<TooN::Vector<3> &> ( *v3Iter ) == ( *( *listListVertexIteratorIterator ) ).v3WorldCoordinate )
        {
          bFound = true;
          if( bFoundFirstVertex ) // bFoundFirstVertex == true
          {
            addLineToPolygonOnPlane( firstVertexIterator, *listListVertexIteratorIterator, se3W2C, image, bGenerateSamplePoints );
          }
          else // bFoundFirstVertex == false
          {
            bFoundFirstVertex = true;
            addPolygonOnPlaneVertex( *listListVertexIteratorIterator, se3W2C, image, bGenerateSamplePoints );
          }
          firstVertexIterator = *listListVertexIteratorIterator;
          break;
        }
      }
      if( bFound ) break;
    }
    if( !bFound )
    {
      bFoundFirstVertex = false;
      addPolygonOnPlaneVertex( *v3Iter, bGenerateSamplePoints );
      //      std::cout <<"Add new points : " << *v3Iter << std::endl;
    }
  }
}

void WireFrame::setPlane( const TooN::Vector<4> & v4Plane )
{
  assert( mActiveObject != mListObject.end() && mActivePolygonOnPlane != mActiveObject->listPolygonOnPlane.end() );
  mActivePolygonOnPlane->v4Plane = v4Plane;
}

void WireFrame::closePolygonOnPlane( bool bGenerateSamplePoints )
{
  assert( mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.size() >= 3 );
  ListListPolygonOnPlaneVertexIteratorIterator tailIterator =
      mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.end();
  --tailIterator;
  addPolygonOnPlaneVertex( ( *tailIterator ), bGenerateSamplePoints );
  /*  TooN::Vector<3> v3direction = mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.back()->v3WorldCoordinate - mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.front()->v3WorldCoordinate;
   REAL_TYPE no_segment = 8.0*sqrt(v3direction*v3direction);
   TooN::Vector<3> v3Step = v3direction/no_segment;
   normalize(v3direction);
   TooN::Vector<3> v3AutoPose = mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.front()->v3WorldCoordinate + v3Step;
   while( no_segment > 1.0f )
   {
   PolygonOnPlaneVertex pgvAuto;
   pgvAuto.v3WorldCoordinate = v3AutoPose;
   pgvAuto.v3WorldCoordinate_plus_direction = v3AutoPose - 0.01f*v3direction;
   pgvAuto.v3Normal = mActivePolygonOnPlane->v4Plane.slice(0,3);
   pgvAuto.bAuto = true;
   pgvAuto.listListPolygonOnPlaneIterator.push_front(mActivePolygonOnPlane);
   addPolygonOnPlaneVertex( pgvAuto );
   no_segment -= 1.f;
   v3AutoPose += v3Step;
   }
   ListListPolygonOnPlaneVertexIteratorIterator tailIterator = mActivePolygonOnPlane->listListPolygonOnPlaneVertexIterator.end();
   --tailIterator;
   addPolygonOnPlaneVertex((*tailIterator)); */
}

/*
 void WireFrame::removeVertex( VertexListIterator itr )
 {

 }

 void WireFrame::selectObject( int noObject )
 {
 if( noObject > miNoObjects )
 return;
 ListObjectIterator itr;
 for( itr = mListObject.begin(); itr != mListObject.end(); ++itr )
 {
 noObject--;
 if( noObject < 0 )
 break;
 }
 mActiveObject = itr;
 }

 void WireFrame::moveActiveVertex( TooN::Vector<3> toPose )
 {
 (*mActiveVertex).v3WorldCoordinate = toPose;
 }
 */

void WireFrame::resetTrackingStatus()
{
  for( ListObjectIterator listObjectIterator = mListObject.begin(); listObjectIterator != mListObject.end(); ++listObjectIterator )
  {
    Object & object = ( *listObjectIterator );
    for( ListFreeVertexIterator itr = object.listFreeVertex.begin(); itr != object.listFreeVertex.end(); ++itr )
    {
      ( *itr ).bTracked = false;
    }
    for( ListPolygonIterator itr = object.listPolygon.begin(); itr != object.listPolygon.end(); ++itr )
    {
      if( (*itr).bTracked )
      {
        ListListPolygonVertexIteratorIterator itrVertex = (*itr).listListPolygonVertexIterator.begin();
        ListListPolygonVertexIteratorIterator itrEndVertex = (*itr).listListPolygonVertexIterator.end();
        for(; itrVertex != itrEndVertex; ++itrVertex )
          (*(*itrVertex)).bTracked = false;
        ( *itr ).bTracked = false;
      }
    }
    for( ListPolygonOnPlaneIterator itr = object.listPolygonOnPlane.begin(); itr != object.listPolygonOnPlane.end(); ++itr )
    {
      if( (*itr).bTracked )
      {
        ListListPolygonOnPlaneVertexIteratorIterator itrVertex = (*itr).listListPolygonOnPlaneVertexIterator.begin();
        ListListPolygonOnPlaneVertexIteratorIterator itrEndVertex = (*itr).listListPolygonOnPlaneVertexIterator.end();
        for(; itrVertex != itrEndVertex; ++itrVertex )
          (*(*itrVertex)).bTracked = false;
        itrVertex = (*itr).listListBoundaryVertexIterator.begin();
        itrEndVertex = (*itr).listListBoundaryVertexIterator.end();
        for(; itrVertex != itrEndVertex; ++itrVertex )
          (*(*itrVertex)).bTracked = false;
        ( *itr ).bTracked = false;
      }
    }
  }
}

void WireFrame::printReport()
{
  std::cout << "Number of Object  : " << miNoObject << std::endl;
  for( ListObjectIterator listObjectIterator = mListObject.begin(); listObjectIterator != mListObject.end(); ++listObjectIterator )
  {
    Object & object = ( *listObjectIterator );
    std::cout << "Number of Free Vertex  : " << object.listFreeVertex.size() << std::endl;
    for( ListFreeVertexIterator itr = object.listFreeVertex.begin(); itr != object.listFreeVertex.end(); ++itr )
    {
      std::cout << ( *itr ).v3WorldCoordinate << std::endl;
    }
    std::cout << "Number of Polygon Vertex  : " << object.listPolygonVertex.size() << std::endl;
    for( ListPolygonVertexIterator itr = object.listPolygonVertex.begin(); itr != object.listPolygonVertex.end(); ++itr )
    {
      std::cout << ( *itr ).v3WorldCoordinate << std::endl;
    }
    std::cout << "Number of Polygon On Plane Vertex  : " << object.listPolygonOnPlaneVertex.size() << std::endl;
    for( ListPolygonOnPlaneVertexIterator itr = object.listPolygonOnPlaneVertex.begin(); itr
        != object.listPolygonOnPlaneVertex.end(); ++itr )
    {
      std::cout << ( *itr ).v3WorldCoordinate << std::endl;
    }
    std::cout << "Number of Polygon : " << object.listPolygon.size() << std::endl;
    int countPolygon = 1;
    for( ListPolygonIterator itr = object.listPolygon.begin(); itr != object.listPolygon.end(); ++itr, ++countPolygon )
    {
      std::cout << "No. Polygon : " << countPolygon << " (" << ( *itr ).listListPolygonVertexIterator.size() << ")"
          << std::endl;
      for( ListListPolygonVertexIteratorIterator itrVert = ( *itr ).listListPolygonVertexIterator.begin(); itrVert
          != ( *itr ).listListPolygonVertexIterator.end(); ++itrVert )
      {
        std::cout << ( *( *itrVert ) ).v3WorldCoordinate << " " << ( *( *itrVert ) ).v3WorldCoordinate_plus_direction
            << std::endl;
      }
    }
    std::cout << "Number of Polygon on Plane : " << object.listPolygonOnPlane.size() << std::endl;
    countPolygon = 1;
    for( ListPolygonOnPlaneIterator itr = object.listPolygonOnPlane.begin(); itr != object.listPolygonOnPlane.end(); ++itr, ++countPolygon )
    {
      std::cout << "No. Polygon on Plane : " << countPolygon << " ("
          << ( *itr ).listListPolygonOnPlaneVertexIterator.size() << ")" << std::endl;
      std::cout << "Control Points " << std::endl;
      for( ListListPolygonOnPlaneVertexIteratorIterator itrVert =
          ( *itr ).controlListListPolygonOnPlaneVertexIterator.begin(); itrVert
          != ( *itr ).controlListListPolygonOnPlaneVertexIterator.end(); ++itrVert )
      {
        std::cout << ( *( *itrVert ) ).v3WorldCoordinate << " " << ( *( *itrVert ) ).v3WorldCoordinate_plus_direction
            << std::endl;
      }
      std::cout << "Sample Points " << std::endl;
      for( ListListPolygonOnPlaneVertexIteratorIterator itrVert = ( *itr ).listListPolygonOnPlaneVertexIterator.begin(); itrVert
          != ( *itr ).listListPolygonOnPlaneVertexIterator.end(); ++itrVert )
      {
        std::cout << ( *( *itrVert ) ).v3WorldCoordinate << " " << ( *( *itrVert ) ).v3WorldCoordinate_plus_direction
            << std::endl;
      }
    }
  }
}

