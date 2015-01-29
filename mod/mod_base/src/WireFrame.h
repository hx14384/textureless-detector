/*
  Computer Vision Group, Department of Computer Science. University of      
  Bristol. This code can not be used or copied without express permission.   
*/

#ifndef WIREFRAME_H
#define WIREFRAME_H

#include <list>
#include <vector>
#include <TooN/TooN.h>
#include <TooN/se3.h>
#include <cvd/image.h>
#include <cvd/byte.h>

// TODO: ? add bourder feature that will be visible when its normal direction perpendicular to viewing direction.
struct Vertex;
struct FreeVertex;
typedef std::list<FreeVertex> ListFreeVertex; // contains actual Vertex object.
typedef std::list<FreeVertex>::iterator ListFreeVertexIterator;
typedef std::list<FreeVertex>::reverse_iterator ReverseListFreeVertexIterator;
typedef std::list<ListFreeVertexIterator> ListListFreeVertexIterator; // contains a pointer to Vertex object.
typedef std::list<ListFreeVertexIterator>::iterator ListListFreeVertexIteratorIterator;

struct PolygonVertex;
typedef std::list<PolygonVertex> ListPolygonVertex;
typedef std::list<PolygonVertex>::iterator ListPolygonVertexIterator;
typedef std::list<PolygonVertex>::reverse_iterator ReverseListPolygonVertexIterator;
typedef std::list<ListPolygonVertexIterator> ListListPolygonVertexIterator;
typedef std::list<ListPolygonVertexIterator>::iterator ListListPolygonVertexIteratorIterator;
typedef std::list<ListPolygonVertexIterator>::reverse_iterator ReverseListListPolygonVertexIteratorIterator;

struct PolygonOnPlaneVertex;
typedef std::list<PolygonOnPlaneVertex> ListPolygonOnPlaneVertex;
typedef std::list<PolygonOnPlaneVertex>::iterator ListPolygonOnPlaneVertexIterator;
typedef std::list<PolygonOnPlaneVertex>::reverse_iterator ReverseListPolygonOnPlaneVertexIterator;
typedef std::list<ListPolygonOnPlaneVertexIterator> ListListPolygonOnPlaneVertexIterator;
typedef std::list<ListPolygonOnPlaneVertexIterator>::iterator ListListPolygonOnPlaneVertexIteratorIterator;
typedef std::list<ListPolygonOnPlaneVertexIterator>::reverse_iterator ReverseListListPolygonOnPlaneVertexIteratorIterator;

struct Polygon;
typedef std::list<Polygon> ListPolygon;
typedef std::list<Polygon>::iterator ListPolygonIterator;
typedef std::list<Polygon>::reverse_iterator ReverseListPolygonIterator;
typedef std::list<ListPolygonIterator> ListListPolygonIterator;
typedef std::list<ListPolygonIterator>::iterator ListListPolygonIteratorIterator;

struct PolygonOnPlane;
typedef std::list<PolygonOnPlane> ListPolygonOnPlane;
typedef std::list<PolygonOnPlane>::iterator ListPolygonOnPlaneIterator;
typedef std::list<PolygonOnPlane>::reverse_iterator ReverseListPolygonOnPlaneIterator;
typedef std::list<ListPolygonOnPlaneIterator> ListListPolygonOnPlaneIterator;
typedef std::list<ListPolygonOnPlaneIterator>::iterator ListListPolygonOnPlaneIteratorIterator;

struct Object;
typedef std::list<Object> ListObject;
typedef std::list<Object>::iterator ListObjectIterator;
typedef std::list<Object>::reverse_iterator ReverseListObjectIterator;

struct Vertex
{
    bool bTracked;
    int iID;
    TooN::Vector<3> v3CameraCoordinate; // in camera coordinate.
    TooN::Vector<3> v3WorldCoordinate; // in world coordinate.
    TooN::Vector<2> v2Image_position; // in 2D image.
    Vertex() :
      bTracked( false )
    {
    }
    ;
    Vertex & operator =( const Vertex & rhs )
    {
      if( this != &rhs )
      {
        bTracked = rhs.bTracked;
        iID = rhs.iID;
        v3CameraCoordinate = rhs.v3CameraCoordinate;
        v3WorldCoordinate = rhs.v3WorldCoordinate;
        v2Image_position = rhs.v2Image_position;
      }
      return *this;
    }
    ;
};

struct FreeVertex: public Vertex
{
    TooN::Vector<3> v3Normal; // a v3Normal direction.
    FreeVertex() :
      Vertex(), v3Normal( TooN::makeVector( 0.f, 0.f, 1.f ) )
    {
    }
    ;
    FreeVertex & operator =( const FreeVertex & rhs )
    {
      if( this != &rhs )
      {
        Vertex::operator=( rhs );
        v3Normal = rhs.v3Normal;
      }
      return *this;
    }
    ;
};

#define EDGE_FEATURE      1
#define EDGE_POLE         2
// x x x x x x pole isSet //
struct EdgeFeature
{
    unsigned char ucBits;
    EdgeFeature() :
      ucBits( 0 )
    {
    }
    ;
    EdgeFeature & operator =( const EdgeFeature & rhs )
    {
      if( this != &rhs )
      {
        ucBits = rhs.ucBits;
      }
      return *this;
    }
    ;
};

struct PolygonVertex: public Vertex
{
    bool bAuto;
    TooN::Vector<3> v3Normal; // a v3Normal direction.
    TooN::Vector<3> v3WorldCoordinate_plus_direction; // to find a 2D normal direction if this vertex belongs to any polygon, but no need for free vertices.
    ListListPolygonIterator listListPolygonIterator; // this vertex is in these listListPolygonIterator.
    EdgeFeature edgeFeature;
    PolygonVertex() :
      Vertex(), bAuto( false ), v3Normal( TooN::makeVector( 0.f, 0.f, 1.f ) )
    {
    }
    ;
    PolygonVertex & operator =( const PolygonVertex & rhs )
    {
      if( this != &rhs )
      {
        Vertex::operator=( rhs );
        bAuto = rhs.bAuto;
        v3WorldCoordinate_plus_direction = rhs.v3WorldCoordinate_plus_direction;
        v3Normal = rhs.v3Normal;
        listListPolygonIterator = rhs.listListPolygonIterator;
        edgeFeature = rhs.edgeFeature;
      }
      return *this;
    }
    ;
};

struct PolygonOnPlaneVertex: public Vertex
{
    bool bAuto;
    TooN::Vector<3> v3Normal; // may not be the same as plane's notmal if it's shared between two planes. // may not need this ;-)
    TooN::Vector<3> v3WorldCoordinate_plus_direction; // to find a 2D normal direction if this vertex belongs to any polygon, but no need for free vertices.
    ListListPolygonOnPlaneIterator listListPolygonOnPlaneIterator;
    EdgeFeature edgeFeature;
    PolygonOnPlaneVertex() :
      Vertex(), bAuto( false ), v3Normal( TooN::makeVector( 0.f, 0.f, 1.f ) )
    {
    }
    ;
    PolygonOnPlaneVertex & operator =( const PolygonOnPlaneVertex & rhs )
    {
      if( this != &rhs )
      {
        Vertex::operator=( rhs );
        bAuto = rhs.bAuto;
        v3WorldCoordinate_plus_direction = rhs.v3WorldCoordinate_plus_direction;
        v3Normal = rhs.v3Normal;
        listListPolygonOnPlaneIterator = rhs.listListPolygonOnPlaneIterator;
        edgeFeature = rhs.edgeFeature;
      }
      return *this;
    }
    ;
};

struct Polygon
{
    bool bTracked;
    int iID;
    int iNoVertex;
    ListListPolygonVertexIterator listListPolygonVertexIterator;
    int iNoControlVertex;
    ListListPolygonVertexIterator controlListListPolygonVertexIterator;
    Polygon() :
      bTracked( false ), iNoVertex(0), iNoControlVertex(0)
    {
    }
    ;
    Polygon & operator =( const Polygon & rhs )
    {
      if( this != &rhs )
      {
        bTracked = rhs.bTracked;
        iID = rhs.iID;
        iNoVertex = rhs.iNoVertex;
        listListPolygonVertexIterator = rhs.listListPolygonVertexIterator;
        iNoControlVertex = rhs.iNoControlVertex;
        controlListListPolygonVertexIterator = rhs.controlListListPolygonVertexIterator;
      }
      return *this;
    }
    ;
};

struct PolygonOnPlaneType
{
  enum Type
  {
    NONE = 0, POLYGON, RECTANGLE, CIRCLE, SQUARE
  }
  ;
}
;

struct PolygonOnPlane
{
    bool bTracked;
    int iID;
    PolygonOnPlaneType::Type type;
    TooN::Vector<4> v4Plane;
    TooN::Vector<3> v3Centroid;
    int iNoVertex;
    ListListPolygonOnPlaneVertexIterator listListPolygonOnPlaneVertexIterator;
    int iNoControlVertex;
    ListListPolygonOnPlaneVertexIterator controlListListPolygonOnPlaneVertexIterator;
    int iNoBoundaryVertex;
    ListListPolygonOnPlaneVertexIterator listListBoundaryVertexIterator;
    PolygonOnPlane() :
      bTracked( false ), v4Plane( TooN::makeVector( 0.f, 1.f, 0.f, 0.f ) ), v3Centroid( TooN::Zeros ), iNoVertex(0), iNoControlVertex(0), iNoBoundaryVertex(0)
    {
    }
    ;
    PolygonOnPlane & operator =( const PolygonOnPlane & rhs )
    {
      if( this != &rhs )
      {
        bTracked = rhs.bTracked;
        iID = rhs.iID;
        type = rhs.type;
        v4Plane = rhs.v4Plane;
        v3Centroid = rhs.v3Centroid;
        iNoVertex = rhs.iNoVertex;
        listListPolygonOnPlaneVertexIterator = rhs.listListPolygonOnPlaneVertexIterator;
        iNoControlVertex = rhs.iNoControlVertex;
        controlListListPolygonOnPlaneVertexIterator = rhs.controlListListPolygonOnPlaneVertexIterator;
        iNoBoundaryVertex = rhs.iNoBoundaryVertex;
        listListBoundaryVertexIterator = rhs.listListBoundaryVertexIterator;
      }
      return *this;
    }
    ;
};

struct Object
{
    bool bTracked;
    int iID;
    ListFreeVertex listFreeVertex;
    ListPolygonVertex listPolygonVertex;
    ListPolygonOnPlaneVertex listPolygonOnPlaneVertex;
    ListPolygonOnPlaneVertex listBoundaryVertex;
    ListPolygon listPolygon;
    ListPolygonOnPlane listPolygonOnPlane;
    TooN::Matrix<4, 4> m4Transform;
    Object( int iObjectID ) :
      bTracked( true ), iID( iObjectID ), m4Transform( TooN::Identity )
    {
    }
    ;
    Object & operator=( const Object & rhs )
    {
      if( this != &rhs )
      {
        bTracked = rhs.bTracked;
        iID = rhs.iID;
        listFreeVertex = rhs.listFreeVertex;
        listPolygonVertex = rhs.listPolygonVertex;
        listPolygonOnPlaneVertex = rhs.listPolygonOnPlaneVertex;
        listBoundaryVertex = rhs.listBoundaryVertex;
        listPolygon = rhs.listPolygon;
        listPolygonOnPlane = rhs.listPolygonOnPlane;
        m4Transform = rhs.m4Transform;
      }
      return *this;
    }
    ;
};
/**
  @file WireFrame.h
  @brief WireFrame will manage all information about 3D objects used in the program. 
*/
class WireFrame
{
  public:
    WireFrame();
    ~WireFrame();
    // higher level features
    void empty();
    void createObject();
    void createPolygon();
    void createPolygonOnPlane( PolygonOnPlaneType::Type type = PolygonOnPlaneType::NONE );
    void setPlane( const TooN::Vector<4> & v4Plane );

    void addFreeVertex( const TooN::Vector<3> & v3pose, const TooN::Vector<3> & v3Normal );
    void addFreeVertex( const FreeVertex & freeVertex );
    //////////////////////////////////////////////////////////////
    void addPolygonVertex( const TooN::Vector<3> & v3pose, const TooN::Vector<3> & v3Normal, bool bGenerateSamplePoints = true );
    void addPolygonVertex( const PolygonVertex & polygonVertex );
    //////////////////////////////////////////////////////////////
    void addPolygonOnPlaneVertex( const TooN::Vector<3> & v3pose, bool bGenerateSamplePoints = true );
    void addPolygonOnPlaneVertex( const TooN::Vector<3> & v3pose, const TooN::SE3<> & se3W2C, 
            const CVD::Image< CVD::byte> & image, bool bGenerateSamplePoints = true );
    void addPolygonOnPlaneVertex( const PolygonOnPlaneVertex & polygonOnPlaneVertex );
    void addPolygonOnPlaneVertex( PolygonOnPlaneVertex & polygonOnPlaneVertex, ListObjectIterator & listObjectIterator, ListPolygonOnPlaneIterator & listPolygonOnPlaneIterator );
    void addPolygonOnPlaneVertex( ListPolygonOnPlaneVertexIterator & listPolygonOnPlaneVertexIterator, bool bGenerateSamplePoints = true );
    void addPolygonOnPlaneVertex( ListPolygonOnPlaneVertexIterator & listPolygonOnPlaneVertexIterator, 
            const TooN::SE3<> & se3W2C, const CVD::Image<CVD::byte> & image, bool bGenerateSamplePoints = true );
    void addBoundaryVertex( const TooN::Vector<3> & v3pose );
    void addBoundaryVertex( const PolygonOnPlaneVertex & polygonOnPlaneVertex );
    //////////////////////////////////////////////////////////////
    void addLineToPolygonOnPlane( ListPolygonOnPlaneVertexIterator & firstVertexIterator, 
            ListPolygonOnPlaneVertexIterator & secondVertexIterator, bool bGenerateSamplePoints = true );
    void addLineToPolygonOnPlane( ListPolygonOnPlaneVertexIterator & firstVertexIterator, 
            ListPolygonOnPlaneVertexIterator & secondVertexIterator, const TooN::SE3<> & se3W2C,
            const CVD::Image<CVD::byte> & image, bool bGenerateSamplePoints = true );
    //////////////////////////////////////////////////////////////
    void addPolygonOnPlane( const std::vector<TooN::Vector<3> > & vv3poses, const TooN::Vector<4> & v4Plane,
            PolygonOnPlaneType::Type type = PolygonOnPlaneType::NONE, bool bGenerateSamplePoints = true );
    void addPolygonOnPlane( const std::vector<TooN::Vector<3> > & vv3poses, const TooN::Vector<4> & v4Plane, 
            const TooN::SE3<> & se3W2C, const CVD::Image<CVD::byte> & image, PolygonOnPlaneType::Type type = PolygonOnPlaneType::NONE, bool bGenerateSamplePoints = true );
    //////////////////////////////////////////////////////////////
    void closePolygon( bool bGenerateSamplePoints = true );
    void closePolygonOnPlane( bool bGenerateSamplePoints = true );
    //////////////////////////////////////////////////////////////
    void resetTrackingStatus();
    /*    void selectVertex( TooN::Vector<2> & v2pose );
     void selectPlane( TooN::Vector<2> & v2pose );
     void selectObject( TooN::Vector<2> & v2pose );
     void moveActiveVertex( TooN::Vector<3> v3pose );
     */
    void printReport();
    //  private:
    int miNoObject;
    int miNoFreeVertex;
    int miNoPolygonVertex;
    int miNoPolygonOnPlaneVertex;
    int miNoBoundaryVertex;
    int miNoPolygon;
    int miNoPolygonOnPlane;
    ListObject mListObject;
    ListFreeVertexIterator mActiveFreeVertex;
    ListPolygonVertexIterator mActivePolygonVertex;
    ListPolygonOnPlaneVertexIterator mActivePolygonOnPlaneVertex;
    ListPolygonOnPlaneVertexIterator mActiveBoundaryVertex;
    ListPolygonIterator mActivePolygon;
    ListPolygonOnPlaneIterator mActivePolygonOnPlane;
    ListObjectIterator mActiveObject;
};
#endif
