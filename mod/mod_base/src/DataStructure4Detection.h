/*
 * =====================================================================================
 *
 *       Filename:  DataStructure4Detection.h
 *
 *    Description:  
 *
 *        Version:  2.0
 *        Created:  27/07/10 01:43:46
 *       Compiler:  gcc
 *
 *        Company:  Computer Science Department, University of Bristol
 *
 * =====================================================================================
 */

#ifndef DATASTRUCTURE4DETECTION_H
#define DATASTRUCTURE4DETECTION_H

#include "DefineTypes.h"
#include <TooN/se3.h>
#include <TooN/TooN.h>
#include <cvd/image.h>
#include <cvd/image_ref.h>
#include <sensor_msgs/Image.h>
#include <vector>

//#define NO_TRAIN_QUESTIONS 3000  // this should be now irrelevant
//#define NO_QUESTIONS 500
#define QUESTIONS_SIZE 4 // this equals the length of the chain - 1. It signifies the number of pairs in the chain
#define NO_MAX_BINS 64 // ???
//#define BIN_RANGE 2 // should be irrelevant now
//#define STEP 0.05
//#define INV_STEP 20

// Contain a long line segment;
struct Line
{
  TooN::Vector<2> v2Start;      // starting point of the line 
  TooN::Vector<2> v2End;        // ending point of the line
  TooN::Vector<2> v2Direction;  // a unit direction pointing from v2Start to v2End
  REAL_TYPE rSlope;		// the slope of the line m in (y = mx + c)
  int no_pixel;                 // number of pixels in the line both in width and height
  REAL_TYPE distanceTransformAverage;
};

// Edgelet is a short line segment;
// Keep the centre and slope of the line;
struct Edgelet
{
  TooN::Vector<2> v2Pose; // do ir(v2pose) before saving to file
//  REAL_TYPE rAngle;
  REAL_TYPE rSlope;
  int regionNo;
  int previousObjNo;
};

struct EdgeLink // ???
{
  TooN::Vector<2> v2Pose; // snd edge
  REAL_TYPE rDir; // may not need
  REAL_TYPE rAngle; // may not need
  REAL_TYPE rDistance;
  TooN::Vector<2> v2Check; // used for checking a side of point
};

// Put EdgeLinks generated from Edgelet into irBins. // ???
struct BinsOfEdgeLinks 
{
  BinsOfEdgeLinks() : mbFound(false) { memset(mBitBins, 1, NO_MAX_BINS*NO_MAX_BINS*sizeof(bool)); };
  TooN::Vector<2> v2Pose; // fst edge
  bool mBitBins[NO_MAX_BINS][NO_MAX_BINS];
  std::vector<EdgeLink> mEdgeLinkBins[NO_MAX_BINS][NO_MAX_BINS];
  bool mbFound;
};

/*struct UCCoordinate
{
  unsigned char r;
  unsigned char c;
};

struct Question
{
  REAL_TYPE rDir; // may not need
  REAL_TYPE rAngle; // may not need
  REAL_TYPE rDistance;
  REAL_TYPE rRelativeDistance;
  bool bSide; // side of 3rd, 4th and 5th points relatived to line of 1st-2nd
  unsigned char ucCentre_r; // bins that might be belonged to this Question.
  unsigned char ucCentre_c;
  std::vector<UCCoordinate> bins;
  Question & operator =( const Question & rhs )
  {
    if( this != &rhs )
    {
      rDir = rhs.rDir;
      rAngle = rhs.rAngle;
      rDistance = rhs.rDistance;
      rRelativeDistance = rhs.rRelativeDistance;
      bSide = rhs.bSide;
      ucCentre_r = rhs.ucCentre_r;
      ucCentre_c = rhs.ucCentre_c;
      bins = rhs.bins;
    }
    return *this;
  };
};

// Table will be in std::vector< EdgeLinkQuestion > form.
struct EdgeLinkQuestion
{
  EdgeLinkQuestion() : vQuestions(4), vPoints(5) {}
  ~EdgeLinkQuestion() {}
  
  std::vector<Question> vQuestions;
  std::vector<TooN::Vector<2> > vPoints;
  EdgeLinkQuestion & operator =( const EdgeLinkQuestion & rhs )
  {
    if( this != &rhs )
    {
      vQuestions = rhs.vQuestions;
      vPoints = rhs.vPoints;
    }
    return *this;
  };
};

struct ObjectDescriptor
{
  ObjectDescriptor() : vse3W2C(2), vv3VisiblePoints(10) {};
  std::vector<EdgeLinkQuestion> vEdgeLinkQuestions;
  TooN::SE3<> se3W2COriginal;
  std::vector<TooN::SE3<> > vse3W2C;
  TooN::Matrix<3> m3H;
  std::vector<TooN::Vector<3> > vv3VisiblePoints;
  std::vector<TooN::Vector<3> > vv3Templates;

  int iObject_ID;
  int iView_ID;
  std::vector<Edgelet> vEdgelets;
};*/

/*
  For a pair of edgelets, the relative direction, distance and angle are stored
*/ 
struct TestDirection
{
  int iIdx_0; // the index of the first edgelet in the pair
  int iIdx_1; // the index of the second edgelet in the pair
  REAL_TYPE rDir_1; // the angle between the first edgelet's orientation and the vector from the first edgelet to the second edgelet
  REAL_TYPE rDir_2; // the angle between the vector between the edgelets and the horizontal line
  REAL_TYPE rDist;  // the distance (L2) between the edgelets
  REAL_TYPE rAngle; // the relative angle between the orientations of both edgelets
};

/*
  For an edgelet and a cell centre, calculate the direction
*/ 
struct CellDirection
{
  int iIdx; // the index of the first edgelet in the pair
  int cellNo; // the index of the second edgelet in the pair
  REAL_TYPE rDir_1; // the angle between the edgelet's orientation and the vector from the first edgelet to the cell centre
  REAL_TYPE rDir_2; // the angle between the vector between the edgelet and the cell centre and the horizontal line
};

/*
  The structure for each line/element in the code book
*/
struct CodeBookElement
{
  REAL_TYPE rAngle_0; // the relative angle between edgelets 1 and 2
  REAL_TYPE rAngle_1; // the relative angle between edgelets 2 and 3
  REAL_TYPE rRelativeDist_1; // the relative distance between vectors 2-3 and 1-2
  REAL_TYPE rAngle_2; // the relative angle between edgelets 3 and 4
  REAL_TYPE rRelativeDist_2; // the relative distance between vectors 3-4 and 2-3
  REAL_TYPE rAngle_3; // the relative angle between edgelets 4 and 5
  REAL_TYPE rRelativeDist_3; // the relative distance between vectors 4-5 and 3-4
  unsigned short iIdx_0;  // the index of edgelet 1
  unsigned short iIdx_1;  // the index of edgelet 2
  unsigned short iIdx_2;  // the index of edgelet 3
  unsigned short iIdx_3;  // the index of edgelet 4
  unsigned short iIdx_4;  // the index of edgelet 5
  unsigned short iView_ID;
  unsigned short iObject_ID;
  unsigned short obj_correct_ID;
};

/*
  This structure stores the cells in a certain direction from a cell in the edge cell matrix
*/
struct CellsInDirection
{
   std::vector<unsigned int> cells;  	// list of cell numbers
   int cellNo;				// number of cell
   int dir;				// binned direction
};

/*
   A list of cells in direction for all directions of the same cell
*/
struct CellDirectStruct
{
   CellsInDirection cellVals [64];
};

/*
  At run time, the edgelets are binned into cells. For each cell, the cellNo, count and cell edgelets are recorded
*/
struct EdgeCellMatrix
{
   int cellNo;
   int count;
   TooN::Vector<2> v2Pose;
   std::vector<int> cellEdgelets;
};

/*
  For each pose/view of the object
*/
struct View
{
  View() : iView_ID(-1), iObject_ID(-1) {};
  int iView_ID;  		// number of the view
  int iObject_ID;		// id of the object
  std::vector<Edgelet> vEdgelets; // set of edgelets in this view
  Edges allEdges;		  // set of edges (image references) ??? 
  std::vector<int> viLineIndexes; // point to the first edgelet of lines ???
};

/*
  For each object
*/
struct ObjectElement
{
  ObjectElement() : iObject_ID(-1), iNoViews(0) {};
  int iObject_ID;		// Id of the object - incremental as loaded from the codebook
  int obj_correct_ID;		// correct id of the object - for visualisation purposes (this is fixed regardless of the number of objects loaded from the codebook or their order)
  int iNoViews; 		// number of views of that object
  std::vector<View> vViews; 	// list of views
};

/*
  For all objects
*/
struct ObjectsTemplate
{
  ObjectsTemplate() : iNoObjects(0) {};
  int iNoObjects; 		// number of loaded objects
  std::vector<ObjectElement> vObjectElements; // list of all objects
};

/*
  For each searched fixed path 
*/
struct ObjDescriptorClass 
{
  TooN::Vector<4> mv4ChainAngles;		// the angles of the path (4 for a path of 5 edgelets)
  ObjectsTemplate mObjectsTemplate;		// the list of objects to be searched
  std::vector<CodeBookElement> mvCodeBook;	// the loaded codebook
  std::vector<CVD::ImageRef> mvirFirstAngleIndexes;  //  hash-table into the codebook. ??? why image-ref?
  std::vector< std::vector<int> > mvirSecondIndexes; // has for second value (angle, angle, distance)
};


struct DetectedObject
{
  DetectedObject() : iObject_ID(-1) {};
  //std::vector<TooN::SE3<> > vse3W2C;
  std::vector<TooN::Vector<2> > vv2Pnts;
  std::vector<TooN::Vector<2> > foundChainPnts;
  CVD::ImageRef irTopLeft;
  CVD::ImageRef irBottomRight;
  double pose_x;
  double pose_y;
  double pose_z;
  double pose_w;
  double detection_time;
  double error;
  int iObject_ID;
  int iView_ID;
  int iObject_No;
  int iRegion_No;
  int missingCount;
};

struct SegmentedMask
{
  sensor_msgs::Image image;
};

struct SegmentedRegion
{
  CVD::ImageRef irTopLeft;
  CVD::ImageRef irBottomRight;
};

struct gtItem
{
  int frameNo;
  int objNo;
  int x1;
  int y1;
  int x2;
  int y2;
};

struct CheckedChains
{
  std::vector< TooN::Vector<2> > chains;
};
#endif
