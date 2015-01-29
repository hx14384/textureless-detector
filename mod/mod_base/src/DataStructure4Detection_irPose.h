#ifndef DATASTRUCTURE4DETECTION_IRPOSE_H
#define DATASTRUCTURE4DETECTION_IRPOSE_H

#include "DefineTypes.h"
#include <TooN/se3.h>
#include <TooN/TooN.h>
#include <cvd/image.h>
#include <cvd/image_ref.h>
#include <vector>

#define NO_TRAIN_QUESTIONS 3000
#define NO_QUESTIONS 500
#define QUESTIONS_SIZE 4
#define NO_MAX_BINS 64
#define BIN_RANGE 2
#define STEP 0.05
#define INV_STEP 20
// (0 ... 3.2)/0.05 , index in rDir(64) x rAngle(64)

// Contain a long line segment;
struct Line
{
  TooN::Vector<2> v2Start;      // starting point of the line 
  TooN::Vector<2> v2End;        // ending point of the line
  TooN::Vector<2> v2Direction;  // a unit direction pointing from v2Start to v2End
  REAL_TYPE rSlope;
  int no_pixel;
};

// Edgelet is a short line segment;
// Keep the centre and slope of the line;
struct Edgelet
{
//  TooN::Vector<2> v2Pose; // do ir(v2pose) before saving to file
//  REAL_TYPE rAngle;
  CVD::ImageRef irPose;
  REAL_TYPE rSlope;
};

struct EdgeLink
{
  TooN::Vector<2> v2Pose; // snd edge
  REAL_TYPE rDir; // may not need
  REAL_TYPE rAngle; // may not need
  REAL_TYPE rDistance;
  TooN::Vector<2> v2Check; // used for checking a side of point
};

// Put EdgeLinks generated from Edgelet into irBins.
struct BinsOfEdgeLinks 
{
  BinsOfEdgeLinks() : mbFound(false) { memset(mBitBins, 1, NO_MAX_BINS*NO_MAX_BINS*sizeof(bool)); };
  TooN::Vector<2> v2Pose; // fst edge
  bool mBitBins[NO_MAX_BINS][NO_MAX_BINS];
  std::vector<EdgeLink> mEdgeLinkBins[NO_MAX_BINS][NO_MAX_BINS];
  bool mbFound;
};

struct UCCoordinate
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
};

struct TestDirection
{
  int iIdx_0;
  int iIdx_1;
  REAL_TYPE rDir_1;
  REAL_TYPE rDir_2;
  REAL_TYPE rDist;
  REAL_TYPE rAngle;
};

struct CodeBookElement
{
  REAL_TYPE rAngle_0;
  REAL_TYPE rAngle_1;
  REAL_TYPE rRelativeDist_1;
  REAL_TYPE rAngle_2;
  REAL_TYPE rRelativeDist_2;
  REAL_TYPE rAngle_3;
  REAL_TYPE rRelativeDist_3;
  unsigned short iIdx_0;
  unsigned short iIdx_1;
  unsigned short iIdx_2;
  unsigned short iIdx_3;
  unsigned short iIdx_4;
  unsigned short iView_ID;
  unsigned short iObject_ID;
};

struct View
{
  View() : iView_ID(-1), iObject_ID(-1) {};
  int iView_ID;
  int iObject_ID;
  std::vector<Edgelet> vEdgelets;
  std::vector<int> viLineIndexes; // point to the first edgelet of lines
};

struct Object
{
  Object() : iObject_ID(-1), iNoViews(0) {};
  int iObject_ID;
  int iNoViews; // size of vViews;
  std::vector<View> vViews;
};

struct ObjectsTemplate
{
  ObjectsTemplate() : iNoObjects(0) {};
  int iNoObjects; // size of vObjects;
  std::vector<Object> vObjects;
};

struct ObjDescriptorClass 
{
  TooN::Vector<4> mv4ChainAngles;
  ObjectsTemplate mObjectsTemplate;
  std::vector<CodeBookElement> mvCodeBook;
  std::vector<CVD::ImageRef> mvirFirstAngleIndexes;
};
#endif
