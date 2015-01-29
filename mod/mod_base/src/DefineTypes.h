/*
  Computer Vision Group, Department of Computer Science. University of      
  Bristol. This code can not be used or copied without express permission.   
*/

#ifndef DEFINE_TYPES_H
#define DEFINE_TYPES_H

#define REAL_TYPE float
//#define REAL_TYPE double
#include <TooN/TooN.h>
#include <vector>
#include <cvd/camera.h>
#include <gvars3/instances.h>
#include <TooN/se3.h>
#include <cvd/image.h>
#include <cvd/image_io.h>
#include "LowLevelImageData.h"
#include "LLImagePyramid.h"
#include <string.h>

struct BasicPlane
{
  std::vector<TooN::Vector<3> > vv3poses;
  TooN::Vector<4> v4Plane;
  void clear() { vv3poses.clear(); v4Plane = TooN::Zeros; };
  BasicPlane & operator =( const BasicPlane & rhs )
  {
    if( this != &rhs )
    {
      vv3poses = rhs.vv3poses;
      v4Plane = rhs.v4Plane;
    }
    return *this;
  }
  ;
};

/*struct CameraView
{
  CVD::Image<REAL_TYPE> * grey;
  CVD::Image<CVD::Rgb<CVD::byte> > * rgb;
  TooN::SE3<> se3W2C;
  CameraView() :
    grey( NULL ), rgb( NULL )
  {
  }
  ;
  ~CameraView()
  {
    if( grey != NULL )
    {
      delete grey;
      delete rgb;
    }
  }
  ;
};*/

/*template <typename T>
struct ImageView
{
  LowLevelImageData<T> llImage;
  TooN::SE3<> se3W2C;
};*/

/* 
  List of all edgelets relative to their image references 
*/
struct Edges
{
  void clear() { virEdges.clear(); viEdgeIdxes.clear(); } ;
  std::vector<CVD::ImageRef> virEdges;
  std::vector<int> viEdgeIdxes;
};

/*
  Image of oriented edges. For each binned orientation (fixed to 11 here), an image of the same size is used for each bin
  Pixels that correspond to edges in each orientation are copied into the image of the specified orientation
  This is used to create oriented distance transform images
*/
struct Orientation
{
  Orientation(): irSize(320, 240) {
    data = NULL;
  };
  ~Orientation(){
     free (data);
  };
  CVD::ImageRef irSize; // size of the image
  REAL_TYPE *data;	// pointer to data
  void resize(CVD::ImageRef irSize1)	// resizes the pointer based on the input image size
  {  irSize.x = irSize1.x;
     irSize.y = irSize1.y;
//     void* tmp = realloc (data, 11*irSize.x*irSize.y*sizeof(REAL_TYPE));
     void* tmp = realloc (data, irSize.x*irSize.y*sizeof(REAL_TYPE));
     if (tmp == NULL)
         std::cout << " CANNOT REALLOCATE MEMORY " << std::endl;
     data = (REAL_TYPE*) tmp;
  };

  inline void reset () 			// resets the memory to zero for the next frame
  {
//    memset(data, 0, 11*irSize.x*irSize.y*sizeof(REAL_TYPE));
    memset(data, 0, irSize.x*irSize.y*sizeof(REAL_TYPE));
  };

  inline REAL_TYPE & operator[](const TooN::Vector<2>& v2Pose)
  {
//    return data[static_cast<int>(v2Pose[0])+static_cast<int>(v2Pose[1])*irSize.y];
    return data[static_cast<int>(v2Pose[0])+static_cast<int>(v2Pose[1])*irSize.x];
  };

  inline REAL_TYPE & operator[](int idx)
  {
    return data[idx];
  };

  inline CVD::ImageRef & size() { return irSize; }; // returns the size
};


// Undistorted camera model ??? it assumes the images are distorted, right???
extern Camera::Linear gUndistortedCameraModel;

// Distorted camera model
extern Camera::Quintic gQuinticCameraModel;

// Camera parameters
extern GVars3::gvar3<TooN::Vector<6> > ggvv6CameraParams;

#endif

