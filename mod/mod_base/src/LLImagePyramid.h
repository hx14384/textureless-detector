/*
  Computer Vision Group, Department of Computer Science. University of      
  Bristol. This code can not be used or copied without express permission.   
*/

#ifndef LLIMAGEPYRAMID_H
#define LLIMAGEPYRAMID_H

#include "LowLevelImageData.h"
#include <vector>
#include <cmath>
/**
 @file LLImagePyramid.h
 @brief LLImagePyramid will handle an image in a different scale. It will keep image in each scale in 
  LowLevelImageData class.
*/

template <typename T>
struct LLImagePyramid
{
  public:
    /**
      @param[in] irSizeAtLevelZero is image size at level zero.
      @param[in] no_levels is a number of image scale.
    */
    LLImagePyramid( CVD::ImageRef irSizeAtLevelZero, unsigned int no_levels = 3 );
    /**
      @param[in] image is a master image at level zero.
      @param[in] no_levels is a number of image scale.
    */
    LLImagePyramid( const CVD::Image<T> & image, unsigned int no_levels = 3 );
    /**
      @param[in] image is a master image at level zero.
    */
    void SetImageAtLevelZero( const CVD::Image<T> & image );
    void SetImageAtLevelZeroOnly (const CVD::Image<T> & image);
    LowLevelImageData<T> & GetImage( unsigned int no_level );
    const LowLevelImageData<T> & GetImage( unsigned int no_level ) const;
    void Reset();
  private:
    unsigned int mNoLevels;
    std::vector<LowLevelImageData<T> > mvImage;
};

template <typename T>
LLImagePyramid<T>::LLImagePyramid( CVD::ImageRef irSizeAtLevelZero, unsigned int no_levels )
{
  mNoLevels = no_levels;
  mvImage.resize( no_levels );
  int devide = 1;
  for( unsigned int i = 0; i < no_levels; ++i )
  {
    mvImage[i].SetImageSize( irSizeAtLevelZero/devide );
    devide *= 2; 
  }
};

template <typename T>
LLImagePyramid<T>::LLImagePyramid( const CVD::Image<T> & image, unsigned int no_levels )
{
  mNoLevels = no_levels;
  mvImage.resize( no_levels );
  mvImage[0].SetImageSize( image.size() );
  mvImage[0].SetImage( image );
  for( unsigned int i = 1; i < mNoLevels; ++i )
  {
    mvImage[i].SetImageSize( mvImage[i-1].mImage.size()/2 );
    mvImage[i].SetImage( mvImage[i-1].mImage, true );
  }
};

template <typename T>
void LLImagePyramid<T>::Reset()
{
  for( unsigned int i = 0; i < mNoLevels; ++i )
    mvImage[i].Reset();
}

template <typename T>
void
LLImagePyramid<T>::SetImageAtLevelZero( const CVD::Image<T> & image )
{
  mvImage[0].SetImage( image );
  for( unsigned int i = 1; i < mNoLevels; ++i )
  {
    mvImage[i].SetImage( mvImage[i-1].mImage, true );
  }
};

template <typename T>
void
LLImagePyramid<T>::SetImageAtLevelZeroOnly( const CVD::Image<T> & image )
{
  mvImage[0].SetImage( image );
};

template <typename T>
LowLevelImageData<T> &
LLImagePyramid<T>::GetImage( unsigned int no_level )
{
  assert( no_level < mNoLevels );
  return mvImage[no_level];
};

template <typename T>
const LowLevelImageData<T> &
LLImagePyramid<T>::GetImage( unsigned int no_level ) const
{
  assert( no_level < mNoLevels );
  return mvImage[no_level];
};

#endif
