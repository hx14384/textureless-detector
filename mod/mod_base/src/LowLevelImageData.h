/*
  Computer Vision Group, Department of Computer Science. University of      
  Bristol. This code can not be used or copied without express permission.   
*/

#ifndef LOWLEVELIMAGEDATA_H
#define LOWLEVELIMAGEDATA_H

#include <cvd/image.h>
#include <cvd/image_ref.h>
#include <cvd/utility.h>
#include <cvd/fast_corner.h>
#include <cvd/integral_image.h>
#include <cvd/vision.h>
#include <cvd/convolution.h>
#include <vector>
#include <algorithm>
#include "DefineTypes.h"
#include <cvd/colourspace_convert.h>

template<typename T>
class LowLevelImageData
{
  public:
    enum CornerState
    {
      NONE, FAST_CORNER, MAX_CORNER, BEST_CORNER, BEST_CORNER_IN_GRID
    };

    LowLevelImageData() :
      mbComputedIntegralImage( false ), mbComputedSmoothedImage( false ), mbComputedGradientImage( false ), mbComputedGradMagnitudeImage( false), 
        mCornerState( NONE )
    {
    }
    ;
    LowLevelImageData( CVD::ImageRef irSize ) :
      mImage( irSize ), mIntegralImage( irSize ), mGradientImage( irSize ), mSmoothedImage( irSize), mGradMagnitudeImage( irSize ),
        mbComputedIntegralImage( false ), mbComputedSmoothedImage( false ), mbComputedGradientImage( false ), mbComputedGradMagnitudeImage( false ), 
          mCornerState( NONE ), mvRowIdx( irSize.y )
    {
    }
    ;
    LowLevelImageData( const CVD::Image<T> & image ) :
      mImage( image.size() ), mIntegralImage( image.size() ), mGradientImage( image.size() ), mSmoothedImage( image.size() ),
        mGradMagnitudeImage( image.size() ), mbComputedIntegralImage( false ), mbComputedSmoothedImage( false ), mbComputedGradientImage( false ), 
          mbComputedGradMagnitudeImage( false ), mbNormalizedGradMagnitudeImage( false), mCornerState( NONE ), mvRowIdx( image.size().y )
    {
      SetImage( image );
    }
    ;
    void SetImageSize( CVD::ImageRef irSize );
    void SetImage( const CVD::Image<T> & image );
    void SetImage( const CVD::Image<T> & image, bool isHalfSample );
    void Reset();
    void ComputeIntegralImage();
    void ComputeSmoothedImage( REAL_TYPE sigma = 1.6f, bool bRealDataType = true );
    void ComputeGradientImage( REAL_TYPE sigma = -1.f );
    void ComputeGradMagnitudeImage( REAL_TYPE sigma = -1.f );
    void NormalizeGradMagnitudeImage();
    void ComputeFastCorner10( int barrier = 10 );
    void ComputeFastCorner9( int barrier = 10 );
    void ComputeMaxCorners( int barrier = 10 );
    void ComputeBestCorners( unsigned int no_best_corners = 1000, int barrier = 10, REAL_TYPE min_score = 50.0 );
    void ComputeBestCornersInGrids( int width = 20, int height = 15, int barrier = 10, REAL_TYPE min_score = 50.0 );
    void ComputeGroupOfBestCornersInGrids( int width = 80, int height = 60, int max_points = 100, int barrier = 10, REAL_TYPE min_score = 50.0 );
    void ComputeRowIndex();

    LowLevelImageData<T> & operator=( const LowLevelImageData<T> & rhs )
    {
      if( this != &rhs )
      {
        mImage = rhs.mImage; mImage.make_unique();
        mIntegralImage = rhs.mIntegralImage; mIntegralImage.make_unique();
#if __GNUC__ == 4 && __GNUC_MINOR__ > 2
        mGradientImage = rhs.mGradientImage; mGradientImage.make_unique();
        mbComputedGradientImage = rhs.mbComputedGradientImage;
#else // there is a bug for gcc 4.2. Have to recompute its gradient.
        mGradientImage.resize( rhs.mGradientImage.size() );
        mbComputedGradientImage = false;
#endif
        mSmoothedImage = rhs.mSmoothedImage; mSmoothedImage.make_unique();
        mGradMagnitudeImage = rhs.mGradMagnitudeImage; mGradMagnitudeImage.make_unique();
        mbComputedIntegralImage = rhs.mbComputedIntegralImage;
        mbComputedSmoothedImage = rhs.mbComputedSmoothedImage;
        mbComputedGradMagnitudeImage = rhs.mbComputedGradMagnitudeImage;
        mbNormalizedGradMagnitudeImage = rhs.mbNormalizedGradMagnitudeImage;
        mvFastCorners = rhs.mvFastCorners;
        mvMaxCorners = rhs.mvMaxCorners;
        mvBestCorners = rhs.mvBestCorners;
        mvRowIdx = rhs.mvRowIdx;
        mCornerState = rhs.mCornerState;
        miMaxWidth = rhs.miMaxWidth;
        miMaxHeight = rhs.miMaxHeight;
      }
      return *this;
    }
    ;

    CVD::Image<T> mImage;
    CVD::Image<REAL_TYPE> mIntegralImage;
    CVD::Image<REAL_TYPE[2]> mGradientImage;
    CVD::Image<REAL_TYPE> mSmoothedImage;
    CVD::Image<REAL_TYPE> mGradMagnitudeImage;
    bool mbComputedIntegralImage;
    bool mbComputedSmoothedImage;
    bool mbComputedGradientImage;
    bool mbComputedGradMagnitudeImage;
    bool mbNormalizedGradMagnitudeImage;
    REAL_TYPE mrMaxGradMagnitude;
    std::vector<CVD::ImageRef> mvFastCorners;
    std::vector<CVD::ImageRef> mvMaxCorners;
    std::vector<CVD::ImageRef> mvBestCorners;
    CornerState mCornerState;
    std::vector<unsigned int> mvRowIdx;
  private:
    int miMaxWidth;
    int miMaxHeight;
    static const int miBorder = 5;
};

template<typename T>
void LowLevelImageData<T>::SetImageSize( CVD::ImageRef irSize )
{
  mImage.resize( irSize );
  mIntegralImage.resize( irSize );
  mGradientImage.resize( irSize );
  mSmoothedImage.resize( irSize );
  mGradMagnitudeImage.resize( irSize );
  miMaxWidth = irSize.x - miBorder;
  miMaxHeight = irSize.y - miBorder;
  mvRowIdx.resize( irSize.y );
}
;

template<typename T>
void LowLevelImageData<T>::SetImage( const CVD::Image<T> & image )
{
  mImage = image;
  mImage.make_unique();
  Reset();
}
;

template<typename T>
void LowLevelImageData<T>::SetImage( const CVD::Image<T> & image, bool isHalfSample )
{
  if( isHalfSample )
  {
    CVD::halfSample( image, mImage );
  }
  else
  {
    mImage = image;
    mImage.make_unique();
  }
  Reset();
}
;

template<typename T>
void LowLevelImageData<T>::Reset()
{
  mvFastCorners.clear();
  mvMaxCorners.clear();
  mvBestCorners.clear();
  mbComputedIntegralImage = false;
  mbComputedSmoothedImage = false;
  mbComputedGradientImage = false;
  mbComputedGradMagnitudeImage = false;
  mbNormalizedGradMagnitudeImage = false;
  mCornerState = NONE;
}
;

template<typename T>
void LowLevelImageData<T>::ComputeGradMagnitudeImage( REAL_TYPE sigma )
{
  if( mbComputedGradMagnitudeImage ) return;
  ComputeGradientImage( sigma );
  mrMaxGradMagnitude = 0.f;
  int w = mGradMagnitudeImage.size().x;
  CVD::BasicImage<REAL_TYPE[2]>::const_iterator beginIterator = mGradientImage.begin() + w + 1; // not include border
  CVD::BasicImage<REAL_TYPE[2]>::const_iterator endIterator = mGradientImage.end() - w - 1; 
  CVD::BasicImage<REAL_TYPE>::iterator resultIterator = mGradMagnitudeImage.begin() + w + 1;
  while( beginIterator != endIterator )
  {
    *resultIterator = CVD::abs((*beginIterator)[0]) + CVD::abs((*beginIterator)[1]);
    if( *resultIterator > mrMaxGradMagnitude )
      mrMaxGradMagnitude = *resultIterator;
    resultIterator++;
    beginIterator++;
  }
  mbComputedGradMagnitudeImage = true;
}
;

template<typename T>
void LowLevelImageData<T>::NormalizeGradMagnitudeImage( )
{
  if( mbNormalizedGradMagnitudeImage ) return;
  CVD::BasicImage<REAL_TYPE>::iterator beginIterator = mGradMagnitudeImage.begin() + mGradMagnitudeImage.size().x + 1;
  CVD::BasicImage<REAL_TYPE>::iterator endIterator = mGradMagnitudeImage.end() - mGradMagnitudeImage.size().x - 1;
  REAL_TYPE invMaxGradMagnitude = 1.f/mrMaxGradMagnitude;
  while( beginIterator != endIterator )
  {
    *beginIterator *= invMaxGradMagnitude;
    ++beginIterator;
  }
  mbNormalizedGradMagnitudeImage = true;
}
;

template<typename T>
void LowLevelImageData<T>::ComputeIntegralImage()
{
  if( mbComputedIntegralImage ) return;
  CVD::integral_image( mImage, mIntegralImage );
  mbComputedIntegralImage = true;
}
;

template<typename T>
void LowLevelImageData<T>::ComputeSmoothedImage( REAL_TYPE sigma, bool bRealDataType )
{
  if( mbComputedSmoothedImage ) return;
  if( bRealDataType )
  {
    mSmoothedImage = convert_image( mImage );
    CVD::convolveGaussian( mSmoothedImage, sigma );
  }
  else
    CVD::convolveGaussian( mImage, sigma );
  mbComputedSmoothedImage = true;
}
;

template<typename T>
void LowLevelImageData<T>::ComputeGradientImage( REAL_TYPE sigma )
{
  if( mbComputedGradientImage ) return;
  if( sigma > 0 )
  {
    ComputeSmoothedImage( sigma );
    CVD::gradient( mSmoothedImage, mGradientImage );
  }
  else
  {
    CVD::gradient( mImage, mGradientImage );
  }
  mbComputedGradientImage = true;
}
;

template<typename T>
void LowLevelImageData<T>::ComputeRowIndex()
{
  unsigned int idx = 0;
  std::vector<unsigned int>::iterator vRowIdxIter = mvRowIdx.begin();
  switch( mCornerState )
  {
    case FAST_CORNER:
      for( int y = 0; y < mImage.size().y; ++y )
      {
        while( y > mvFastCorners[idx].y )
        {
          ++idx;
          if( idx == mvFastCorners.size() )
          {
            std::fill(vRowIdxIter, mvRowIdx.end(), 0);
            goto STOP;
          }
        }
        *vRowIdxIter++ = idx;
//        mvRowIdx.push_back( idx );
      }
      break;
    case MAX_CORNER:
      for( int y = 0; y < mImage.size().y; ++y )
      {
        while( y > mvMaxCorners[idx].y )
        {
          ++idx;
          if( idx == mvMaxCorners.size() )
          {
            std::fill(vRowIdxIter, mvRowIdx.end(), 0);
            goto STOP;
          }
        }
        *vRowIdxIter++ = idx;
//        mvRowIdx.push_back( idx );
      }
      break;
    default:
      break;
  }
  STOP:;
}
;

template<typename T>
void LowLevelImageData<T>::ComputeFastCorner10( int barrier )
{
  if( mCornerState == FAST_CORNER ) return;
  CVD::fast_corner_detect_10( mImage, mvFastCorners, barrier );
  mCornerState = FAST_CORNER;
}
;

template<typename T>
void LowLevelImageData<T>::ComputeFastCorner9( int barrier )
{
  if( mCornerState == FAST_CORNER ) return;
  CVD::fast_corner_detect_9( mImage, mvFastCorners, barrier );
  mCornerState = FAST_CORNER;
}
;

template<typename T>
void LowLevelImageData<T>::ComputeMaxCorners( int barrier )
{
  if( mCornerState == MAX_CORNER ) return;
  ComputeFastCorner9( barrier ); // if ComputeFastCorner 9 or 10 haven't called, call ComputeFastCorner9 as a default.
  CVD::fast_nonmax( mImage, mvFastCorners, barrier, mvMaxCorners );
  mCornerState = MAX_CORNER;
}
;

template<typename T>
void LowLevelImageData<T>::ComputeBestCorners( unsigned int no_best_corners, int barrier, REAL_TYPE min_score )
{
  if( mCornerState == BEST_CORNER ) return;
  std::vector<std::pair<REAL_TYPE, CVD::ImageRef> > vScoreAndCorners;
  ComputeMaxCorners( barrier );
  for( std::vector<CVD::ImageRef>::iterator it = mvMaxCorners.begin(); it != mvMaxCorners.end(); it++ )
  {
    CVD::ImageRef & ir = *it;
    if( ir.x < miBorder || ir.x >= miMaxWidth || ir.y < miBorder || ir.y >= miMaxHeight ) continue;
    REAL_TYPE score = CalculateTrackingScore( mImage, ir );
    if( score > min_score ) vScoreAndCorners.push_back( std::make_pair( -1.f * score, ir ) );
  }
  std::sort( vScoreAndCorners.begin(), vScoreAndCorners.end() );
  if( vScoreAndCorners.size() > no_best_corners ) vScoreAndCorners.resize( no_best_corners );
  for( std::vector<std::pair<REAL_TYPE, CVD::ImageRef> >::iterator it = vScoreAndCorners.begin(); it
      != vScoreAndCorners.end(); it++ )
    mvBestCorners.push_back( ( *it ).second );
  mCornerState = BEST_CORNER;
}
;

template<typename T>
void LowLevelImageData<T>::ComputeBestCornersInGrids( int width, int height, int barrier, REAL_TYPE min_score )// width of grid in pixel
{
  if( mCornerState == BEST_CORNER_IN_GRID ) return;
  ComputeMaxCorners( barrier );
  ComputeRowIndex();
  REAL_TYPE max_score;
  CVD::ImageRef irMax;
  int grid_x = mImage.size().x / width;
  int grid_y = mImage.size().y / height;
  int left, right;
// the last row of grid; j = 0 
  for( int i = 0; i < grid_x; ++i )
  {
    max_score = min_score;
    for( int row = miBorder; row < height; ++row )
    {
      if( mvRowIdx[ row ] == mvRowIdx[ row + 1 ] )
        continue;
      left = i*width ; right = (i+1)*width; 
      left = ( left < miBorder ) ? miBorder : left; right = ( right > miMaxWidth ) ? miMaxWidth : right;
      std::vector<CVD::ImageRef>::iterator irIter = mvMaxCorners.begin() + mvRowIdx[ row ];
      while( irIter != mvMaxCorners.end() && irIter->y == row && irIter->x < left )
        ++irIter;
      while( irIter != mvMaxCorners.end() && irIter->y == row && irIter->x < right )
      {
        REAL_TYPE score = CalculateTrackingScore( mImage, *irIter );
        if( score > max_score )
        {
          max_score = score;
          irMax = *irIter; 
        }
        ++irIter;
      }
    }
    if( max_score > min_score )
    {
      mvBestCorners.push_back( irMax );
    }
  }
// normal cases  
  for( int j = 1; j < grid_y - 1; ++j )
  {
    for( int i = 0; i < grid_x; ++i )
    {
      max_score = min_score;
      for( int row = 0; row < height; ++row )
      {
        int atRow = row + j*height;
        if( mvRowIdx[ atRow ] == mvRowIdx[ atRow + 1 ] )
          continue;
        left = i*width ; right = (i+1)*width; 
        left = ( left < miBorder ) ? miBorder : left; right = ( right > miMaxWidth ) ? miMaxWidth : right;
        std::vector<CVD::ImageRef>::iterator irIter = mvMaxCorners.begin() + mvRowIdx[ atRow ];
        while( irIter != mvMaxCorners.end() && irIter->y == atRow && irIter->x < left )
          ++irIter;
        while( irIter != mvMaxCorners.end() && irIter->y == atRow && irIter->x < right )
        {
          REAL_TYPE score = CalculateTrackingScore( mImage, *irIter );
          if( score > max_score )
          {
            max_score = score;
            irMax = *irIter; 
          }
          ++irIter;
        }
      }
      if( max_score > min_score )
      {
        mvBestCorners.push_back( irMax );
      }
    }
  }
// the last row of grid; j = grid_y - 1 
  for( int i = 0; i < grid_x; ++i )
  {
    max_score = min_score;
    for( int row = height*(grid_y - 1); row < miMaxHeight; ++row )
    {
      if( mvRowIdx[ row ] == mvRowIdx[ row + 1 ] )
        continue;
      left = i*width; right = (i+1)*width;
      left = ( left < miBorder ) ? miBorder : left; right = ( right > miMaxWidth ) ? miMaxWidth : right;
      std::vector<CVD::ImageRef>::iterator irIter = mvMaxCorners.begin() + mvRowIdx[ row ];
      while( irIter != mvMaxCorners.end() && irIter->y == row && irIter->x < left ) ++irIter;
      while( irIter != mvMaxCorners.end() && irIter->y == row && irIter->x < right )
      {
        REAL_TYPE score = CalculateTrackingScore( mImage, *irIter );
        if( score > max_score )
        {
          max_score = score;
          irMax = *irIter; 
        }
        ++irIter;
      }
    }
    if( max_score > min_score )
    {
      mvBestCorners.push_back( irMax );
    }
  }
  mCornerState = BEST_CORNER_IN_GRID;
}
;

template<typename T>
void LowLevelImageData<T>::ComputeGroupOfBestCornersInGrids( int width, int height, int max_points, int barrier, REAL_TYPE min_score )
{
  if( mCornerState == BEST_CORNER_IN_GRID ) return;
  std::vector<std::pair<REAL_TYPE, CVD::ImageRef> > vScoreAndCorners;
  ComputeMaxCorners( barrier );
  ComputeRowIndex();
  int grid_x = mImage.size().x / width;
  int grid_y = mImage.size().y / height;
  int left, right;
 // the last row of grid; j = 0 
  for( int i = 0; i < grid_x; ++i )
  {
    vScoreAndCorners.clear();
    for( int row = miBorder; row < height; ++row )
    {
      if( mvRowIdx[ row ] == mvRowIdx[ row + 1 ] )
        continue;
      left = i*width; right = (i+1)*width;
      left = ( left < miBorder ) ? miBorder : left; right = ( right > miMaxWidth ) ? miMaxWidth : right;
      std::vector<CVD::ImageRef>::iterator irIter = mvMaxCorners.begin() + mvRowIdx[ row ];
      while( irIter != mvMaxCorners.end() && irIter->y == row && irIter->x < left ) ++irIter;
      while( irIter != mvMaxCorners.end() && irIter->y == row && irIter->x < right )
      {
        REAL_TYPE score = CalculateTrackingScore( mImage, *irIter );
        if( score > min_score )
          vScoreAndCorners.push_back( std::make_pair(-1.f*score, *irIter) );
        ++irIter;
      }
    }
    if( vScoreAndCorners.size() > 0 )
    {
      std::sort( vScoreAndCorners.begin(), vScoreAndCorners.end() );
      if( vScoreAndCorners.size() > max_points )
        vScoreAndCorners.resize( max_points );
      for( std::vector<std::pair<REAL_TYPE, CVD::ImageRef> >::iterator itr = vScoreAndCorners.begin();
        itr != vScoreAndCorners.end(); ++itr )
          mvBestCorners.push_back( (*itr).second );
    }
  }
  
// normal cases  
  for( int j = 1; j < grid_y - 1; ++j )
  {
    for( int i = 0; i < grid_x; ++i )
    {
      vScoreAndCorners.clear();
      for( int row = 0; row < height; ++row )
      {
        int atRow = row + j*height;
        if( mvRowIdx[ atRow ] == mvRowIdx[ atRow + 1 ] )
          continue;
        left = i*width; right = (i+1)*width;
        left = ( left < miBorder ) ? miBorder : left; right = ( right > miMaxWidth ) ? miMaxWidth : right;
        std::vector<CVD::ImageRef>::iterator irIter = mvMaxCorners.begin() + mvRowIdx[ atRow ];
        while( irIter != mvMaxCorners.end() && irIter->y == atRow && irIter->x < left ) ++irIter;
        while( irIter != mvMaxCorners.end() && irIter->y == atRow && irIter->x < right )
        {
          REAL_TYPE score = CalculateTrackingScore( mImage, *irIter );
          if( score > min_score )
            vScoreAndCorners.push_back( std::make_pair(-1.f*score, *irIter) );
          ++irIter;
        }
      }
      if( vScoreAndCorners.size() > 0 )
      {
        std::sort( vScoreAndCorners.begin(), vScoreAndCorners.end() );
        if( vScoreAndCorners.size() > max_points )
          vScoreAndCorners.resize( max_points );
        for( std::vector<std::pair<REAL_TYPE, CVD::ImageRef> >::iterator itr = vScoreAndCorners.begin();
          itr != vScoreAndCorners.end(); ++itr )
            mvBestCorners.push_back( (*itr).second );
      }
    }
  }
// the last row of grid; j = grid_y - 1 
  for( int i = 0; i < grid_x; ++i )
  {
    vScoreAndCorners.clear();
    for( int row = height*(grid_y - 1); row < miMaxHeight; ++row )
    {
      if( mvRowIdx[ row ] == mvRowIdx[ row + 1 ] )
        continue;
      left = i*width; right = (i+1)*width;
      left = ( left < miBorder ) ? miBorder : left; right = ( right > miMaxWidth ) ? miMaxWidth : right;
      std::vector<CVD::ImageRef>::iterator irIter = mvMaxCorners.begin() + mvRowIdx[ row ];
      while( irIter != mvMaxCorners.end() && irIter->y == row && irIter->x < left ) ++irIter;
      while( irIter != mvMaxCorners.end() && irIter->y == row && irIter->x < right )
      {
        REAL_TYPE score = CalculateTrackingScore( mImage, *irIter );
        if( score > min_score )
          vScoreAndCorners.push_back( std::make_pair(-1.f*score, *irIter) );
        ++irIter;
      }
    }
    if( vScoreAndCorners.size() > 0 )
    {
      std::sort( vScoreAndCorners.begin(), vScoreAndCorners.end() );
      if( vScoreAndCorners.size() > max_points )
        vScoreAndCorners.resize( max_points );
      for( std::vector<std::pair<REAL_TYPE, CVD::ImageRef> >::iterator itr = vScoreAndCorners.begin();
        itr != vScoreAndCorners.end(); ++itr )
          mvBestCorners.push_back( (*itr).second );
    }
  }
  mCornerState = BEST_CORNER_IN_GRID;
}
;

#endif
