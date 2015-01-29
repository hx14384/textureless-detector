#ifndef GETTRANSFORMDIST_H
#define GETTRANSFORMDIST_H
#include <assert.h>
#include <cvd/image.h>
#include "DefineTypes.h"
#define MY_MIN(a,b) ((a<b) ? a : b )

struct GetTransformDist
{
  inline static int compute( CVD::Image<REAL_TYPE> & imTransformDist, const CVD::Image<CVD::byte> & image )
  {
    assert( image.size() == imTransformDist.size() );
    int imDist [image.size().y] [image.size().x];
    int imDist2 [image.size().y] [image.size().x];
    //CVD::Image<REAL_TYPE> imDist( image.size(), 0 );
    int width = image.size().x;
    int height = image.size().y;
    int widthPlusheight = width+height;
    for( int row = 0; row < height; ++row )
    {
      for( int col = 0; col < width; ++col )
      {
       //   std::cerr << "Before adding " << (int) image[row][col] << std::endl;
          if ((int)image[row][col] > 100)
             imDist[row][col] = 0;
          else
          {
             imDist[row][col] = ( row == 0 ) ? ( col == 0 ) ? widthPlusheight : imDist[row][col-1]+1 : ( col == 0 ) ? imDist[row-1][col]+1 : MY_MIN( imDist[row-1][col]+1, imDist[row][col-1]+1 );
          }
          
      }
    }
    
    height -= 1; width -= 1;
    int maxValue = 0;
    for( int row = height; row > -1; --row )
    {
      for( int col = width; col > -1; --col )
      {
        imDist2[row][col] = ( col == width ) ? ( row == height ) ? imDist[row][col] : MY_MIN( imDist[row][col], imDist2[row+1][col]+1 ) : (row == height ) ? MY_MIN( imDist[row][col], imDist2[row][col+1]+1 ) : MY_MIN( MY_MIN( imDist[row][col], imDist2[row][col+1]+1 ), imDist2[row+1][col]+1 );
        if (imDist2[row][col] > maxValue)
          maxValue = imDist2[row][col];
      }
    }

    imTransformDist.fill(0);

    for( int row = height; row > -1; --row )
    {
      for( int col = width; col > -1; --col )
      {
        imTransformDist[row][col] = (REAL_TYPE)imDist2[row][col]/(float)maxValue;
      }
    }    
    return maxValue;
  };


inline static void integral( CVD::Image<REAL_TYPE> & integralImage, const CVD::Image<REAL_TYPE> & inputImage )
  {
    assert( inputImage.size() == integralImage.size() );
    //CVD::Image<REAL_TYPE> imDist( image.size(), 0 );
    int width = inputImage.size().x;
    int height = inputImage.size().y;
    for( int row = 0; row < height; ++row )
    {
      for( int col = 0; col < width; ++col )
      {
         integralImage[row][col] = ( row == 0 ) ? ( col == 0 ) ? inputImage[row][col] : inputImage[row][col]+integralImage[row][col-1] : ( col == 0 ) ? inputImage[row][col]+integralImage[row-1][col] : inputImage[row][col]+integralImage[row-1][col]+integralImage[row][col-1]-integralImage[row-1][col-1];
          
      }
    }
  };

inline static void integral( CVD::Image<REAL_TYPE> & integralImage, const CVD::Image<CVD::byte> & inputImage )
  {
    assert( inputImage.size() == integralImage.size() );
    //CVD::Image<REAL_TYPE> imDist( image.size(), 0 );
    int width = inputImage.size().x;
    int height = inputImage.size().y;
    for( int row = 0; row < height; ++row )
    {
      for( int col = 0; col < width; ++col )
      {
         integralImage[row][col] = ( row == 0 ) ? ( col == 0 ) ? inputImage[row][col] : inputImage[row][col]+integralImage[row][col-1] : ( col == 0 ) ? inputImage[row][col]+integralImage[row-1][col] : inputImage[row][col]+integralImage[row-1][col]+integralImage[row][col-1]-integralImage[row-1][col-1];
          
      }
    }
  };

  inline static void compute( Orientation & transformDist, Orientation & image)
  {
    assert( image.size() == transformDist.size() );
    Orientation imDist;
    CVD::ImageRef irSize1;
    irSize1.x = image.size().x;
    irSize1.y = image.size().y;
    imDist.resize(irSize1);
    imDist.reset();
    int width = imDist.size().x;
    int height = imDist.size().y;
    int widthPlusheight = width+height;
    for( int row = 0; row < height; ++row )
    {
      for( int col = 0; col < width; ++col )
      {
        int idx = row*width+col;
        int idx_1 = (row-1)*width+col;
        if( image[idx] == 0 )
          imDist[idx] = ( row == 0 ) ? ( col == 0 ) ? widthPlusheight : imDist[idx-1]+1 : ( col == 0 ) ? imDist[idx_1]+1 : MY_MIN( imDist[idx_1]+1, imDist[idx-1]+1 );
      }
    }
    transformDist.reset();
    height -= 1; width -= 1;
    for( int row = height; row > -1; --row )
    {
      for( int col = width; col > -1; --col )
      {
        int idx = row*(width+1)+col;
        int idx_1 = (row+1)*(width+1)+col;
        transformDist[idx] = ( col == width ) ? ( row == height ) ? imDist[idx] : MY_MIN( imDist[idx], transformDist[idx_1]+1 ) : (row == height ) ? MY_MIN( imDist[idx], transformDist[idx+1]+1 ) : MY_MIN( MY_MIN( imDist[idx], transformDist[idx+1]+1 ), transformDist[idx_1]+1 );
      }
    }
  };
};
#endif

