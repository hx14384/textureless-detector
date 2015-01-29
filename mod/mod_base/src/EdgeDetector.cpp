#include "EdgeDetector.h"
#include "BitMacros.h"

//#define _SOBEL_LIKE_
#define _SOBEL_
//#define _SYM_
#ifdef _SOBEL_LIKE_
const int MaskX[5][5] =
{
  { 1, 1, 0, -1, -1},
  { 2, 2, 0, -2, -2},
  { 3, 3, 0, -3, -3},
  { 2, 2, 0, -2, -2},
  { 1, 1, 0, -1, -1}
};

const int MaskY[5][5] =
{
  { 1, 2, 3, 2, 1},
  { 1, 2, 3, 2, 1},
  { 0, 0, 0, 0, 0},
  { -1, -2, -3, -2, -1},
  { -1, -2, -3, -2, -1}
};
#elif defined _SYM_
const int MaskX[5][5] =
{
  { 10, 10, 0, -10, -10},
  { 10, 10, 0, -10, -10},
  { 10, 10, 0, -10, -10},
  { 10, 10, 0, -10, -10},
  { 10, 10, 0, -10, -10}
};

const int MaskY[5][5] =
{
  { 10, 10, 10, 10, 10},
  { 10, 10, 10, 10, 10},
  { 0, 0, 0, 0, 0},
  { -10, -10, -10, -10, -10},
  { -10, -10, -10, -10, -10}
};
#elif defined _SOBEL_
const int MaskX[3][3] = { { 1, 0, -1 }, { 2, 0, -2 }, { 1, 0, -1 }, };
const int MaskY[3][3] = { { 1, 2, 1 }, { 0, 0, 0 }, { -1, -2, -1 } };
#endif

EdgeDetector::EdgeDetector( int iMaxSearchRange, REAL_TYPE fEdgeThreshold ): miMaxSearchRange(iMaxSearchRange), mfEdgeThreshold(fEdgeThreshold), mvrPMAX( 2*miMaxSearchRange + 1, 0)
{

}

EdgeDetector::~EdgeDetector()
{

}

void EdgeDetector::mark_edge_point_along_search_path( std::vector<int>& viEdgeIndex,
    std::vector<TooN::Vector<3> > & vv3EdgeStrs, const CVD::Image<REAL_TYPE>& image,
    const TooN::Vector<2>& v2ImagePoint, const TooN::Vector<2>& v2Normal )
{
  TooN::Vector<3> edgeStr;
  int i, j;
  register int k, chkpt, nSearch, nSearch_2nw, nSearch_nw, TwoiMaxSearchRange;
  TooN::Vector<2> v2Xn;
  viEdgeIndex.clear();
  TwoiMaxSearchRange = 2*miMaxSearchRange;
  // apply mask to find edgeStr in both magnitude and direction for a whole range of search path.
  for( nSearch = 0; nSearch <= TwoiMaxSearchRange; ++nSearch )
  {
    v2Xn = v2ImagePoint + ( nSearch - miMaxSearchRange ) * v2Normal;
    if( !image.in_image_with_border( CVD::ir( v2Xn ), 10 ) ) break;
    edgeStrength( edgeStr, image, int( v2Xn[0] ), int( v2Xn[1] ) );
    vv3EdgeStrs[nSearch] = edgeStr;
  }
  // 1D NMS for (2nw+1)-Neighborhood
  // Efficient Non-Maximum Suppression, Alexander Neubeck and Luc Van Gool, ICPR, 2006
  const int nw = 2; // test nw = 2, must be more than or equal to 2;
  i = nw;
  CompPartialMax( vv3EdgeStrs, 0, i - 1 );
  chkpt = -1;
  nSearch_2nw = nSearch - 2*nw; // nSearch may not be equal to TwoiMaxSearchRange if it is reach an image border.
  nSearch_nw = nSearch - nw;
  while( i <= nSearch_2nw )
  {
    j = CompPartialMax( vv3EdgeStrs, i, i + nw );
    k = CompPartialMax( vv3EdgeStrs, i + nw + 1, j + nw );
    if( i == j || vv3EdgeStrs[j][0] > vv3EdgeStrs[k][0] )
    {
      if( ( chkpt <= j - nw || vv3EdgeStrs[j][0] >= mvrPMAX[chkpt] )
          && ( j - nw == i || vv3EdgeStrs[j][0] >= mvrPMAX[j - nw] ) )
      {
        if( vv3EdgeStrs[j][0] > mfEdgeThreshold ) viEdgeIndex.push_back( j );
      }
      if( i < j ) chkpt = i + nw + 1;
      i = j + nw + 1;
    }
    else
    {
      i = k;
      chkpt = j + nw + 1;
      while( i <= nSearch_nw )
      {
        j = CompPartialMax( vv3EdgeStrs, chkpt, i + nw );
        if( vv3EdgeStrs[i][0] > vv3EdgeStrs[j][0] )
        {
          if( vv3EdgeStrs[i][0] > mfEdgeThreshold ) viEdgeIndex.push_back( i );
          i = i + nw + 1;
          break;
        }
        else
        {
          chkpt = i + nw - 1;
          i = j;
        }
      }
    }
  }
}

const REAL_TYPE COSANGLE=0.98481;
// select only the edeg with normal angle is less than COSANGLE cos( 10 degree=10*pi/180 )
// A.B = |A||B|cos(theta) --> x0*x1 + y0*y1 = cos(theta)

unsigned int EdgeDetector::select_edge_point( DetectedEdge& result, unsigned int no_edge,
    const std::vector<int>& viEdgeIndex, const std::vector<TooN::Vector<3> > & vv3EdgeStrs,
    const TooN::Vector<2>& v2Normal )
{
  REAL_TYPE subpixel;
  const int & nSearch = viEdgeIndex[no_edge];
  result.bNotFound = true;
  if( fabs(v2Normal*vv3EdgeStrs[nSearch].slice(1,2) ) > COSANGLE ) // select the first edge passing the condition
  {
    // subpixel m = a-c/(2*(a-2b-c) where the parabola passing though (-1,a), (0,b), (1,c)
    // A Non-Maxima Suppression Method for Edge Detection with Sub-Pixel Accuracy.
    float den = ( 2.0 * ( vv3EdgeStrs[nSearch - 1][0] - 2.0 * vv3EdgeStrs[nSearch][0] + vv3EdgeStrs[nSearch + 1][0] ) );
    if( nSearch >= 1 && den != 0 )
      subpixel = ( vv3EdgeStrs[nSearch - 1][0] - vv3EdgeStrs[nSearch + 1][0] ) / den;
    else subpixel = 0.0;
    result.rDistance = nSearch - miMaxSearchRange + subpixel;
    result.bNotFound = false;
  }
  ++no_edge;
  if( no_edge < viEdgeIndex.size() ) return no_edge;
  else return 0;
}

int EdgeDetector::CompPartialMax( const std::vector<TooN::Vector<3> > & vv3EdgeStrs, int from, int to )
{
  mvrPMAX[to] = vv3EdgeStrs[to][0];
  int best = to;
  while( to > from )
  {
    to = to - 1;
    if( vv3EdgeStrs[to][0] <= vv3EdgeStrs[best][0] ) mvrPMAX[to] = vv3EdgeStrs[best][0];
    else
    {
      mvrPMAX[to] = vv3EdgeStrs[to][0];
      best = to;
    }
  }
  return best;
}

void EdgeDetector::edgeStrength( TooN::Vector<3> &edgeStr, const CVD::Image<REAL_TYPE> &image, int x, int y )
{
#ifndef _SOBEL_
  REAL_TYPE valx = 0, valy = 0;
  for (int r = -2; r < 3; ++r)
  {
    for (int c = -2; c < 3; ++c)
    {
      REAL_TYPE intensity = (REAL_TYPE)image[y+r][x+c];
      valx += intensity * MaskX[r + 2][c + 2];
      valy += intensity * MaskY[r + 2][c + 2];
    }
  }
#ifndef _SYM_
  valx /= 18;
  valy /= 18;
#else
  valx /= 100;
  valy /= 100;
#endif
  edgeStr[0] = sqrtf (valx * valx + valy * valy);
  REAL_TYPE inv_div = 1.f/edgeStr[0];
  edgeStr[1] = valx*inv_div;
  edgeStr[2] = valy*inv_div;
#else
  REAL_TYPE valx = 0, valy = 0;
  for( int r = -1; r < 2; ++r )
  {
    for( int c = -1; c < 2; ++c )
    {
      REAL_TYPE intensity = (REAL_TYPE) image[y + r][x + c];
      valx += intensity * MaskX[r + 1][c + 1];
      valy += intensity * MaskY[r + 1][c + 1];
    }
  }
  valx /= 4;
  valy /= 4;
  edgeStr[0] = sqrtf( valx * valx + valy * valy );
  REAL_TYPE inv_div = 1.f/edgeStr[0];
  edgeStr[1] = valx*inv_div;
  edgeStr[2] = valy*inv_div;
#endif
}

void EdgeDetector::gen_Normal( TooN::Vector<2>& v2Normal, TooN::Vector<2>& v2Direction )
{
  TooN::normalize( v2Direction );
  v2Normal = TooN::makeVector( -1.f * v2Direction[1], v2Direction[0] );
}

void EdgeDetector::search_for_closest_edge_point( ClosestEdge& result, std::vector<TooN::Vector<3> > & vv3EdgeStrs,
    const CVD::Image<REAL_TYPE>& image, const TooN::Vector<2> & v2ImagePoint, const TooN::Vector<2>& v2Normal )
{
  TooN::Vector<3> edgeStr;
  TooN::Vector<2> v2Xn;
  register int i, j, k, chkpt, nSearch, nSearch_2nw, nSearch_nw, TwoiMaxSearchRange;
  TwoiMaxSearchRange = 2*miMaxSearchRange;
  // apply mask to find edgeStr in both magnitude and direction for a whole range of search path.
  for( nSearch = 0; nSearch <= TwoiMaxSearchRange; ++nSearch )
  {
    v2Xn = v2ImagePoint + ( nSearch - miMaxSearchRange ) * v2Normal;
    if( !image.in_image_with_border( CVD::ir( v2Xn ), 10 ) ) break;
    edgeStrength( edgeStr, image, int( v2Xn[0] ), int( v2Xn[1] ) );
    vv3EdgeStrs[nSearch] = edgeStr;
  }
  // 1D NMS for (2nw+1)-Neighborhood
  // Efficient Non-Maximum Suppression, Alexander Neubeck and Luc Van Gool, ICPR, 2006
  result.bNotFound = true;
  REAL_TYPE subpixel, distance;
  result.rDistance = TwoiMaxSearchRange + 1.f; 
  const int nw = 2; // test nw = 2, must be more than or equal to 2;
  i = nw;
  CompPartialMax( vv3EdgeStrs, 0, i - 1 );
  chkpt = -1;
  nSearch_2nw = nSearch - 2*nw;
  nSearch_nw = nSearch - nw;
  while( i <= nSearch_2nw )
  {
    j = CompPartialMax( vv3EdgeStrs, i, i + nw );
    k = CompPartialMax( vv3EdgeStrs, i + nw + 1, j + nw );
    if( i == j || vv3EdgeStrs[j][0] > vv3EdgeStrs[k][0] )
    {
      if( ( chkpt <= j - nw || vv3EdgeStrs[j][0] >= mvrPMAX[chkpt] )
          && ( j - nw == i || vv3EdgeStrs[j][0] >= mvrPMAX[j - nw] ) )
      {
        if( vv3EdgeStrs[j][0] > mfEdgeThreshold ) // found edge
        {
          if( fabs(v2Normal*vv3EdgeStrs[j].slice(1,2) ) > COSANGLE ) // check angle
          {
            distance = j - miMaxSearchRange;
            if( fabs( result.rDistance ) >= fabs( distance ) )
            {
              result.rDistance = distance;
              result.iEdgeIndex = j;
              result.bNotFound = false;
            }
            else 
              goto exit_label;
          }
        }
      }
      if( i < j ) chkpt = i + nw + 1;
      i = j + nw + 1;
    }
    else
    {
      i = k;
      chkpt = j + nw + 1;
      while( i <= nSearch_nw )
      {
        j = CompPartialMax( vv3EdgeStrs, chkpt, i + nw );
        if( vv3EdgeStrs[i][0] > vv3EdgeStrs[j][0] )
        {
          if( vv3EdgeStrs[i][0] > mfEdgeThreshold ) // found edge
          {
            if( fabs(v2Normal*vv3EdgeStrs[i].slice(1,2)) > COSANGLE ) // check angle
            {
              distance = i - miMaxSearchRange;
              if( fabs( result.rDistance ) >= fabs( distance) )
              {
                result.rDistance = distance;
                result.iEdgeIndex = i;
                result.bNotFound = false;
              }
              else
                goto exit_label;
            }
          }
          i = i + nw + 1;
          break;
        }
        else
        {
          chkpt = i + nw - 1;
          i = j;
        }
      }
    }
  }
  exit_label:
  if( !result.bNotFound ) 
  {
// subpixel m = a-c/(2*(a-2b-c) where the parabola passing though (-1,a), (0,b), (1,c)
// A Non-Maxima Suppression Method for Edge Detection with Sub-Pixel Accuracy.
    REAL_TYPE den = ( 2.0 * ( vv3EdgeStrs[result.iEdgeIndex - 1][0] - 2.0 * vv3EdgeStrs[result.iEdgeIndex][0] + vv3EdgeStrs[result.iEdgeIndex + 1][0] ) );
    if( result.iEdgeIndex >= 1 && den != 0 )
      subpixel = ( vv3EdgeStrs[result.iEdgeIndex - 1][0] - vv3EdgeStrs[result.iEdgeIndex + 1][0] ) / den;
    else subpixel = 0.0;
    result.rDistance += subpixel;
  }
}

void EdgeDetector::search_for_closest_edge_point_with_pole( ClosestEdge& result, std::vector<TooN::Vector<3> > & vv3EdgeStrs, 
    const EdgeFeature & edgeFeature, const CVD::Image<REAL_TYPE>& image, const TooN::Vector<2> & v2ImagePoint, const TooN::Vector<2>& v2Normal )
{
  TooN::Vector<3> edgeStr;
  TooN::Vector<2> v2Xn;
  register int i, j, k, chkpt, nSearch, nSearch_2nw, nSearch_nw, TwoiMaxSearchRange;
  // apply mask to find edgeStr in both magnitude and direction for a whole range of search path.
  TwoiMaxSearchRange = 2*miMaxSearchRange;
  for( nSearch = 0; nSearch <= TwoiMaxSearchRange; ++nSearch )
  {
    v2Xn = v2ImagePoint + ( nSearch - miMaxSearchRange ) * v2Normal;
    if( !image.in_image_with_border( CVD::ir( v2Xn ), 10 ) ) break;
    edgeStrength( edgeStr, image, int( v2Xn[0] ), int( v2Xn[1] ) );
    vv3EdgeStrs[nSearch] = edgeStr;
  }
  // 1D NMS for (2nw+1)-Neighborhood
  // Efficient Non-Maximum Suppression, Alexander Neubeck and Luc Van Gool, ICPR, 2006
  result.bNotFound = true;
  REAL_TYPE subpixel, distance, test;
  result.rDistance = 2 * miMaxSearchRange + 1.f; 
  const int nw = 2; // test nw = 2, must be more than or equal to 2;
  i = nw;
  CompPartialMax( vv3EdgeStrs, 0, i - 1 );
  chkpt = -1;
  nSearch_2nw = nSearch - 2*nw;
  nSearch_nw = nSearch - nw;
  while( i <= nSearch_2nw )
  {
    j = CompPartialMax( vv3EdgeStrs, i, i + nw );
    k = CompPartialMax( vv3EdgeStrs, i + nw + 1, j + nw );
    if( i == j || vv3EdgeStrs[j][0] > vv3EdgeStrs[k][0] )
    {
      if( ( chkpt <= j - nw || vv3EdgeStrs[j][0] >= mvrPMAX[chkpt] )
          && ( j - nw == i || vv3EdgeStrs[j][0] >= mvrPMAX[j - nw] ) )
      {
        if( vv3EdgeStrs[j][0] > mfEdgeThreshold ) // found edge
        {
          test = v2Normal*vv3EdgeStrs[j].slice(1,2);
          if( ( IS_SET_BIT( edgeFeature.ucBits, EDGE_POLE ) && test > 0 ) || ( !IS_SET_BIT( edgeFeature.ucBits, EDGE_POLE ) && test <= 0 ) )
          {
            if( fabs(test) > COSANGLE ) // check angle
            {
              distance = j - miMaxSearchRange;
              if( fabs( result.rDistance ) >= fabs( distance ) )
              {
                result.rDistance = distance;
                result.iEdgeIndex = j;
                result.bNotFound = false;
              }
              else
                goto exit_label;
            }
          }
        }
      }
      if( i < j ) chkpt = i + nw + 1;
      i = j + nw + 1;
    }
    else
    {
      i = k;
      chkpt = j + nw + 1;
      while( i <= nSearch_nw )
      {
        j = CompPartialMax( vv3EdgeStrs, chkpt, i + nw );
        if( vv3EdgeStrs[i][0] > vv3EdgeStrs[j][0] )
        {
          if( vv3EdgeStrs[i][0] > mfEdgeThreshold ) // found edge
          {
            test = v2Normal*vv3EdgeStrs[i].slice(1,2);
            if( ( IS_SET_BIT( edgeFeature.ucBits, EDGE_POLE ) && test > 0 ) || ( !IS_SET_BIT( edgeFeature.ucBits, EDGE_POLE ) && test <= 0 ) )
            {
              if( fabs(test) > COSANGLE ) // check angle
              {
                distance = i - miMaxSearchRange;
                if( fabs( result.rDistance ) >= fabs( distance) )
                {
                  result.rDistance = distance;
                  result.iEdgeIndex = i;
                  result.bNotFound = false;
                }
                else
                  goto exit_label;
              }
            }
          }
          i = i + nw + 1;
          break;
        }
        else
        {
          chkpt = i + nw - 1;
          i = j;
        }
      }
    }
  }
  exit_label:
  if( !result.bNotFound )
  {
// subpixel m = a-c/(2*(a-2b-c) where the parabola passing though (-1,a), (0,b), (1,c)
// A Non-Maxima Suppression Method for Edge Detection with Sub-Pixel Accuracy.
    REAL_TYPE den = ( 2.0 * ( vv3EdgeStrs[result.iEdgeIndex - 1][0] - 2.0 * vv3EdgeStrs[result.iEdgeIndex][0] + vv3EdgeStrs[result.iEdgeIndex + 1][0] ) );
    if( result.iEdgeIndex >= 1 && den != 0 )
      subpixel = ( vv3EdgeStrs[result.iEdgeIndex - 1][0] - vv3EdgeStrs[result.iEdgeIndex + 1][0] ) / den;
    else subpixel = 0.0;
    result.rDistance += subpixel;
  }
}

void EdgeDetector::search_for_edge_point( EdgeData& ed, std::vector<TooN::Vector<3> > & vv3EdgeStrs, const CVD::Image<REAL_TYPE> & image,
        const TooN::Vector<2> & v2ImagePoint, const TooN::Vector<2>& v2Normal )
{
  TooN::Vector<3> edgeStr;
  TooN::Vector<2> v2Xn;
  ed.v2ExpectedImagePoint = v2ImagePoint;
  ed.v2Normal = v2Normal;
  register int i, j, k, chkpt, nSearch, nSearch_2nw, nSearch_nw, TwoiMaxSearchRange;
  TwoiMaxSearchRange = 2*miMaxSearchRange;
  // apply mask to find edgeStr in both magnitude and direction for a whole range of search path.
  for( nSearch = 0; nSearch <= TwoiMaxSearchRange; ++nSearch )
  {
    v2Xn = v2ImagePoint + ( nSearch - miMaxSearchRange ) * v2Normal;
    if( !image.in_image_with_border( CVD::ir( v2Xn ), 10 ) ) break; // reach boundary.
    edgeStrength( edgeStr, image, int( v2Xn[0] ), int( v2Xn[1] ) );
    vv3EdgeStrs[nSearch] = edgeStr;
  }
  // 1D NMS for (2nw+1)-Neighborhood
  // Efficient Non-Maximum Suppression, Alexander Neubeck and Luc Van Gool, ICPR, 2006
  REAL_TYPE subpixel, distance;
  const int nw = 2; // test nw = 2, must be more than or equal to 2;
  i = nw;
  CompPartialMax( vv3EdgeStrs, 0, i - 1 );
  chkpt = -1;
  nSearch_2nw = nSearch - 2*nw;
  nSearch_nw = nSearch - nw;
  while( i <= nSearch_2nw )
  {
    j = CompPartialMax( vv3EdgeStrs, i, i + nw );
    k = CompPartialMax( vv3EdgeStrs, i + nw + 1, j + nw );
    if( i == j || vv3EdgeStrs[j][0] > vv3EdgeStrs[k][0] )
    {
      if( ( chkpt <= j - nw || vv3EdgeStrs[j][0] >= mvrPMAX[chkpt] )
          && ( j - nw == i || vv3EdgeStrs[j][0] >= mvrPMAX[j - nw] ) )
      {
        if( vv3EdgeStrs[j][0] > mfEdgeThreshold ) // found edge
        {
          if( fabs(v2Normal*vv3EdgeStrs[j].slice(1,2) ) > COSANGLE ) // check angle
          {
            distance = j - miMaxSearchRange;
// subpixel m = a-c/(2*(a-2b-c) where the parabola passing though (-1,a), (0,b), (1,c)
// A Non-Maxima Suppression Method for Edge Detection with Sub-Pixel Accuracy.
            REAL_TYPE den = ( 2.0 * ( vv3EdgeStrs[j - 1][0] - 2.0 * vv3EdgeStrs[j][0] + vv3EdgeStrs[j + 1][0] ) );
            if( j >= 1 && den != 0 )
              subpixel = ( vv3EdgeStrs[j - 1][0] - vv3EdgeStrs[j + 1][0] ) / den;
            else subpixel = 0.0;
            distance += subpixel;
            ed.vrDistance.push_back( distance );
            ed.vrDistanceSquared.push_back( distance*distance );
            ++ed.no_candidate;
          }
        }
      }
      if( i < j ) chkpt = i + nw + 1;
      i = j + nw + 1;
    }
    else
    {
      i = k;
      chkpt = j + nw + 1;
      while( i <= nSearch_nw )
      {
        j = CompPartialMax( vv3EdgeStrs, chkpt, i + nw );
        if( vv3EdgeStrs[i][0] > vv3EdgeStrs[j][0] )
        {
          if( vv3EdgeStrs[i][0] > mfEdgeThreshold ) // found edge
          {
            if( fabs(v2Normal*vv3EdgeStrs[i].slice(1,2)) > COSANGLE ) // check angle
            {
              distance = i - miMaxSearchRange;
// subpixel m = a-c/(2*(a-2b-c) where the parabola passing though (-1,a), (0,b), (1,c)
// A Non-Maxima Suppression Method for Edge Detection with Sub-Pixel Accuracy.
              REAL_TYPE den = ( 2.0 * ( vv3EdgeStrs[i - 1][0] - 2.0 * vv3EdgeStrs[i][0] + vv3EdgeStrs[i + 1][0] ) );
              if( i >= 1 && den != 0 )
                subpixel = ( vv3EdgeStrs[i - 1][0] - vv3EdgeStrs[i + 1][0] ) / den;
              else subpixel = 0.0;
              distance += subpixel;
              ed.vrDistance.push_back( distance );
              ed.vrDistanceSquared.push_back( distance*distance );
              ++ed.no_candidate;
            }
          }
          i = i + nw + 1; // in the paper, it is i = i + nw - 1, it is wrong.
          break;
        }
        else
        {
          chkpt = i + nw - 1;
          i = j;
        }
      }
    }
  }
}

