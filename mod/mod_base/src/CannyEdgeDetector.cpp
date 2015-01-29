#include "CannyEdgeDetector.h"

void CannyEdgeDetector::compute( CVD::Image<CVD::byte> & imEdge, Edges & edges, LowLevelImageData<CVD::byte> & llimSource, unsigned int uiShortEdge, REAL_TYPE maxThreshold, REAL_TYPE minThreshold )
{
  assert( imEdge.size() == llimSource.mImage.size() );
// smooth image with simga = 1.6
// and compute gradient
// and compute magnitude
  llimSource.ComputeGradMagnitudeImage( 1.6 );
  llimSource.NormalizeGradMagnitudeImage();
// estimate threshold
#if 0
  TooN::Vector<2> threshold;
  estimateThreshold( threshold, llimSource.mGradMagnitudeImage, maxThreshold, minThreshold );
// non maximum suppression
  std::vector< EdgeProperties > vPossibleEdges;
  non_max_suppression( imEdge, vPossibleEdges, llimSource, threshold[0] );
  hysteresis( imEdge, edges, vPossibleEdges, uiShortEdge, threshold[1] );
#else
  std::vector< EdgeProperties > vPossibleEdges;
  non_max_suppression( imEdge, vPossibleEdges, llimSource, 0.05 );
  hysteresis( imEdge, edges, vPossibleEdges, uiShortEdge, 0.2 );
#endif
}

void CannyEdgeDetector::non_max_suppression( CVD::Image<CVD::byte> & imEdge, std::vector<EdgeProperties> & vPossibleEdges, const LowLevelImageData<CVD::byte> & llimSource, REAL_TYPE lowThreshold )
{
  vPossibleEdges.clear();
  imEdge.fill(0);
  int w = imEdge.size().x;
  CVD::BasicImage<REAL_TYPE>::const_iterator magIter = llimSource.mGradMagnitudeImage.begin() + w;
  CVD::BasicImage<REAL_TYPE>::const_iterator magEndIter = llimSource.mGradMagnitudeImage.begin() + (imEdge.size().y-1)*w;
  CVD::BasicImage<REAL_TYPE[2]>::const_iterator gradIter = llimSource.mGradientImage.begin() + w;
  CVD::ImageRef irPose = CVD::ImageRef(0,1);
  CVD::ImageRef irMax = imEdge.size();
  while( magIter != magEndIter )
  {
    REAL_TYPE rMag = *magIter;
    if( rMag > lowThreshold )
    {
      REAL_TYPE dx = (*gradIter)[0];
      REAL_TYPE dy = (*gradIter)[1];
      REAL_TYPE m1, m2;
      if( dx >= 0 )
      {
        if( dy >= 0 ) // 0 - 90
        {
          if( dx >= dy ) // 0 - 45
          {
            m1 = -( rMag - *(magIter-1) )*dx + ( *(magIter-w-1) - *(magIter-1) )*dy;
            m2 = -( rMag - *(magIter+1) )*dx + ( *(magIter+w+1) - *(magIter+1) )*dy;
          }
          else // 45 - 90
          {
            m1 = -( *(magIter-w) - *(magIter-w-1) )*dx + ( *(magIter-w) - rMag )*dy;
            m2 = -( *(magIter+w) - *(magIter+w+1) )*dx + ( *(magIter+w) - rMag )*dy;
          }
        }
        else // 270 - 360
        {
          if( dx >= -dy ) // 315 - 360
          {
            m1 = -( rMag - *(magIter-1) )*dx + ( *(magIter-1) - *(magIter+w-1) )*dy;
            m2 = -( rMag - *(magIter+1) )*dx + ( *(magIter+1) - *(magIter-w+1) )*dy;
          }
          else // 270 - 315
          {
            m1 = -( *(magIter+w) - *(magIter+w-1) )*dx + ( rMag - *(magIter+w) )*dy;
            m2 = -( *(magIter-w) - *(magIter-w+1) )*dx + ( rMag - *(magIter-w) )*dy;
          }
        }
      }
      else
      {
        if( dy >= 0 ) // 90 - 180
        {
          if( -dx >= dy ) // 135 - 180
          {
            m1 = -( *(magIter+1) - rMag )*dx + ( *(magIter-w+1) - *(magIter+1) )*dy;
            m2 = -( *(magIter-1) - rMag )*dx + ( *(magIter+w-1) - *(magIter-1) )*dy;
          }
          else // 90 - 135
          {
            m1 = -( *(magIter-w+1) - *(magIter-w) )*dx + ( *(magIter-w) - rMag )*dy;
            m2 = -( *(magIter+w-1) - *(magIter+w) )*dx + ( *(magIter+w) - rMag )*dy;
          }
        }
        else // 180 - 270
        {
          if( dx <= dy ) // 180 - 225
          {
            m1 = -( *(magIter+1) - rMag )*dx + ( *(magIter+1) - *(magIter+w+1) )*dy;
            m2 = -( *(magIter-1) - rMag )*dx + ( *(magIter-1) - *(magIter-w-1) )*dy;
          }
          else // 225 - 270
          {
            m1 = -( *(magIter+w+1) - *(magIter+w) )*dx + ( rMag - *(magIter+w) )*dy;
            m2 = -( *(magIter-w-1) - *(magIter-w) )*dx + ( rMag - *(magIter-w) )*dy;
          }
        }
      }
      if( m1 < 0 && m2 < 0 )
      {
        imEdge[irPose] = 32;
        EdgeProperties edgeProperties;
        edgeProperties.irPossibleEdge = irPose;
        edgeProperties.pGradMagnitude = magIter;
        vPossibleEdges.push_back( edgeProperties );
      }
    } // end if rMag more than threshold 
    irPose.next(irMax);
    ++magIter;
    ++gradIter;
  }
  CVD::zeroBorders( imEdge );
}

void CannyEdgeDetector::estimateThreshold( TooN::Vector<2> & threshold, const CVD::Image<REAL_TYPE> & magImage, REAL_TYPE maxThreshold, REAL_TYPE minThreshold )
{
  CVD::Image<CVD::byte> mag = CVD::convert_image(magImage);
  int noPixels = magImage.size().x * magImage.size().y;
  int histogram[256];
  memset(histogram, 0, 256*sizeof(int));
  CVD::BasicImage<CVD::byte>::iterator beginIter = mag.begin();
  CVD::BasicImage<CVD::byte>::iterator endIter = mag.end();
  while( beginIter != endIter )
  {
    histogram[*beginIter] += 1;
    ++beginIter;
  }
  REAL_TYPE sumCount=0;
  int highThreshold=7;
  for( int i = 0; i < highThreshold; ++i )
    sumCount += histogram[i];
  REAL_TYPE noPixelsMaxThreshold = maxThreshold*(REAL_TYPE)(noPixels - sumCount);
  sumCount = 0;
  while( sumCount < noPixelsMaxThreshold )
  {
    sumCount += histogram[highThreshold++];
  }
  threshold[1] = ((REAL_TYPE)highThreshold)/255.f;
  threshold[0] = minThreshold*threshold[1];
}

void CannyEdgeDetector::hysteresis( CVD::Image<CVD::byte> & imEdge, Edges & edges, std::vector<EdgeProperties> & vPossibleEdges, unsigned int uiShortEdge, REAL_TYPE highThreshold )
{
  std::vector<EdgeProperties>::iterator beginIter = vPossibleEdges.begin();
  std::vector<EdgeProperties>::iterator endIter = vPossibleEdges.end();
  std::vector<CVD::ImageRef> virEdges;
  while( beginIter != endIter )
  {
    if( *((*beginIter).pGradMagnitude) > highThreshold && imEdge[(*beginIter).irPossibleEdge] == 32 )
    {
      imEdge[(*beginIter).irPossibleEdge] = 255;
      virEdges.push_back((*beginIter).irPossibleEdge);
      follow_edges( imEdge, virEdges, (*beginIter).irPossibleEdge );
    }
    ++beginIter;
  }
//  for( int i = 0; i < virEdges.size(); ++i )
//    std::cout << virEdges[i].x << " " << virEdges[i].y << std::endl;
  thinning( imEdge, virEdges );
//  std::cout <<"Thinning " << virEdges.size() << std::endl;
//  for( int i = 0; i < virEdges.size(); ++i )
//    std::cout << virEdges[i].x << " " << virEdges[i].y << std::endl;
  link_edges( imEdge, edges.virEdges, edges.viEdgeIdxes, virEdges, uiShortEdge );
//  std::cout <<"Linking " << edges.virEdges.size() << std::endl;
//  for( int i = 0; i < edges.virEdges.size(); ++i )
//    std::cout << edges.virEdges[i].x << " " << edges.virEdges[i].y << std::endl;
//  for( int i = 0; i < edges.viEdgeIdxes.size(); ++i )
//    std::cout << edges.virEdges[ edges.viEdgeIdxes[i] ].x << " " <<  edges.virEdges[ edges.viEdgeIdxes[i] ].y << std::endl;
/*
  edges.virEdges.reserve( virEdges.size() );
  std::vector<CVD::ImageRef>::iterator irBeginIter = virEdges.begin();
  std::vector<CVD::ImageRef>::iterator irEndIter = virEdges.end();
  std::vector<CVD::ImageRef>::iterator irTemp = irBeginIter;
  while( irBeginIter != irEndIter )
  {
    int noPoints = 0;
    std::vector<CVD::ImageRef>::iterator irTempPrevious = irTemp++;
    while( irTemp != irEndIter )
    {
      int x = abs( (*irTemp).x - (*irTempPrevious).x );
      int y = abs( (*irTemp).y - (*irTempPrevious).y );
      if( x + y < 3 )
      {
        irTempPrevious = irTemp++;
        ++noPoints;
      }
      else
        break;
    }
    if( noPoints > uiShortEdge )
    {
      edges.viEdgeIdxes.push_back(edges.virEdges.size());
      edges.virEdges.insert(edges.virEdges.end(), irBeginIter, irTemp);
    }
    irBeginIter = irTemp;
  } */
  edges.viEdgeIdxes.push_back(edges.virEdges.size());
/*// No need
  beginIter = vPossibleEdges.begin();
  // clear non-edge
  while( beginIter != endIter )
  {
    if( imEdge[(*beginIter).irPossibleEdge] == 32 )
      imEdge[(*beginIter).irPossibleEdge] = 0;
    ++beginIter;
  }
*/
}

void CannyEdgeDetector::follow_edges( CVD::Image<CVD::byte> & imEdge, std::vector<CVD::ImageRef> & virEdges, CVD::ImageRef & irPossibleEdge )
{
  CVD::ImageRef irNext;
  irNext = irPossibleEdge+CVD::ImageRef(1,0);
  if( imEdge[irNext] == 32 )
  {
    imEdge[irNext] = 255;
    virEdges.push_back(irNext);
    follow_edges( imEdge, virEdges, irNext );
  }
  irNext = irPossibleEdge+CVD::ImageRef(1,1);
  if( imEdge[irNext] == 32 )
  {
    imEdge[irNext] = 255;
    virEdges.push_back(irNext);
    follow_edges( imEdge, virEdges, irNext );
  }
  irNext = irPossibleEdge+CVD::ImageRef(0,1);
  if( imEdge[irNext] == 32 )
  {
    imEdge[irNext] = 255;
    virEdges.push_back(irNext);
    follow_edges( imEdge, virEdges, irNext );
  }
  irNext = irPossibleEdge+CVD::ImageRef(-1,1);
  if( imEdge[irNext] == 32 )
  {
    imEdge[irNext] = 255;
    virEdges.push_back(irNext);
    follow_edges( imEdge, virEdges, irNext );
  }
  irNext = irPossibleEdge+CVD::ImageRef(-1,0);
  if( imEdge[irNext] == 32 )
  {
    imEdge[irNext] = 255;
    virEdges.push_back(irNext);
    follow_edges( imEdge, virEdges, irNext );
  }
  irNext = irPossibleEdge+CVD::ImageRef(-1,-1);
  if( imEdge[irNext] == 32 )
  {
    imEdge[irNext] = 255;
    virEdges.push_back(irNext);
    follow_edges( imEdge, virEdges, irNext );
  }
  irNext = irPossibleEdge+CVD::ImageRef(0,-1);
  if( imEdge[irNext] == 32 )
  {
    imEdge[irNext] = 255;
    virEdges.push_back(irNext);
    follow_edges( imEdge, virEdges, irNext );
  }
  irNext = irPossibleEdge+CVD::ImageRef(1,-1);
  if( imEdge[irNext] == 32 )
  {
    imEdge[irNext] = 255;
    virEdges.push_back(irNext);
    follow_edges( imEdge, virEdges, irNext );
  }
}

void CannyEdgeDetector::thinning( CVD::Image<CVD::byte> & imEdge, std::vector<CVD::ImageRef> & virEdges )
{
  std::vector<CVD::ImageRef> virMasks;
  std::vector<CVD::ImageRef> virEdgesProcess;
  unsigned char eightNeighbours[8];
  bool bCheck = true;
  while( bCheck )
  {
    bCheck = false;
    virMasks.clear(); virMasks.reserve( virEdges.size() );
    virEdgesProcess.clear(); virEdgesProcess.reserve( virEdges.size() );
    std::vector<CVD::ImageRef>::iterator irBeginIter = virEdges.begin();
    std::vector<CVD::ImageRef>::iterator irEndIter = virEdges.end();

// first pass  
    while( irBeginIter != irEndIter )
    {
      int np = 0, sr = 0;
      for( int i = 0; i < 7; ++i )
      {
        get_eightNeighbours( eightNeighbours, imEdge, *irBeginIter );
        np += eightNeighbours[i];
        if( eightNeighbours[i] == 0 && eightNeighbours[i+1] == 1 )
          ++sr;
      }
      np += eightNeighbours[7];
      if( eightNeighbours[7] == 0 && eightNeighbours[0] == 1 )
        ++sr;
      int p024=eightNeighbours[0]*eightNeighbours[2]*eightNeighbours[4];
      int p246=eightNeighbours[2]*eightNeighbours[4]*eightNeighbours[6];
      if( sr == 1 && np >= 2 && np <= 6 && p024 == 0 && p246 == 0 && eightNeighbours[5] != 0 )
      {
        virMasks.push_back(*irBeginIter);
        bCheck = true;
      }
      else
        virEdgesProcess.push_back(*irBeginIter);
      ++irBeginIter;
    }
    irBeginIter = virMasks.begin(); irEndIter = virMasks.end();
    while( irBeginIter != irEndIter )
    {
      imEdge[*irBeginIter] = 64;
      ++irBeginIter;
    }
  // second pass  
    irBeginIter = virEdgesProcess.begin(); irEndIter = virEdgesProcess.end();
    virMasks.clear(); virMasks.reserve( virEdgesProcess.size() );
    virEdges.clear(); virEdges.reserve( virEdgesProcess.size() );
    while( irBeginIter != irEndIter )
    {
      get_eightNeighbours( eightNeighbours, imEdge, *irBeginIter );
      int np = 0, sr = 0;
      for( int i = 0; i < 7; ++i )
      {
        np += eightNeighbours[i];
        if( eightNeighbours[i] == 0 && eightNeighbours[i+1] == 1 )
          ++sr;
      }
      np += eightNeighbours[7];
      if( eightNeighbours[7] == 0 && eightNeighbours[0] == 1 )
        ++sr;
      int p026=eightNeighbours[0]*eightNeighbours[2]*eightNeighbours[6];
      int p046=eightNeighbours[0]*eightNeighbours[4]*eightNeighbours[6];
      if( sr == 1 && np >= 2 && np <= 6 && p026 == 0 && p046 == 0 && eightNeighbours[1] != 0 )
      {
        virMasks.push_back(*irBeginIter);
        bCheck = true;
      }
      else
        virEdges.push_back(*irBeginIter);
      ++irBeginIter;
    }
    irBeginIter = virMasks.begin(); irEndIter = virMasks.end();
    while( irBeginIter != irEndIter )
    {
      imEdge[*irBeginIter] = 64;
      ++irBeginIter;
    }
  }
}

void CannyEdgeDetector::link_edges( CVD::Image<CVD::byte> & imEdge, std::vector<CVD::ImageRef> & virLinkedEdges, std::vector<int> & viLinkedEdgesIndexes, std::vector<CVD::ImageRef> & virEdges, unsigned int uiShortEdge )
{
  virLinkedEdges.clear();
  virLinkedEdges.reserve( virEdges.size() );
  viLinkedEdgesIndexes.clear();
  std::vector<CVD::ImageRef> virLeft;
  std::vector<CVD::ImageRef> virRight;
  std::vector<CVD::ImageRef> virEightNeighbours(8);
  virEightNeighbours[0] = CVD::ImageRef(1,0); virEightNeighbours[1] = CVD::ImageRef(1,1);
  virEightNeighbours[2] = CVD::ImageRef(0,1); virEightNeighbours[3] = CVD::ImageRef(-1,1);
  virEightNeighbours[4] = CVD::ImageRef(-1,0); virEightNeighbours[5] = CVD::ImageRef(-1,-1);
  virEightNeighbours[6] = CVD::ImageRef(0,-1); virEightNeighbours[7] = CVD::ImageRef(1,-1);
  int no_possible_edges, prev_no_possible_edges;
  prev_no_possible_edges = virEdges.size();
  do
  {
    std::vector<CVD::ImageRef>::iterator irBeginIter = virEdges.begin();
    std::vector<CVD::ImageRef>::iterator irEndIter = virEdges.end();
    no_possible_edges = virEdges.size();
    while( irBeginIter != irEndIter )
    {
      virLeft.clear(); virLeft.reserve( virEdges.size() ); virRight.clear(); virRight.reserve( virEdges.size() );
      unsigned char eightNeighbours[8];
      CVD::ImageRef irCurrent = *irBeginIter;
      ++irBeginIter;
      if( imEdge[irCurrent] == 255 )
      {
        int no_cross, no_edges;
        get_eightNeighbours( eightNeighbours, imEdge, irCurrent );
        no_crossing0to1( no_cross, no_edges, eightNeighbours );
        if( no_cross == 1 || no_cross >= 3 ) // this is an end point or junction
        {
          if( no_cross >= 3 )
          {
            imEdge[irCurrent] = 128;
            virLeft.push_back( irCurrent );
            unsigned char idx = 0; while( eightNeighbours[idx] == 0 ) ++idx;
            if( ( (1 << idx) & 0b10101010 ) == 0 )
              irCurrent += virEightNeighbours[idx];
            else if( eightNeighbours[idx+1] ) // if ( 1 << idx & 0b10101010 != 0 )
              irCurrent += virEightNeighbours[idx+1];
            else
              irCurrent += virEightNeighbours[idx];
            get_eightNeighbours( eightNeighbours, imEdge, irCurrent );
            no_crossing0to1( no_cross, no_edges, eightNeighbours );
            if( no_cross != 1 )
              continue;
          }
          do
          {
            imEdge[irCurrent] = 128;
            virLeft.push_back( irCurrent );
            unsigned char idx = 0; while( eightNeighbours[idx] == 0 ) ++idx;
            if( ( (1 << idx) & 0b10101010 ) == 0 )
              irCurrent += virEightNeighbours[idx];
            else if( eightNeighbours[idx+1] )
              irCurrent += virEightNeighbours[idx+1];
            else
              irCurrent += virEightNeighbours[idx];
            get_eightNeighbours( eightNeighbours, imEdge, irCurrent );
            no_crossing0to1( no_cross, no_edges, eightNeighbours );
          } while( no_cross == 1 );
          if( no_edges == 0 ) // end point
            virLeft.push_back( irCurrent );
          if( virLeft.size() > uiShortEdge )
          {
            viLinkedEdgesIndexes.push_back(virLinkedEdges.size());
            virLinkedEdges.insert(virLinkedEdges.end(), virLeft.begin(), virLeft.end());
          }
        }
        else if( no_cross == 0 ) // isolate point
          imEdge[irCurrent] = 128;
#if 1          
        else if( no_cross == 2 ) 
        {
          CVD::ImageRef irCurrentCopy = irCurrent;
          imEdge[irCurrent] = 128;
          virLeft.push_back( irCurrent );
          unsigned char idx_fst = 0; while( eightNeighbours[idx_fst] == 0 ) ++idx_fst;
          unsigned char idx_snd = idx_fst+1; while( eightNeighbours[idx_snd] == 1 ) ++idx_snd;
          while( eightNeighbours[idx_snd] == 0 ) ++idx_snd;
          // left side
          if( ( (1 << idx_fst) & 0b10101010 ) == 0 )
            irCurrent += virEightNeighbours[idx_fst];
          else if( eightNeighbours[idx_fst+1] )
            irCurrent += virEightNeighbours[idx_fst+1];
          else
            irCurrent += virEightNeighbours[idx_fst];
          get_eightNeighbours( eightNeighbours, imEdge, irCurrent );
          no_crossing0to1( no_cross, no_edges, eightNeighbours );
          if( no_cross == 1 )
          {
            do
            {
              imEdge[irCurrent] = 128;
              virLeft.push_back( irCurrent );
              unsigned char idx=0; while( eightNeighbours[idx] == 0 ) ++idx;
              if( ( (1 << idx) & 0b10101010 ) == 0 )
                irCurrent += virEightNeighbours[idx];
              else if( eightNeighbours[idx+1] )
                irCurrent += virEightNeighbours[idx+1];
              else
                irCurrent += virEightNeighbours[idx];
              get_eightNeighbours( eightNeighbours, imEdge, irCurrent );
              no_crossing0to1( no_cross, no_edges, eightNeighbours );
            } while( no_cross == 1 );
            if( no_edges == 0 )
              virLeft.push_back( irCurrent );
          }
          // right side
          irCurrent = irCurrentCopy;
          if( ( (1 << idx_snd) & 0b10101010 ) == 0 )
            irCurrent += virEightNeighbours[idx_snd];
          else if( eightNeighbours[idx_snd+1] )
            irCurrent += virEightNeighbours[idx_snd+1];
          else
            irCurrent += virEightNeighbours[idx_snd];
          get_eightNeighbours( eightNeighbours, imEdge, irCurrent );
          no_crossing0to1( no_cross, no_edges, eightNeighbours );
          if( no_cross == 1 )
          {
            do
            {
              imEdge[irCurrent] = 128;
              virRight.push_back( irCurrent );
              unsigned char idx=0; while( eightNeighbours[idx] == 0 ) ++idx;
              if( ( (1 << idx) & 0b10101010 ) == 0 )
                irCurrent += virEightNeighbours[idx];
              else if( eightNeighbours[idx+1] )
                irCurrent += virEightNeighbours[idx+1];
              else
                irCurrent += virEightNeighbours[idx];
              get_eightNeighbours( eightNeighbours, imEdge, irCurrent );
              no_crossing0to1( no_cross, no_edges, eightNeighbours );
            } while( no_cross == 1 );
            if( no_edges == 0 )
              virRight.push_back( irCurrent );
          }
          // link left and right side
          if( virRight.size() + virLeft.size() > uiShortEdge )
          {
            viLinkedEdgesIndexes.push_back(virLinkedEdges.size());
            virLinkedEdges.insert( virLinkedEdges.end(), virRight.rbegin(), virRight.rend() );
            virLinkedEdges.insert( virLinkedEdges.end(), virLeft.begin(), virLeft.end() );
          }
        }
#endif        
      }
      else
        --no_possible_edges;
    }
    if( no_possible_edges == prev_no_possible_edges )
      break;
    else
      prev_no_possible_edges = no_possible_edges;
  }while( no_possible_edges );
}

