#include "EdgeExtraction.h"
#include <cvd/vector_image_ref.h>
#include <cvd/timer.h>
#include <fcntl.h>
#include <algorithm>
#include "GetTransformDist.h"
#include "BitMacros.h"

#include <boost/lexical_cast.hpp>
#include <cvd/image_io.h>

#define PRL_THRESHOLD 0.0524
//#define LSD
image_double lsd_image2;
int data_size2;
bool isLSD2 = true;
double start_time2 = 0;

const int EDGELET_LENGTH = 10;

int img_h2, img_w2;

EdgeExtraction::EdgeExtraction( CVD::ImageRef irSize, bool lsdV) : mvEdgeOrientationImages(11), mvEdgeOrientationDistTransformImages(11), mvEdgeOrientations(11), mvEdgeOrientationDistTransforms(11), mEdgeImage( irSize ), mirImageSize( irSize), mLLImagePyramid(irSize, 2), mGrayImage(irSize), mGrayImageF(irSize)
{

  img_h2 = irSize.y;
  img_w2 = irSize.x;
  isLSD2 = lsdV;

  for( int i = 0; i < 11; ++i )
  {
    mvEdgeOrientationImages[i].resize( irSize );
    mvEdgeOrientations[i].resize(irSize);
    mvEdgeOrientations[i].reset();
    mvEdgeOrientationDistTransformImages[i].resize( irSize );
    mvEdgeOrientationDistTransforms[i].resize(irSize);
    mvEdgeOrientationDistTransforms[i].reset();
  }
 lsd_image2 = new_image_double( irSize.x, irSize.y );
 data_size2 = lsd_image2->xsize * lsd_image2->ysize;
}

EdgeExtraction::~EdgeExtraction()
{
   free_image_double(lsd_image2);
}

void EdgeExtraction::preCalculateEdgeLnks( LowLevelImageData<CVD::byte> & llImage)
{
  std::srand( time(NULL) );
  start_time2 = CVD::timer.get_time();
  std::vector<Line> vLines;
  double line_time;
  //std::cout << " rewriting master_file_id to " << master_file_id << std::endl;
if (!isLSD2)
{
  //Edges edges;
  mCannyEdgeDetector.compute( mEdgeImage, edges, llImage, 10, 0.5, 0.2 );
//  for( int i = 0; i < edges.virEdges.size(); ++i )
//    std::cout <<  CVD::vec(edges.virEdges[i]) << std::endl;
  get_lines( vLines, edges );
  line_time = CVD::timer.get_time();
  mEdgeImage.fill(0);
  char filename[100];
// SAVING MASTER FILE
/*  sprintf(filename,"%s/master_%d.xvy", folderName, master_file_id);
  std::ofstream out_xvy(filename);
  for( unsigned int i = 0; i < edges.virEdges.size(); ++i )
  {
    mEdgeImage[edges.virEdges[i]] = 255;
    out_xvy << CVD::vec(edges.virEdges[i]) << std::endl;
  }
  out_xvy.close();
  sprintf(filename,"%s/master_%d.png",folderName, master_file_id);
  //++master_file_id;
  CVD::img_save(mEdgeImage, filename);
*/
}
else
{
  int no_data = 0;
  double * img_ptr = lsd_image2->data;
  CVD::byte * int_img_ptr = llImage.mImage.data();
  while( no_data != data_size2 )
  {
    *img_ptr = static_cast<double>(*int_img_ptr);
    ++ no_data;
    ++ img_ptr;
    ++ int_img_ptr;
  }
  ntuple_list out;
  out = lsd(lsd_image2);
  lsd_to_lines( vLines, out );
  lines_to_edges (edges, vLines);
  //lines_to_points(
#if 0
  mEdgeImage.fill(0);
  char filename[100];
  sprintf(filename,"%s/master_%d.xvy",folderName, master_file_id);
  std::ofstream out_xvy(filename);
  for( unsigned int i = 0; i < edges.virEdges.size(); ++i )
  {
    mEdgeImage[edges.virEdges[i]] = 255;
    out_xvy << CVD::vec(edges.virEdges[i]) << std::endl;
  }
  sprintf(filename,"%s/master_%d.png",folderName,master_file_id);
  //++master_file_id;
  CVD::img_save(mEdgeImage, filename);
#endif
  //std::cout << "no of lines = " << vLines.size() << std::endl;
  free_ntuple_list(out);
  line_time = CVD::timer.get_time();
}
  genDistanceTransform (vLines);
//  gen_Edgelets_with_Directional( mvEdgelets, mviLineIndexes, vLines, EDGELET_LENGTH );
  double edgelet_time = CVD::timer.get_time();
  double td_time = CVD::timer.get_time();

}

void EdgeExtraction::getLSDEdges (Edges & thisEdges)
{
   for (int i = 0; i < edges.virEdges.size(); i++){
      edges.virEdges[i].x *= 4;
      edges.virEdges[i].y *= 4;
      thisEdges.virEdges.push_back (edges.virEdges[i]);
   }
}

void EdgeExtraction::get_lines( std::vector<Line> & vLines, const Edges & edges )
{
  int max_LinkedEdges = edges.viEdgeIdxes.size(), startIdx, endIdx;
  for( int no_LinkedEdges = 1; no_LinkedEdges < max_LinkedEdges ; ++ no_LinkedEdges )
  {
    Line line;
    REAL_TYPE a,b,c,d;
    startIdx = edges.viEdgeIdxes[ no_LinkedEdges - 1 ];
    int maxIdx = edges.viEdgeIdxes[ no_LinkedEdges ] - 1;
    do 
    {
      REAL_TYPE max = 4.f;
      while( max > 2 )// 2 pixel deviation
      {
        endIdx = maxIdx;
        line.v2Start = vec(edges.virEdges[ startIdx ]);
        line.v2End = vec(edges.virEdges[ endIdx ]);
        a = line.v2End[0] - line.v2Start[0];
        b = line.v2Start[1] - line.v2End[1];
        c = line.v2End[1]*line.v2Start[0] - line.v2End[0]*line.v2Start[1];
        d = sqrt( a*a + b*b );
        max = -1.f;
        for( int i = startIdx; i <= endIdx; ++i )
        {
          const CVD::ImageRef & irPoint = edges.virEdges[i];
          REAL_TYPE rDistance = fabs( a*irPoint.y + b*irPoint.x + c );
          if( rDistance > max )
          {
            max = rDistance;
            maxIdx = i;
          }
        }
        max /= d;
      }
      line.v2Direction = line.v2End - line.v2Start;
      line.no_pixel = norm(line.v2Direction);
      if( line.no_pixel > 4 )
      {
//        line.rSlope = ( a != 0 ) ? -b/a: -b/0.01;
        line.rSlope = ( (line.v2Direction[0]) != 0 ) ? (line.v2Direction[1])/(line.v2Direction[0]) : (line.v2Direction[1])/0.01; 
        normalize(line.v2Direction);
        vLines.push_back(line);
      }
      startIdx = endIdx;
      maxIdx = edges.viEdgeIdxes[ no_LinkedEdges ] - 1;
    } while ( startIdx < maxIdx );
  }
}

void EdgeExtraction::genDistanceTransform(const std::vector<Line> &vLines, long oriSize)
{
  for( int i = 0; i < 11; ++ i )
  {
    mvEdgeOrientations[i].reset();
  }
  for( unsigned int i = 0; i < vLines.size(); ++i )
  {
    const Line & line = vLines[i];
    int binNo = static_cast<int>( (atan( line.rSlope )/M_PI + 0.5f )*10.f );
    mvEdgeOrientations[binNo][line.v2End] = 1;    
    for( int no_pixel = 0; no_pixel <= line.no_pixel; ++no_pixel )
    {
      TooN::Vector<2> v2Pose = line.v2Start + no_pixel*line.v2Direction;
        mvEdgeOrientations[binNo][v2Pose] = 1;      
    }
  }
  double distTransform1 = 0, distTransform2 = 0, distTransform3 = 0, distTransform4 = 0;
// Generate Edge Orientation Image 
  for( unsigned int i = 0; i < mvEdgeOrientationImages.size(); ++i )
  {
    distTransform1 = CVD::timer.get_time();
    GetTransformDist::compute( mvEdgeOrientationDistTransforms[i], mvEdgeOrientations[i]); 

    distTransform2 = CVD::timer.get_time();
    REAL_TYPE max = 0;
    int idx = 0;
    REAL_TYPE *data = mvEdgeOrientationDistTransforms[i].data;
    while( idx < oriSize )
    {
      if( max < *data )
        max = *data;
      ++data;
      ++idx;
    }
    distTransform3 = CVD::timer.get_time();
    data = mvEdgeOrientationDistTransforms[i].data;
    idx = 0;
    while( idx < oriSize )
    {
      *data /= max; // normalise;
      ++data;
      ++idx;
    }
    distTransform4 = CVD::timer.get_time();
  }
}

void EdgeExtraction::gen_Edgelets_with_Directional( std::vector<Edgelet> & vEdgelets, std::vector<int> & viLineIndexes, const std::vector<Line> & vLines, long oriSize, int iMaxLength)
{
  double local_start_time = CVD::timer.get_time();
  vEdgelets.clear();
  viLineIndexes.clear();
  for( int i = 0; i < 11; ++ i )
  {
    mvEdgeOrientations[i].reset();
  }
  REAL_TYPE halfMaxLength = static_cast<REAL_TYPE>(iMaxLength)*0.5f;
  for( unsigned int i = 0; i < vLines.size(); ++i )
  {
    const Line & line = vLines[i];
    int binNo = static_cast<int>( (atan( line.rSlope )/M_PI + 0.5f )*10.f );
    mvEdgeOrientations[binNo][line.v2End] = 1;    
    for( int no_pixel = 0; no_pixel <= line.no_pixel; ++no_pixel )
    {
      TooN::Vector<2> v2Pose = line.v2Start + no_pixel*line.v2Direction;
        mvEdgeOrientations[binNo][v2Pose] = 1;      
    }
    int no_segment = line.no_pixel/iMaxLength;// iMaxLength pixel per edgelet
    viLineIndexes.push_back( vEdgelets.size() );
    if( no_segment == 0 )
    {
      if( line.v2Start[0] < 0 || line.v2Start[1] < 0 || line.v2End[0] < 0 || line.v2End[1] < 0 )
        continue;
      Edgelet edgeLet;
      edgeLet.v2Pose = (line.v2Start + line.v2End)*0.5;
      edgeLet.v2Pose[0] = round(edgeLet.v2Pose[0] ); edgeLet.v2Pose[1] = round(edgeLet.v2Pose[1]);

      edgeLet.rSlope = line.rSlope;
      edgeLet.regionNo = 0;
      vEdgelets.push_back(edgeLet);
    }
    else
    {
      for( int j = 0; j < no_segment; ++ j )
      {
        Edgelet edgeLet;
        edgeLet.v2Pose = line.v2Start + (iMaxLength*j+halfMaxLength)*line.v2Direction;
        edgeLet.v2Pose[0] = round(edgeLet.v2Pose[0] ); edgeLet.v2Pose[1] = round(edgeLet.v2Pose[1]);

        edgeLet.rSlope = line.rSlope;
        edgeLet.regionNo = 0;
        vEdgelets.push_back(edgeLet);
      }
      if( norm(vEdgelets.back().v2Pose - line.v2End) > 1 )
      {
        if( line.v2End[0] < 0 || line.v2End[1] < 0 )
          continue;
        Edgelet edgeLet;
        edgeLet.v2Pose = ( line.v2Start + iMaxLength*no_segment*line.v2Direction + line.v2End )*0.5f;
        edgeLet.v2Pose[0] = round(edgeLet.v2Pose[0] ); edgeLet.v2Pose[1] = round(edgeLet.v2Pose[1]);
        if (edgeLet.v2Pose[0] > 0 && edgeLet.v2Pose[0] <= img_h2 && edgeLet.v2Pose[1] > 0 && edgeLet.v2Pose[1] <= img_w2)
        {
           edgeLet.rSlope = line.rSlope;
           edgeLet.regionNo = 0;
           vEdgelets.push_back(edgeLet);
        }
      }
    }
  }
  viLineIndexes.push_back( vEdgelets.size() );
  double edgelet_time = CVD::timer.get_time();
  double distTransform1 = 0, distTransform2 = 0, distTransform3 = 0, distTransform4 = 0;
// Generate Edge Orientation Image 
  for( unsigned int i = 0; i < mvEdgeOrientationImages.size(); ++i )
  {
    distTransform1 = CVD::timer.get_time();
    GetTransformDist::compute( mvEdgeOrientationDistTransforms[i], mvEdgeOrientations[i]); 

    distTransform2 = CVD::timer.get_time();
    REAL_TYPE max = 0;
    int idx = 0;
    REAL_TYPE *data = mvEdgeOrientationDistTransforms[i].data;
    while( idx < oriSize )
    {
      if( max < *data )
        max = *data;
      ++data;
      ++idx;
    }
    distTransform3 = CVD::timer.get_time();
    data = mvEdgeOrientationDistTransforms[i].data;
    idx = 0;
    while( idx < oriSize )
    {
      *data /= max; // normalise;
      ++data;
      ++idx;
    }
    distTransform4 = CVD::timer.get_time();
  }
}

void EdgeExtraction::lsd_to_lines( std::vector<Line> & vLines, ntuple_list lsd_format )
{
  vLines.reserve( lsd_format->size );
  for( unsigned int i = 0; i < lsd_format->size; ++i )
  {
    Line line;
    int linePose = i*lsd_format->dim;
    line.v2Start = TooN::makeVector( (lsd_format->values[linePose]), (lsd_format->values[linePose+1]) );
    if( line.v2Start[0] < 0 || line.v2Start[1] < 0 )//|| line.v2Start[0] >= lsd_image->xsize || line.v2Start[1] >= lsd_image->ysize )
      continue;
    line.v2End = TooN::makeVector( (lsd_format->values[linePose+2]), (lsd_format->values[linePose+3]));
    if( line.v2End[0] < 0 || line.v2End[1] < 0 )//|| line.v2End[0] >= lsd_image->xsize || line.v2End[1] >= lsd_image->ysize )
      continue;
    line.v2Direction = line.v2End - line.v2Start;
    line.no_pixel = norm(line.v2Direction);
    line.rSlope = ( (line.v2Direction[0]) != 0 ) ? (line.v2Direction[1])/(line.v2Direction[0]) : (line.v2Direction[1])/0.01; 
    normalize( line.v2Direction );
    vLines.push_back(line);
  }
}

void EdgeExtraction::lines_to_edges (Edges & edges, std::vector<Line> vLines)
{
   edges.virEdges.clear();
   edges.virEdges.reserve (vLines.size()*EDGELET_LENGTH);
  for( unsigned int no_line = 0; no_line < vLines.size(); ++ no_line )
  {
    Line & line = vLines[no_line];
    for( int idx = 0; idx < line.no_pixel; ++idx )
    {
      TooN::Vector<2> v2 = line.v2Start + static_cast<REAL_TYPE>(idx)*line.v2Direction;
      if( v2[0] >= 0 && v2[1] >0  && v2[0] < lsd_image2->xsize && v2[1] < lsd_image2->ysize )
        edges.virEdges.push_back( CVD::ir(v2) );
    }
  }

}

void EdgeExtraction::project_edgelets( std::vector<Edgelet> & vProjectedEdgelets, const std::vector<Edgelet> & vEdgelets, const TooN::Matrix<3> & h )
{
  vProjectedEdgelets.clear();
  vProjectedEdgelets.reserve( vEdgelets.size() );
  std::vector<Edgelet>::const_iterator beginIter = vEdgelets.begin();
  std::vector<Edgelet>::const_iterator endIter = vEdgelets.end();

  for( ; beginIter != endIter; ++beginIter )
  {
    TooN::Vector<3> v3Pose = unproject( beginIter->v2Pose );
    TooN::Vector<3> v3Project = h*v3Pose;
    TooN::Vector<3> v3Slope = h*TooN::makeVector( 1.f, beginIter->rSlope, 0 );
    Edgelet ed;
    ed.v2Pose = project( v3Project );
    if( ed.v2Pose[0] > 0 && ed.v2Pose[0] < mirImageSize.x && ed.v2Pose[1] > 0 && ed.v2Pose[1] < mirImageSize.y )
    {
      vProjectedEdgelets.push_back( ed );
      vProjectedEdgelets.back().rSlope = (v3Project[1]*v3Slope[2] - v3Project[2]*v3Slope[1])/(v3Project[0]*v3Slope[2] - v3Project[2]*v3Slope[0]);
    }
  }
}

void EdgeExtraction::project_edges( Edges & vProjectedEdges, std::vector<TooN::Vector<2> >& vv2Pnts, const Edges & vEdges, const TooN::Matrix<3> & h )
{
  std::vector<CVD::ImageRef> vEdgelets = vEdges.virEdges;
  std::vector<CVD::ImageRef> vProjectedEdgelets = vProjectedEdges.virEdges;
  vProjectedEdgelets.clear();
  vProjectedEdgelets.reserve( vEdgelets.size() );
  std::vector<CVD::ImageRef>::const_iterator beginIter = vEdgelets.begin();
  std::vector<CVD::ImageRef>::const_iterator endIter = vEdgelets.end();

 // TooN::Vector<2> v2Pose;
  for( ; beginIter != endIter; ++beginIter )
  {
    TooN::Vector<2> v2Pose;
    v2Pose[0] = beginIter->x;
    v2Pose[1] = beginIter->y;
    TooN::Vector<3> v3Pose = unproject(v2Pose );
    TooN::Vector<3> v3Project = h*v3Pose;
    v2Pose = project(v3Project);
    if( v2Pose[0] > 0 && v2Pose[0] < mirImageSize.x && v2Pose[1] > 0 && v2Pose[1] < mirImageSize.y )
    {  vProjectedEdgelets.push_back(CVD::ir(v2Pose));
    }
    v2Pose[0] = v2Pose[0]*4;
    v2Pose[1] = v2Pose[1]*4;
    vv2Pnts.push_back (v2Pose);
  }
}

double EdgeExtraction::getConfidenceScore (View v, TooN::Matrix<3> h)
{
    std::vector<Edgelet> vProjectedEdgelets;
    project_edgelets( vProjectedEdgelets, v.vEdgelets, h );
    std::vector<REAL_TYPE> vrDistanceSum( vProjectedEdgelets.size() );

    for( unsigned int no_edgelet = 0; no_edgelet < vProjectedEdgelets.size(); ++ no_edgelet )
    {
         int binNo = static_cast<int>( (atan( vProjectedEdgelets[no_edgelet].rSlope )/M_PI + 0.5f ) * 10.f );
         vrDistanceSum[no_edgelet] = mvEdgeOrientationDistTransforms[binNo][vProjectedEdgelets[no_edgelet].v2Pose];
    }
    std::sort( vrDistanceSum.begin(), vrDistanceSum.end() );
    REAL_TYPE rMedian = vrDistanceSum[vrDistanceSum.size()/2];
    return rMedian;
}

double EdgeExtraction::getConfidenceScoreModified (View v, TooN::Matrix<3> h, std::vector<Edgelet>& projectedEdgelets)
{
    projectedEdgelets.clear();
    project_edgelets( projectedEdgelets, v.vEdgelets, h );

    std::vector<REAL_TYPE> vrDistanceSum( projectedEdgelets.size() );
    REAL_TYPE sum = 0.0;

    for( unsigned int no_edgelet = 0; no_edgelet < projectedEdgelets.size(); ++ no_edgelet )
    {
         int binNo = static_cast<int>( (atan( projectedEdgelets[no_edgelet].rSlope )/M_PI + 0.5f ) * 10.f );
         vrDistanceSum[no_edgelet] = mvEdgeOrientationDistTransforms[binNo][projectedEdgelets[no_edgelet].v2Pose];
         sum += vrDistanceSum[no_edgelet];
    }

    // Add max score for each unprojected edgelet
    for (unsigned int i = 0; i < v.vEdgelets.size() - projectedEdgelets.size(); ++i)
    {
		vrDistanceSum.push_back(1.0);
		sum += 1.0;
    }

    std::sort( vrDistanceSum.begin(), vrDistanceSum.end() );

    REAL_TYPE rMedian = vrDistanceSum[vrDistanceSum.size()/2];
    REAL_TYPE rMean = sum / vrDistanceSum.size();

    return rMedian;
}

void EdgeExtraction::load_codebook(std::string filename)
{
  FILE * if_codebook = fopen( filename.c_str(), "r+b" );
  mvObjDescriptorClasses.clear();
  unsigned int class_size;
  int img_read;
  img_read = fread( &class_size, 1, sizeof(class_size), if_codebook );
  mvObjDescriptorClasses.resize( class_size );
  for( unsigned int no_class = 0; no_class < class_size; ++no_class )
  {
    ObjectsTemplate & mObjectsTemplate = mvObjDescriptorClasses[no_class].mObjectsTemplate;
    std::vector<CodeBookElement> & mvCodeBook = mvObjDescriptorClasses[no_class].mvCodeBook;
    TooN::Vector<4> & mv4ChainAngles = mvObjDescriptorClasses[no_class].mv4ChainAngles;
    for( int i = 0; i < 4; ++i ) 
    {
        img_read = fread( &mv4ChainAngles[i], 1, sizeof(mv4ChainAngles[i]), if_codebook );
    }
    unsigned int codebook_size;
    img_read = fread( &codebook_size, 1, sizeof(codebook_size), if_codebook );
    mvCodeBook.resize( codebook_size );
    for( unsigned int i = 0; i < codebook_size; ++i )
    {
      img_read = fread( &mvCodeBook[i].rAngle_0, 1, sizeof(REAL_TYPE), if_codebook );
      img_read = fread( &mvCodeBook[i].rAngle_1, 1, sizeof(REAL_TYPE), if_codebook );
      img_read = fread( &mvCodeBook[i].rRelativeDist_1, 1, sizeof(REAL_TYPE), if_codebook );
      img_read = fread( &mvCodeBook[i].rAngle_2, 1, sizeof(REAL_TYPE), if_codebook );
      img_read = fread( &mvCodeBook[i].rRelativeDist_2, 1, sizeof(REAL_TYPE), if_codebook );
      img_read = fread( &mvCodeBook[i].rAngle_3, 1, sizeof(REAL_TYPE), if_codebook );
      img_read = fread( &mvCodeBook[i].rRelativeDist_3, 1, sizeof(REAL_TYPE), if_codebook );
      img_read = fread( &mvCodeBook[i].iIdx_0, 1, sizeof(unsigned short), if_codebook );
      img_read = fread( &mvCodeBook[i].iIdx_1, 1, sizeof(unsigned short), if_codebook );
      img_read = fread( &mvCodeBook[i].iIdx_2, 1, sizeof(unsigned short), if_codebook );
      img_read = fread( &mvCodeBook[i].iIdx_3, 1, sizeof(unsigned short), if_codebook );
      img_read = fread( &mvCodeBook[i].iIdx_4, 1, sizeof(unsigned short), if_codebook );
      img_read = fread( &mvCodeBook[i].iView_ID, 1, sizeof(unsigned short), if_codebook );
      img_read = fread( &mvCodeBook[i].iObject_ID, 1, sizeof(unsigned short), if_codebook );
      img_read = fread( &mvCodeBook[i].obj_correct_ID, 1, sizeof(unsigned short), if_codebook );
    }
    std::vector<CVD::ImageRef> & mvirFirstAngleIndexes = mvObjDescriptorClasses[no_class].mvirFirstAngleIndexes;
    mvirFirstAngleIndexes.resize(64);
    for( int i = 0; i < 64; ++i )
    {
      img_read = fread( &mvirFirstAngleIndexes[i].x, 1, sizeof(mvirFirstAngleIndexes[i].x), if_codebook );
      img_read = fread( &mvirFirstAngleIndexes[i].y, 1, sizeof(mvirFirstAngleIndexes[i].y), if_codebook );
    }
    std::vector< std::vector<int> > & mvirSecondIndexes = mvObjDescriptorClasses[no_class].mvirSecondIndexes;
    mvirSecondIndexes.resize(163840);
    for (int i = 0; i < 163840; i++) 
    {
       int sx;
       img_read = fread (&sx, 1, sizeof (sx), if_codebook);
       std::vector<int> &vectors = mvirSecondIndexes[i];
       vectors.resize(sx);
       for (int j = 0; j < sx; j++)
       {
          fread (&vectors[j], 1, sizeof (vectors[j]), if_codebook);
       }
    }

    unsigned int no_objects;
    img_read = fread( &no_objects, 1, sizeof(no_objects), if_codebook );
    mObjectsTemplate.vObjectElements.resize(no_objects);
    mObjectsTemplate.iNoObjects = no_objects;
    std::cout << "LOADING OBJECTS OF SIZE " << no_objects << std::endl;
    for( unsigned int obj_ID = 0; obj_ID < no_objects; ++ obj_ID )
    {
      ObjectElement & obj = mObjectsTemplate.vObjectElements[obj_ID];
      obj.iObject_ID = obj_ID;
      img_read = fread( &obj.obj_correct_ID, 1, sizeof (obj.obj_correct_ID), if_codebook);
      img_read = fread( &obj.iNoViews, 1, sizeof(obj.iNoViews), if_codebook );
      obj.vViews.resize( obj.iNoViews );
      for( int view_ID = 0; view_ID < obj.iNoViews; ++view_ID )
      {
        View & view = obj.vViews[view_ID];
        view.iView_ID = view_ID;
        view.iObject_ID = obj_ID;
        std::vector<Edgelet> & vEdgelets = view.vEdgelets;
        unsigned int no_edgelets, no_edgelets2;
        img_read = fread( &no_edgelets, 1, sizeof(no_edgelets), if_codebook );
        vEdgelets.resize( no_edgelets );
        for( unsigned int edg_ID = 0; edg_ID < no_edgelets; ++ edg_ID )
        {
          img_read = fread( &vEdgelets[edg_ID].v2Pose[0], 1, sizeof( vEdgelets[edg_ID].v2Pose[0] ), if_codebook );
          img_read = fread( &vEdgelets[edg_ID].v2Pose[1], 1, sizeof( vEdgelets[edg_ID].v2Pose[1] ), if_codebook );
          img_read = fread( &vEdgelets[edg_ID].rSlope, 1, sizeof(REAL_TYPE), if_codebook );
        }
        img_read = fread( &no_edgelets2, 1, sizeof(no_edgelets2), if_codebook);
        std::vector<CVD::ImageRef> & allEdges = view.allEdges.virEdges;
        allEdges.resize (no_edgelets2);
        for (unsigned int edg_ID = 0; edg_ID < allEdges.size(); edg_ID++){
           img_read = fread(&allEdges[edg_ID].y, 1, sizeof (allEdges[edg_ID].y), if_codebook);
           img_read = fread(&allEdges[edg_ID].x, 1, sizeof (allEdges[edg_ID].x), if_codebook);
        }

      }
    }
    // if a new object is added (no prior detection), it adds it at the end
    std::cout <<"Size of Code Book : "<< mvCodeBook.size() << std::endl;
  }
  fclose( if_codebook );
}
    
View EdgeExtraction::getView (int chainNo, int classNo, int viewNo) const
{
   if (chainNo >= mvObjDescriptorClasses.size())
   {
       std::cout << "ERROR: chainNo more than the maximum" << std::endl;
       exit(-1);
   }
   if (classNo >= mvObjDescriptorClasses[chainNo].mObjectsTemplate.iNoObjects){
      std::cout << "ERROR: classNo more than the maximum" << std::endl;
      exit(-1);
   }
   if (viewNo >= mvObjDescriptorClasses[chainNo].mObjectsTemplate.vObjectElements[classNo].iNoViews){
       std::cout << "ERROR: viewNo more than the maximum" << std::endl;
       exit(-1);
   }
   return mvObjDescriptorClasses[chainNo].mObjectsTemplate.vObjectElements[classNo].vViews[viewNo];
}

LLImagePyramid<CVD::byte> EdgeExtraction::getImageFromKinect (const sensor_msgs::ImageConstPtr& msg)
{
    rgb_image = cv_bridge::toCvShare(msg, sensor_msgs::image_encodings::MONO8);
    CVD::BasicImage<CVD::byte> cvd_rgb_image_b (rgb_image->image.data, CVD::ImageRef(rgb_image->image.cols, rgb_image->image.rows));
    convert_image (cvd_rgb_image_b, mGrayImageF);
    copy (mGrayImageF, mGrayImage, mGrayImage.size(), CVD::ImageRef(1,1));
    mLLImagePyramid.SetImageAtLevelZero( mGrayImage );
    mLLImagePyramid.Reset();
    return mLLImagePyramid;
}

void EdgeExtraction::saveDistanceTransformImages(const std::string& directory)
{
	for( unsigned int i = 0; i < mvEdgeOrientationImages.size(); ++i )
	{
		std::string name = directory + "edgeOrientation" + boost::lexical_cast<std::string>(i) + ".png";
		CVD::BasicImage<REAL_TYPE> image(mvEdgeOrientations[i].data, mvEdgeOrientations[i].size());
		CVD::img_save(image, name);
	}

	for (unsigned int i = 0; i < mvEdgeOrientationDistTransforms.size(); ++i)
	{
		std::string name = directory + "distanceTransform" + boost::lexical_cast<std::string>(i) + ".png";
		CVD::BasicImage<REAL_TYPE> image(mvEdgeOrientationDistTransforms[i].data, mvEdgeOrientationDistTransforms[i].size());
		CVD::img_save(image, name);
	}
}
