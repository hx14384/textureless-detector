#include "EdgeObjectDetection.h"
#include <cvd/vector_image_ref.h>
#include <cvd/timer.h>
#include <fcntl.h>
#include <algorithm>
#include "GetTransformDist.h"
#include "BitMacros.h"

#define PRL_THRESHOLD 0.0524
//#define ECM
//#define LSD
image_double lsd_image;
int data_size;
int master_file_id=1;
int template_file_id=1;
int projected_file_id=1;
unsigned short obj_correct_ID=0;
double start_time = 0;
bool isLSD = true;

#define TEST_RESTART
const REAL_TYPE RESTART_TIME = 5.f;
const REAL_TYPE SAMPLING_FACTOR = 2.f;
const REAL_TYPE EPS_DIR    = 0.04; // 0.1 without depth
REAL_TYPE EPS_DIR_TEST = 0.08; //0.06 ICCV 0.15 without depth
REAL_TYPE EPS_DIST   = 0.15; //0.06 ICCV 0.15 without depth
REAL_TYPE EPS_ANGLE  = 0.15; //0.06 ICCV 0.15 without depth
const int FIRST_DIST_TRAIN = 100; //100

int FIRST_DIST_TEST  = 20; //60
//MULTISCALE
const int MAX_MULTISCALE_LEVEL = 1;
int MAX_DIST_TEST [MAX_MULTISCALE_LEVEL+1]  = {60, 100};
int MIN_DIST_TEST [MAX_MULTISCALE_LEVEL+1]  = {2, 50};
const int EDGELET_LENGTH_TRAIN = 15; //5
//int stepNo [100000];
//double timestamp [100000];
const int NumberOfSteps = 14;
double timestampsteps [NumberOfSteps];
int counterSteps [NumberOfSteps];
int EDGELET_LENGTH_TEST = 5; // was 5
const REAL_TYPE MAX_DT = 6; //    = 0.04; 15 without depth

const int LINE_MIN_LEN = 4;
const int EDGELET_MIN_LEN = 5;



//to write to a folder
char* folderName;
char outputFileName [100];

const int fstRing[8][2]={{-1,-1},{-1,0},{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1}};
const int sndRing[16][2]={{-2,-2},{-2,-1},{-2,0},{-2,1},{-2,2},{-1,2},{0,2},{1,2},{2,2},{2,1},{2,0},{2,-1},{2,-2},{1,-2},{0,-2},{-1,-2}};

int img_h, img_w, cellCount, cellRow;

/*typedef std::pair<std::vector<EdgeLinkQuestion>::iterator , int> my_pair;
bool sort_pred( const my_pair& left, const my_pair& right )
{
  return left.second < right.second;
}*/

EdgeObjectDetection::EdgeObjectDetection( CVD::ImageRef irSize, char* fN, char* ofN, bool isS, bool lsdV, double maxTime, int isr) : mbEmpty (true), miTotalObjectViews(0), mvEdgeOrientationImages(11), mvEdgeOrientationDistTransformImages(11), mvEdgeOrientations(11), mvEdgeOrientationDistTransforms(11), mEdgeImage( irSize ), mirImageSize( irSize), viewAdded(false), learntObjNo (0), learntViewNo (0), learntConstNo (0)
{

  char folderString [100];
  isStar = isS;
  folderName = fN;
  img_h = irSize.y;
  isMask = false;
  isDepth = false;
  imgSizeRatio = isr;
  lineMaxLen = 4;
  // in case no detection, this sets it for an empty codebook
  latestObjID = 11;
  latest_correctObjID = 11;
  img_w = irSize.x;
  cellCount = (img_h/10*img_w/10);
  ecm = (EdgeCellMatrix*) malloc (img_h*img_w*sizeof(EdgeCellMatrix));
  cellRow = img_w/10;
 for (int i = 0; i < cellCount; i++)
 {
    ecm[i].cellNo = i;
    ecm[i].count = 0;
    ecm[i].v2Pose[0] = (i / cellRow) * 10 + 1;
    ecm[i].v2Pose[1] = (i % cellRow) * 10 + 1;
    ecm[i].cellEdgelets.reserve (40);
 }

//  generateOfflineCellDirectStruct(240, 320, 10, 10, 0.1);
//  exit(-1);
  isLSD = lsdV;
  sprintf(outputFileName, ofN);
  sprintf(folderString, "mkdir %s", folderName);
  system (folderString);
  MAX_TIME = maxTime;

  for( int i = 0; i < 11; ++i )
  {
    mvEdgeOrientationImages[i].resize( irSize );
    mvEdgeOrientations[i].resize(irSize);
    mvEdgeOrientations[i].reset();
    mvEdgeOrientationDistTransformImages[i].resize( irSize );
    mvEdgeOrientationDistTransforms[i].resize(irSize);
    mvEdgeOrientationDistTransforms[i].reset();
  }
  mvEdgeDistanceTransformImage.resize( irSize );
  mvEdgeDistanceTransformImageIntegral.resize( irSize );
 lsd_image = new_image_double( irSize.x, irSize.y );
 data_size = lsd_image->xsize * lsd_image->ysize;
 mvCheckedChains.reserve (200000);
 //EdgeCellMatrix ecmM [img_h*img_w];
 //ecm = ecmM;
 #ifdef ECM
 cds = (CellDirectStruct*) malloc (768*sizeof(CellDirectStruct));
 char cdsFilename[100];
 sprintf(cdsFilename, "cellOfflineStruct3_%d_%d.dat", img_h, img_w);
 load_cellDirectStructure (cdsFilename);
#endif
}

EdgeObjectDetection::~EdgeObjectDetection()
{
   free_image_double(lsd_image);
   free(ecm);
   #ifdef ECM
   free(cds);
   #endif
}

bool comp_testdirection( const TestDirection & left, REAL_TYPE value )
{
  return left.rDir_1 < value;
}

bool comp_pair_r_i ( const std::pair<REAL_TYPE, int> & left, const std::pair<REAL_TYPE, int> & right )
{
  return left.first < right.first;
}

bool comp_lower_pair( const std::pair<REAL_TYPE, int> & left, REAL_TYPE value )
{
  return left.first < value;
}

bool comp_upper_pair( REAL_TYPE value, const std::pair<REAL_TYPE, int> & right )
{
  return value < right.first;
}

bool comp_lower_codebook( const CodeBookElement & left, REAL_TYPE value )
{
  return left.rAngle_0 < value;
}

bool comp_upper_codebook( REAL_TYPE value, const CodeBookElement & right )
{
  return value < right.rAngle_0;
}

bool EdgeObjectDetection::detect_from_previous(int subtract_time, int previousObjNo)
{
  bool obj_detected = false;
  double current_time;
  int beginIdx = 0;
  int endIdx;

  std::vector<MatchedPair> v5MatchedPairs(5);
  
  int prevObjCorrectNo;
  if (mvDetectedObjects_previousFrame[previousObjNo].iObject_ID <= 2)
      prevObjCorrectNo = mvDetectedObjects_previousFrame[previousObjNo].iObject_ID - 1; 
  else
      prevObjCorrectNo = mvDetectedObjects_previousFrame[previousObjNo].iObject_ID - 3;
  //std::cout << "we are here with prevObjCorrectNo " << previousObjNo << " " << prevObjCorrectNo << std::endl;
  TooN::Vector<4> & mv4ChainAngles = mvIndObjDescriptorClasses[prevObjCorrectNo].mv4ChainAngles;
  std::vector<CodeBookElement> & mvCodeBook = mvIndObjDescriptorClasses[prevObjCorrectNo].mvCodeBook;
  std::vector<CVD::ImageRef> & mvirFirstAngleIndexes = mvIndObjDescriptorClasses[prevObjCorrectNo].mvirFirstAngleIndexes;
  std::vector< std::vector<int> > & mvirSecondIndexes = mvIndObjDescriptorClasses[prevObjCorrectNo].mvirSecondIndexes;
  //std::cout << "correctly read the first phase " << mvCodeBook.size() << " " << mvirSecondIndexes.size() << " " << mvTestDirectionsPreviousObjects[previousObjNo].size() << std::endl;

  std::vector<bool> vBinTestDirections(mvTestDirectionsPreviousObjects[previousObjNo].size(),false);
  std::vector<int> vRandomSamples(mvTestDirectionsPreviousObjects[previousObjNo].size());
  for( unsigned int i = 0; i < mvTestDirectionsPreviousObjects[previousObjNo].size(); ++ i ) vRandomSamples[i]=i;
  endIdx = mvTestDirectionsPreviousObjects[previousObjNo].size() - 1;
  std::vector<TestDirection> & mvTD = mvTestDirectionsPreviousObjects[previousObjNo];

  mvCheckedChains.clear();
  double delta_time;
  double step_time = 0;

  RESTART:
  std::random_shuffle( vRandomSamples.begin(), vRandomSamples.end() );

  for( ; beginIdx <= endIdx; ++ beginIdx )
  {
    delta_time = CVD::timer.get_time() - start_time; 
    if( delta_time > MAX_TIME ) goto END;
#ifdef TEST_RESTART    
    if( delta_time > step_time + RESTART_TIME )
    {
      //std::cout << delta_time <<" " << step_time << std::endl;
      step_time += RESTART_TIME;
      goto RESTART;
    }
#endif   
    if( vBinTestDirections[vRandomSamples[beginIdx]] )
      continue;
    TestDirection & firstTestDirectionPairs = mvTD[vRandomSamples[beginIdx]];
    if(!( (fabs( firstTestDirectionPairs.rDir_1-mv4ChainAngles[0] ) < EPS_DIR_TEST || fabs( firstTestDirectionPairs.rDir_1+mv4ChainAngles[0] ) < EPS_DIR_TEST ) && firstTestDirectionPairs.rDist < FIRST_DIST_TEST ) ){
        continue;
      }
     /*   CheckedChains ch;
        ch.chains.reserve (2);
         TooN::Vector<2> v1, v2;
         v1[0] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[0]*imgSizeRatio;
         v1[1] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[1]*imgSizeRatio;
         v2[0] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
         v2[1] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
        ch.chains.push_back (v1);
        ch.chains.push_back (v2);
        mvCheckedChains.push_back (ch);*/
      // looking for a second pair that fst.e2 == snd.e1
      int beginSndIndex = 0;
      if( beginSndIndex == -1 ) continue;
      int endSndIndex = mvTD.size();
      // searching for a first pair inside codebook.
      REAL_TYPE rLower_Angle = firstTestDirectionPairs.rAngle - EPS_ANGLE;
      REAL_TYPE rUpper_Angle = firstTestDirectionPairs.rAngle + EPS_ANGLE;
      returnedEdgelets.clear();
      returnedCB.clear();
      for( ; beginSndIndex <= endSndIndex; ++beginSndIndex ) // search through second pair
      {
        delta_time = CVD::timer.get_time() - start_time; 

        if( delta_time > MAX_TIME ) goto END;

        #ifdef TEST_RESTART        
        if( delta_time > step_time + RESTART_TIME )
        {
          step_time += RESTART_TIME;
          goto RESTART;
        }
        #endif        

        if( vBinTestDirections[beginSndIndex] )
          continue;
        const TestDirection & sndTestDirectionPairs = mvTD[beginSndIndex];
        if(firstTestDirectionPairs.iIdx_1 == sndTestDirectionPairs.iIdx_0 && fabs( calculateAngle(sndTestDirectionPairs.rDir_2,firstTestDirectionPairs.rDir_2) - mv4ChainAngles[1] ) < EPS_DIR_TEST ) // second edge lin k
        {
        // Check with codebook betweeen lowerCodeBook and upperCodeBook 
          std::vector< std::vector<CodeBookElement>::iterator > vSecondPassCodeBook;
          // trying to read from the codebook

          int indexNo = 64.f*(firstTestDirectionPairs.rAngle/M_PI + 0.5f); if( indexNo < 0 ) indexNo = 0;
          int indexNo2 = 64.f*(sndTestDirectionPairs.rAngle/M_PI + 0.5f); if( indexNo2 < 0 ) indexNo2 = 0;
          int indexNo3;
          REAL_TYPE ratio = sndTestDirectionPairs.rDist/firstTestDirectionPairs.rDist;
          if (ratio < 1)
              indexNo3 = int(ratio/0.1);
          else if (ratio > 1 && ratio < 29.5)
              indexNo3 = int (ratio)+9;
          else
              indexNo3 = 38;
          std::vector<int> & thisList = mvirSecondIndexes[(indexNo)*2496+(indexNo2)*39+indexNo3];
          if (((indexNo)*2496+(indexNo2)*39+indexNo3) > mvirSecondIndexes.size())
          { 
             std::cout << "trying to reserve 1 " << mvirSecondIndexes.size() << " " << ((indexNo)*2496+(indexNo2)*39+indexNo3) << " " << thisList.size() << std::endl;
             std::cout << indexNo << " " << indexNo2 << " " << indexNo3 << std::endl;
             std::cout << sndTestDirectionPairs.rAngle << std::endl;
             exit (-1);
          }
          vSecondPassCodeBook.reserve( thisList.size() );
          for(int k = 0; k < thisList.size(); k++)
          {
            std::vector<CodeBookElement>::iterator lx = mvCodeBook.begin()+thisList[k];
            if( fabs( (*lx).rAngle_1 - sndTestDirectionPairs.rAngle ) < EPS_ANGLE 
               && fabs( (*lx).rRelativeDist_1 - sndTestDirectionPairs.rDist/firstTestDirectionPairs.rDist ) < EPS_DIST )
                 vSecondPassCodeBook.push_back( lx );
          }     


          if( vSecondPassCodeBook.size() == 0 )
          {
             continue; 
          }
          int beginTrdIndex = 0;
          if( beginTrdIndex == -1 ) continue;
          int endTrdIndex = mvTD.size();
          
          for( ; beginTrdIndex <= endTrdIndex; ++beginTrdIndex ) // search through third pair
          {
            delta_time = CVD::timer.get_time() - start_time;

            if( delta_time > MAX_TIME ) goto END;
#ifdef TEST_RESTART            
            if( delta_time > step_time + RESTART_TIME )
            {
              step_time += RESTART_TIME;
              goto RESTART;
            }
#endif            
            if( vBinTestDirections[beginTrdIndex] )
              continue;
            TestDirection & trdTestDirectionPairs = mvTD[beginTrdIndex];
            if(sndTestDirectionPairs.iIdx_1 == trdTestDirectionPairs.iIdx_0 &&  fabs( calculateAngle(trdTestDirectionPairs.rDir_2, sndTestDirectionPairs.rDir_2) - mv4ChainAngles[2] ) < EPS_DIR_TEST ) // third edge link
            {
            // Check with codebook inside vSecondPassCodeBook;
              std::vector< std::vector<CodeBookElement>::iterator > vThirdPassCodeBook;
              vThirdPassCodeBook.reserve( vSecondPassCodeBook.size() );
              for( unsigned int snd_pass = 0; snd_pass < vSecondPassCodeBook.size(); ++snd_pass )
              {
                if( fabs( vSecondPassCodeBook[snd_pass]->rAngle_2 - trdTestDirectionPairs.rAngle ) < EPS_ANGLE
                   && fabs( vSecondPassCodeBook[snd_pass]->rRelativeDist_2 - trdTestDirectionPairs.rDist/sndTestDirectionPairs.rDist ) < EPS_DIST)
                  vThirdPassCodeBook.push_back( vSecondPassCodeBook[snd_pass] );
              }
              if( vThirdPassCodeBook.size() == 0 ) continue;

              int beginFrtIndex = 0;
              if( beginFrtIndex == -1 ) continue;
              int endFrtIndex = mvTD.size();

              for( ; beginFrtIndex <= endFrtIndex; ++beginFrtIndex ) // search through fourth pair
              {
                delta_time = CVD::timer.get_time() - start_time;
 
                if( delta_time > MAX_TIME ) goto END;
#ifdef TEST_RESTART                
                if( delta_time > step_time + RESTART_TIME )
                {
                  step_time += RESTART_TIME;
                  goto RESTART;
                }
#endif                
                if( vBinTestDirections[beginFrtIndex] )
                  continue;
                TestDirection & fourtTestDirectionPairs = mvTD[beginFrtIndex];
                if(trdTestDirectionPairs.iIdx_1 == fourtTestDirectionPairs.iIdx_0 &&  fabs( calculateAngle(fourtTestDirectionPairs.rDir_2,trdTestDirectionPairs.rDir_2) - mv4ChainAngles[3] ) < EPS_DIR_TEST ) // fourth edge link
                {
                  for( unsigned int trd_pass = 0; trd_pass < vThirdPassCodeBook.size(); ++trd_pass )
                  {
                    delta_time = CVD::timer.get_time() - start_time; 
                    if( delta_time > MAX_TIME ) goto END;
#ifdef TEST_RESTART                    
                    if( delta_time > step_time + RESTART_TIME )
                    {
                      step_time += RESTART_TIME;
                      goto RESTART;
                    }
#endif                 
                    if( fabs( vThirdPassCodeBook[trd_pass]->rAngle_3 - fourtTestDirectionPairs.rAngle ) < EPS_ANGLE 
                       && fabs( vThirdPassCodeBook[trd_pass]->rRelativeDist_3 - fourtTestDirectionPairs.rDist/trdTestDirectionPairs.rDist ) < EPS_DIST )
                    { 
                      CodeBookElement & cd = *(vThirdPassCodeBook[trd_pass]);
                      View & view = mvIndObjDescriptorClasses[prevObjCorrectNo].mObjectsTemplate.vObjectElements[ cd.iObject_ID ].vViews[ cd.iView_ID ];
                      std::vector<Edgelet> & vEdgelets = view.vEdgelets;
                      Edges & allEdges = view.allEdges;
                      v5MatchedPairs[0].v2FstView = vEdgelets[cd.iIdx_0].v2Pose; v5MatchedPairs[0].v2SndView = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose;
                      v5MatchedPairs[1].v2FstView = vEdgelets[cd.iIdx_1].v2Pose; v5MatchedPairs[1].v2SndView = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose;
                      v5MatchedPairs[2].v2FstView = vEdgelets[cd.iIdx_2].v2Pose; v5MatchedPairs[2].v2SndView = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose;
                      v5MatchedPairs[3].v2FstView = vEdgelets[cd.iIdx_3].v2Pose; v5MatchedPairs[3].v2SndView = mvEdgelets[trdTestDirectionPairs.iIdx_1].v2Pose;
                      v5MatchedPairs[4].v2FstView = vEdgelets[cd.iIdx_4].v2Pose; v5MatchedPairs[4].v2SndView = mvEdgelets[fourtTestDirectionPairs.iIdx_1].v2Pose;
                      // add the chain
                   /*   TooN::Vector<2> v1, v2, v3, v4, v5;
                      v1[0] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[0]*imgSizeRatio;
                      v1[1] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[1]*imgSizeRatio;
                      v2[0] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
                      v2[1] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
                      v3[0] = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
                      v3[1] = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
                      v4[0] = mvEdgelets[trdTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
                      v4[1] = mvEdgelets[trdTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
                      v5[0] = mvEdgelets[fourtTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
                      v5[1] = mvEdgelets[fourtTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
                      CheckedChains ch;
                      ch.chains.reserve (5);
                      ch.chains.push_back (v1);
                      ch.chains.push_back (v2);
                      ch.chains.push_back (v3);
                      ch.chains.push_back (v4);
                      ch.chains.push_back (v5);
                      mvCheckedChains.push_back (ch);*/
                      TooN::Matrix<3> h;
                      if ( mHomography.estimateHomography_group_of_5(h, v5MatchedPairs) )
                      {
                        //std::cout << "AND HOMOGRAPHY ESTIMATED " << std::endl;
                        std::vector<Edgelet> vProjectedEdgelets;
                        project_edgelets( vProjectedEdgelets, vEdgelets, h );
                        std::vector<REAL_TYPE> vrDistanceSum( vProjectedEdgelets.size() );

                        for( unsigned int no_edgelet = 0; no_edgelet < vProjectedEdgelets.size(); ++ no_edgelet )
                        {
                            int binNo = static_cast<int>( (atan( vProjectedEdgelets[no_edgelet].rSlope )/M_PI + 0.5f ) * 10.f );
                            vrDistanceSum[no_edgelet] = mvEdgeOrientationDistTransforms[binNo][vProjectedEdgelets[no_edgelet].v2Pose];
                        }
                        std::sort( vrDistanceSum.begin(), vrDistanceSum.end() );
                        REAL_TYPE rMedian = vrDistanceSum[vrDistanceSum.size()/2];
                        
                        if( rMedian < (MAX_DT*2) )
                        {
                          project_edgelets( vProjectedEdgelets, vEdgelets, h, false );

                          TooN::Matrix<3> newH;
                          std::vector<int> viFoundEdgeIndexes;
                          std::vector<int> viFullIndexes;
                          std::vector<MatchedPair> vMatchedPairs;
                          double angle;
                          REAL_TYPE rErr1 = iterativeClosestEdges2( newH, viFullIndexes, viFoundEdgeIndexes, vProjectedEdgelets, mvEdgelets, vMatchedPairs, 25.f, 30.f, 0.5f);
                          //angle = estimateEdgeletRotation (vEdgelets, viFullIndexes, mvEdgelets);

                          std::vector<bool> vbBins( mvEdgelets.size(), true );
                          if( viFoundEdgeIndexes.size() > 5 )
                          {
                            int no_unique = 0;
                            std::vector<int>::iterator iIter = viFoundEdgeIndexes.begin();
                            std::vector<int>::iterator iEndIter = viFoundEdgeIndexes.end();
                            while( iIter != iEndIter )
                            {
                              if( vbBins[*iIter] )
                              {
                                vbBins[*iIter] = false;
                                ++no_unique;
                              }
                              ++iIter;
                            }
                              
                            REAL_TYPE rErr2 = rErr1*static_cast<REAL_TYPE>(viFoundEdgeIndexes.size())/(SAMPLING_FACTOR*static_cast<REAL_TYPE>(no_unique));
                            if( rErr2 < 0.15 &&  static_cast<REAL_TYPE>(no_unique)/static_cast<REAL_TYPE>(vEdgelets.size()) > 0.3)
                            {//FOUND
                              std::vector<MatchedPair> matchedCorrespondences;
                              for (int i = 0; i < viFullIndexes.size(); i++)
                              {
                                if (viFullIndexes[i] == -1)
                                   continue;
                                MatchedPair mp;
                                mp.v2FstView[0] = vEdgelets[i].v2Pose[0];
                                mp.v2FstView[1] = vEdgelets[i].v2Pose[1];
                                mp.v2SndView[0] = mvEdgelets[viFullIndexes[i]].v2Pose[0];
                                mp.v2SndView[1] = mvEdgelets[viFullIndexes[i]].v2Pose[1];
                                matchedCorrespondences.push_back (mp);
                              }
                              double stdDev;
                              angle = mHomography.estimateRotationAngle (stdDev, matchedCorrespondences);
                              if (stdDev > 0.2)
                                   continue;
                         //     std::cout << "View ID " << cd.iView_ID << " standard dev is " << stdDev << std::endl;
                              double lastRecordedTime = CVD::timer.get_time(); 
                              current_time = CVD::timer.get_time();
                              std::string objName;
                              obj_correct_ID = static_cast<int>(cd.obj_correct_ID);
                              switch( obj_correct_ID )
                              {
                                case 3:
                                  objName = "Plier";
                                  break;
                                case 5:
                                  objName = "Box";
                                  break;
                                case 10:
                                  objName = "Charger";
                                  break;
                                case 4:
                                  objName = "Driver";
                                  break;
                                case 2:
                                  objName = "Hammer";
                                  break;
                                case 1:
                                  objName = "Screw Driver";
                                  break;
                                case 7:
                                  objName = "Stapler";
                                  break;
                                case 6:
                                  objName = "Wood";
                                  break;
                                case 8:
                                  objName = "Wrench";
                                  break;
                                case 9:
                                  objName = "Yellow";
                                  break;
                                default:
                                  objName = "Other";
                              }
                              h = newH*h;
                              obj_detected = true;
                              project_edgelets( vProjectedEdgelets, vEdgelets, h, false );
                              mvDetectedObjects.resize( mvDetectedObjects.size() + 1);
                              mvDetectedObjects.back().iObject_ID = obj_correct_ID;
                              mvDetectedObjects.back().iView_ID = cd.iView_ID + 1;
                              Edges vProjectedEdges;
                              std::vector<TooN::Vector<2> > vv2Pnts;
                              project_edges (vProjectedEdges, vv2Pnts, allEdges, h, false);
                              for (int ipt = 0; ipt < vv2Pnts.size(); ipt++){      
                                  mvDetectedObjects.back().vv2Pnts.push_back (vv2Pnts[ipt]);
                              }
                              // change to the bigger size

                              TooN::Vector<2> y1, y2, y3, y4, y5;
           	              y1[0] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[0]*imgSizeRatio;
                    	      y1[1] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[1]*imgSizeRatio;
	                      y2[0] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
        	              y2[1] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
                	      y3[0] = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
	                      y3[1] = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
        	              y4[0] = mvEdgelets[trdTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
                	      y4[1] = mvEdgelets[trdTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
           	              y5[0] = mvEdgelets[fourtTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
                	      y5[1] = mvEdgelets[fourtTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
                              mvDetectedObjects.back().foundChainPnts.push_back(y1);
                              mvDetectedObjects.back().foundChainPnts.push_back(y2);
                              mvDetectedObjects.back().foundChainPnts.push_back(y3);
                              mvDetectedObjects.back().foundChainPnts.push_back(y4);
                              mvDetectedObjects.back().foundChainPnts.push_back(y5);
                              mvDetectedObjects.back().detection_time = lastRecordedTime - start_time;
                              mvDetectedObjects.back().error = rErr2;

                              mvDetectedObjects.back().pose_x = 0;
                              mvDetectedObjects.back().pose_y = 0;
//                              std::cout << "VIEW ID " << cd.iView_ID << " " << angle << std::endl;
                              mvDetectedObjects.back().pose_z = sin(angle/2);
                              mvDetectedObjects.back().pose_w = cos(angle/2);
                        //      std::cout << "VIEW ID " << cd.iView_ID << " " << angle << std::endl;
                              std::vector<Edgelet>::iterator tedIter = vProjectedEdgelets.begin();
                              CVD::ImageRef irTopLeft(img_h,img_w), irButtomRight(-1,-1);
                              while( tedIter != vProjectedEdgelets.end() )
                              {
                                if( (*tedIter).v2Pose[0] < irTopLeft.x )
                                  irTopLeft.x = (*tedIter).v2Pose[0];
                                else if( (*tedIter).v2Pose[0] > irButtomRight.x )
                                  irButtomRight.x = (*tedIter).v2Pose[0];
                                if( (*tedIter).v2Pose[1] < irTopLeft.y )
                                  irTopLeft.y = (*tedIter).v2Pose[1];
                                else if( (*tedIter).v2Pose[1] > irButtomRight.y )
                                  irButtomRight.y = (*tedIter).v2Pose[1];
                                ++tedIter;
                              }
                              mvDetectedObjects.back().irTopLeft = irTopLeft;
                              mvDetectedObjects.back().irBottomRight = irButtomRight;
                              if( current_time - subtract_time - start_time > MAX_TIME )
                                goto END;
                              // mark detected edges
                              #if 1
                              std::vector<TestDirection>::iterator bTestIter = mvTestDirectionsPreviousObjects[previousObjNo].begin();
                              std::vector<TestDirection>::iterator eTestIter = mvTestDirectionsPreviousObjects[previousObjNo].end();
                              std::vector<bool>::iterator bBinIter = vBinTestDirections.begin();
                              while( bTestIter != eTestIter )
                              {
                                if( !vbBins[bTestIter->iIdx_0] )
                                {
                                  *bBinIter = true;
                                  mvEdgelets[bTestIter->iIdx_0].v2Pose = TooN::Zeros;
                                }
                                else if ( !vbBins[bTestIter->iIdx_1] )
                                {
                                  *bBinIter = true;
                                  mvEdgelets[bTestIter->iIdx_1].v2Pose = TooN::Zeros;
                                }
                                ++bTestIter;
                                ++bBinIter;
                              }
                              #endif

                           //   goto NEXT_ROUND;
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      NEXT_ROUND:;
    }
  END:
  if (!obj_detected && mvDetectedObjects_previousFrame[previousObjNo].missingCount < 5)
  {
     //std::cout << "MISSING COUNT: " << mvDetectedObjects_previousFrame[previousObjNo].missingCount << std::endl;
     mvDetectedObjects_undetected.push_back (mvDetectedObjects_previousFrame[previousObjNo]);
  }
  return obj_detected;
}

int EdgeObjectDetection::detect(int subtract_time)
{
  std::ofstream logFile;
  logFile.open( outputFileName, std::ios::app );
  double current_time;
  int count=0;
  int count_2 = 0;
  int beginIdx = 0;
  int endIdx;
  int numberOfStartingEdges = 0;
  int numberOf2LongChains = 0;
  int numberOf3LongChains = 0;
  int numberOf4LongChains = 0;
  int numberOfCompletedChains = 0;
  double chainDuration = 0;
  int numberOfPairs;
  for (int i = 0; i < NumberOfSteps; i++)
  {
     counterSteps[i] = 0;
     timestampsteps[i] = 0;
  }

  std::vector<MatchedPair> v5MatchedPairs(5);
  for( unsigned int no_class = 0; no_class < mvObjDescriptorClasses.size(); ++no_class )
  {
    TooN::Vector<4> & mv4ChainAngles = mvObjDescriptorClasses[no_class].mv4ChainAngles;
    std::vector<CodeBookElement> & mvCodeBook = mvObjDescriptorClasses[no_class].mvCodeBook;
    std::vector<CVD::ImageRef> & mvirFirstAngleIndexes = mvObjDescriptorClasses[no_class].mvirFirstAngleIndexes;
    std::vector< std::vector<int> > & mvirSecondIndexes = mvObjDescriptorClasses[no_class].mvirSecondIndexes;

#if 0
// This block may help speed up the program if there are more than one class ( fixed angle chains ).
    std::vector<std::pair<REAL_TYPE, int> > mvpairDir_1_Indexes( mvTestDirections.size() ); // having mvpairDir_1_Indexes will help if there are more classes ( more v4ChainAngles ) to be used.
    for( unsigned int no_testdirection = 0; no_testdirection < mvTestDirections.size(); ++ no_testdirection )
      mvpairDir_1_Indexes[no_testdirection] = std::make_pair( mvTestDirections[no_testdirection].rDir_1, no_testdirection );
    std::sort( mvpairDir_1_Indexes.begin(), mvpairDir_1_Indexes.end(), comp_pair_r_i );
    // | dir - mv4ChainAngles[0] | < EPS_DIR_TEST;
    // first part
    REAL_TYPE rUpper_bound = mv4ChainAngles[0] + EPS_DIR_TEST; // first edge link
    REAL_TYPE rLower_bound = mv4ChainAngles[0] - EPS_DIR_TEST;

    std::vector<std::pair<REAL_TYPE, int> >::iterator fstBeginPairIter = std::lower_bound( mvpairDir_1_Indexes.begin(), mvpairDir_1_Indexes.end(), rLower_bound,   comp_lower_pair );
    std::vector<std::pair<REAL_TYPE, int> >::iterator fstEndPairIter = std::upper_bound( fstBeginPairIter, mvpairDir_1_Indexes.end(), rUpper_bound, comp_upper_pair );
    for( ; fstBeginPairIter != fstEndPairIter; ++fstBeginPairIter ) // after passing the first check from mv4ChainAngles[0], 
    {
//      std::cout << fstBeginPairIter->first << " " << fstBeginPairIter->second << std::endl;
    // then check for angle between the pair and elements inside codebooks.
      TestDirection & firstTestDirectionPairs = mvTestDirections[(*fstBeginPairIter).second];
      if( firstTestDirectionPairs.rDist > FIRST_DIST_TRAIN )
        continue;
      count_2 ++;
#else
  #ifdef ECM
  int sum2 = 0, posTD, edgeNo2;

 // double howLong1 = CVD::timer.get_time();
  // here I need to search through the edges only not the pairs
  for (int i = 0; i < mvCellDirections.size(); i++)
  {     
     CellDirection thisCD = mvCellDirections[i];
     if (fabs(thisCD.rDir_1 - mv4ChainAngles[0]) < (EPS_DIR_TEST) || fabs(thisCD.rDir_1 + mv4ChainAngles[0]) < (EPS_DIR_TEST))
     {
        for (int thisCellEdgelets = 0; thisCellEdgelets < ecm[thisCD.cellNo].cellEdgelets.size(); thisCellEdgelets++)
        {
           int edgeNo2 = ecm[thisCD.cellNo].cellEdgelets[thisCellEdgelets];
           if (thisCD.iIdx == edgeNo2)
              continue;
           if (edgeNo2 < thisCD.iIdx)
              posTD = thisCD.iIdx*(mvEdgelets.size()-1) + edgeNo2; 
           else
              posTD = thisCD.iIdx*(mvEdgelets.size()-1) + edgeNo2 - 1;
           vCellTestDirections.push_back (mvTestDirections[posTD]);
           sum2++;
        }
     }
  }
  std::vector<bool> vBinTestDirections(vCellTestDirections.size(),false);
  std::vector<int> vRandomSamples(vCellTestDirections.size());
  for( unsigned int i = 0; i < vCellTestDirections.size(); ++ i ) vRandomSamples[i]=i;  
  endIdx = sum2 - 1;
  numberOfPairs = sum2;
  #else
  std::vector<bool> vBinTestDirections(mvTestDirections.size(),false);
  std::vector<int> vRandomSamples(mvTestDirections.size());
  for( unsigned int i = 0; i < mvTestDirections.size(); ++ i ) vRandomSamples[i]=i;
  endIdx = mvTestDirections.size() - 1;
  numberOfPairs = mvTestDirections.size();
  #endif

  mvCheckedChains.clear();
  double delta_time;
  double step_time = 0;
 // double lastRecordedTime = start_time;

  RESTART:
 // std::cout <<"Shuffle ..." << std::endl;
  std::random_shuffle( vRandomSamples.begin(), vRandomSamples.end() );
  for (int multiScale_level = 0; multiScale_level <= MAX_MULTISCALE_LEVEL; multiScale_level++)
  {
  //std::cerr << "Entering Multi-Scale level = " << multiScale_level << std::endl;
  beginIdx = 0;
  for( ; beginIdx <= endIdx; ++ beginIdx )
  {
    chainDuration += CVD::timer.get_time() - delta_time;
    numberOfStartingEdges++;
    delta_time = CVD::timer.get_time() - start_time; 
  /*  counterSteps[1]++;
    timestampsteps[1] += CVD::timer.get_time() - lastRecordedTime; 
    lastRecordedTime = CVD::timer.get_time();*/
    if( delta_time > MAX_TIME ) goto END;
#ifdef TEST_RESTART    
    if( delta_time > step_time + RESTART_TIME )
    {
      //std::cout << delta_time <<" " << step_time << std::endl;
      step_time += RESTART_TIME;
      goto RESTART;
    }
#endif    
    if( vBinTestDirections[vRandomSamples[beginIdx]] )
      continue;
    #ifdef ECM
    TestDirection & firstTestDirectionPairs = vCellTestDirections[vRandomSamples[beginIdx]];
    #else
    TestDirection & firstTestDirectionPairs = mvTestDirections[vRandomSamples[beginIdx]];
    #endif
    count_2 ++;


    if(!( (fabs( firstTestDirectionPairs.rDir_1-mv4ChainAngles[0] ) < EPS_DIR_TEST || fabs( firstTestDirectionPairs.rDir_1+mv4ChainAngles[0] ) < EPS_DIR_TEST ) && firstTestDirectionPairs.rDist < MAX_DIST_TEST[multiScale_level] &&  firstTestDirectionPairs.rDist > MIN_DIST_TEST[multiScale_level]) )
//    if( !((fabs( firstTestDirectionPairs.rDir_1 - mv4ChainAngles[0] ) < EPS_DIR_TEST) && firstTestDirectionPairs.rDist < FIRST_DIST_TRAIN ))
      continue;
   /* CheckedChains ch;
    ch.chains.reserve (2);
    TooN::Vector<2> v1, v2;
    v1[0] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[0]*imgSizeRatio;
    v1[1] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[1]*imgSizeRatio;
    v2[0] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
    v2[1] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;

    ch.chains.push_back (v1);
    ch.chains.push_back (v2);
    mvCheckedChains.push_back (ch);*/

     // try to get rid of parallel edges?
     if ((firstTestDirectionPairs.rAngle < 0.05 && firstTestDirectionPairs.rAngle > -0.05) ||  (fabs(firstTestDirectionPairs.rAngle - M_PI) < 0.05))
     {
       // std::cerr << "dropping pair with angle " << firstTestDirectionPairs.rAngle << std::endl;
        continue;
     }
#endif
      // looking for a second pair that fst.e2 == snd.e1
      
      int beginSndIndex = mvirDirectionIndexes[firstTestDirectionPairs.iIdx_1].x;
      if( beginSndIndex == -1 ) continue;
      int endSndIndex = mvirDirectionIndexes[firstTestDirectionPairs.iIdx_1].y;
      // searching for a first pair inside codebook.
    /*  counterSteps[2]++;
      timestampsteps[2] += CVD::timer.get_time() - lastRecordedTime; 
      lastRecordedTime = CVD::timer.get_time();*/
      REAL_TYPE rLower_Angle = firstTestDirectionPairs.rAngle - EPS_ANGLE;
      REAL_TYPE rUpper_Angle = firstTestDirectionPairs.rAngle + EPS_ANGLE;
      returnedEdgelets.clear();
      returnedCB.clear();

      for( ; beginSndIndex <= endSndIndex; ++beginSndIndex ) // search through second pair
      {
        delta_time = CVD::timer.get_time() - start_time; 
      /*  counterSteps[3]++;
        timestampsteps[3] += CVD::timer.get_time() - lastRecordedTime; 
        lastRecordedTime = CVD::timer.get_time();*/

        if( delta_time > MAX_TIME ) goto END;
#ifdef TEST_RESTART        
        if( delta_time > step_time + RESTART_TIME )
        {
          //std::cout << delta_time <<" " << step_time << std::endl;
          step_time += RESTART_TIME;
          goto RESTART;
        }
#endif        
        if( vBinTestDirections[beginSndIndex] )
          continue;
        const TestDirection & sndTestDirectionPairs = mvTestDirections[beginSndIndex];
        if( fabs( calculateAngle(sndTestDirectionPairs.rDir_2,firstTestDirectionPairs.rDir_2) - mv4ChainAngles[1] ) < EPS_DIR_TEST && sndTestDirectionPairs.rDist < MAX_DIST_TEST[multiScale_level] &&  sndTestDirectionPairs.rDist > MIN_DIST_TEST[multiScale_level]) // second edge lin k
        {
/*CheckedChains ch;
		ch.chains.reserve (3);
		 TooN::Vector<2> v1, v2, v3;
		 v1[0] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[0]*imgSizeRatio;
		 v1[1] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[1]*imgSizeRatio;
		 v2[0] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
		 v2[1] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
		 v3[0] = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
		 v3[1] = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;

		ch.chains.push_back (v1);
		ch.chains.push_back (v2);
		ch.chains.push_back (v3);
		mvCheckedChains.push_back (ch);*/
        // Check with codebook betweeen lowerCodeBook and upperCodeBook 
          std::vector< std::vector<CodeBookElement>::iterator > vSecondPassCodeBook;
      /*    vSecondPassCodeBook.reserve( upperCodeBook - lowerCodeBook );
          lowerCodeBook = minCodeBook;
          for( ; lowerCodeBook != upperCodeBook; ++lowerCodeBook )
          {
            thisCB.push_back(lowerCodeBook);
            if( fabs( (*lowerCodeBook).rAngle_1 - sndTestDirectionPairs.rAngle ) < EPS_ANGLE 
               && fabs( (*lowerCodeBook).rRelativeDist_1 - sndTestDirectionPairs.rDist/firstTestDirectionPairs.rDist ) < EPS_DIST )
              vSecondPassCodeBook.push_back( lowerCodeBook );
          }*/
          // trying to read from the codebook

          int indexNo = 64.f*(firstTestDirectionPairs.rAngle/M_PI + 0.5f); if( indexNo < 0 ) indexNo = 0;
          int indexNo2 = 64.f*(sndTestDirectionPairs.rAngle/M_PI + 0.5f); if( indexNo2 < 0 ) indexNo2 = 0;
          int indexNo3;
          REAL_TYPE ratio = sndTestDirectionPairs.rDist/firstTestDirectionPairs.rDist;
          if (ratio < 1)
              indexNo3 = int(ratio/0.1);
          else if (ratio > 1 && ratio < 29.5)
              indexNo3 = int (ratio)+9;
          else
              indexNo3 = 38;
          std::vector<int> & thisList = mvirSecondIndexes[(indexNo)*2496+(indexNo2)*39+indexNo3];
          vSecondPassCodeBook.reserve( thisList.size() );
         // int distBin [25] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
          for(int k = 0; k < thisList.size(); k++)
          {
            std::vector<CodeBookElement>::iterator lx = mvCodeBook.begin()+thisList[k];
            //thisCB.push_back(lx);
            if( fabs( (*lx).rAngle_1 - sndTestDirectionPairs.rAngle ) < EPS_ANGLE 
               && fabs( (*lx).rRelativeDist_1 - sndTestDirectionPairs.rDist/firstTestDirectionPairs.rDist ) < EPS_DIST )
             {
                 vSecondPassCodeBook.push_back( lx );
                 REAL_TYPE lookingforDist = (*lx).rRelativeDist_2 * sndTestDirectionPairs.rDist;
              //   int binNo = (int)(lookingforDist/10);
              //   distBin[binNo]++;
             }
          }   

          if( vSecondPassCodeBook.size() == 0 )
          {
             continue;
          }
          //numberOf3LongChains++;
          int beginTrdIndex = mvirDirectionIndexes[sndTestDirectionPairs.iIdx_1].x;
          if( beginTrdIndex == -1 ) continue;
          int endTrdIndex = mvirDirectionIndexes[sndTestDirectionPairs.iIdx_1].y;
          
          for( ; beginTrdIndex <= endTrdIndex; ++beginTrdIndex ) // search through third pair
          {
            delta_time = CVD::timer.get_time() - start_time;
          /*  counterSteps[5]++;
           timestampsteps[5] += CVD::timer.get_time() - lastRecordedTime; 
           lastRecordedTime = CVD::timer.get_time();*/

            if( delta_time > MAX_TIME ) goto END;
#ifdef TEST_RESTART            
            if( delta_time > step_time + RESTART_TIME )
            {
              step_time += RESTART_TIME;
              goto RESTART;
            }
#endif            
            if( vBinTestDirections[beginTrdIndex] )
              continue;
            TestDirection & trdTestDirectionPairs = mvTestDirections[beginTrdIndex];
            if( fabs( calculateAngle(trdTestDirectionPairs.rDir_2, sndTestDirectionPairs.rDir_2) - mv4ChainAngles[2] ) < EPS_DIR_TEST && trdTestDirectionPairs.rDist < MAX_DIST_TEST[multiScale_level] &&  trdTestDirectionPairs.rDist > MIN_DIST_TEST[multiScale_level]) // third edge link
            {
/*              int binNo = (int)(trdTestDirectionPairs.rDist/10);
              int countBins = distBin[binNo];
              if (binNo > 0)
                  binNo += distBin[binNo-1];
              if (binNo < 24)
                  binNo += distBin[binNo+1];*/
              
            //  std::cerr << ((binNo > 0) ? distBin[binNo-1] : "") << " " << distBin[binNo] << " " << ((binNo < 24) ? distBin[binNo+1] : "") << std::endl;
              std::vector< std::vector<CodeBookElement>::iterator > vThirdPassCodeBook;
              vThirdPassCodeBook.reserve( vSecondPassCodeBook.size() );
              for( unsigned int snd_pass = 0; snd_pass < vSecondPassCodeBook.size(); ++snd_pass )
              {
                if( fabs( vSecondPassCodeBook[snd_pass]->rAngle_2 - trdTestDirectionPairs.rAngle ) < EPS_ANGLE
                   && fabs( vSecondPassCodeBook[snd_pass]->rRelativeDist_2 - trdTestDirectionPairs.rDist/sndTestDirectionPairs.rDist ) < EPS_DIST)
                  vThirdPassCodeBook.push_back( vSecondPassCodeBook[snd_pass] );
                  // here we can identify which distances and relative orientations are we looking at?
                  //std::cerr << "Should be looking for distance: " << vSecondPassCodeBook[snd_pass]->rRelativeDist_2 << " " << sndTestDirectionPairs.rDist << " " << std::endl;
              }
              if( vThirdPassCodeBook.size() == 0 ) continue;

              //numberOf4LongChains++;
              int beginFrtIndex = mvirDirectionIndexes[trdTestDirectionPairs.iIdx_1].x;
              if( beginFrtIndex == -1 ) continue;
              int endFrtIndex = mvirDirectionIndexes[trdTestDirectionPairs.iIdx_1].y;

              for( ; beginFrtIndex <= endFrtIndex; ++beginFrtIndex ) // search through fourth pair
              {
                delta_time = CVD::timer.get_time() - start_time;
               /* counterSteps[7]++;
                timestampsteps[7] += CVD::timer.get_time() - lastRecordedTime; 
                lastRecordedTime = CVD::timer.get_time();*/
 
                if( delta_time > MAX_TIME ) goto END;
#ifdef TEST_RESTART                
                if( delta_time > step_time + RESTART_TIME )
                {
                  //std::cout << delta_time <<" " << step_time << std::endl;
                  step_time += RESTART_TIME;
                  goto RESTART;
                }
#endif                
                if( vBinTestDirections[beginFrtIndex] )
                  continue;
                TestDirection & fourtTestDirectionPairs = mvTestDirections[beginFrtIndex];
                if( fabs( calculateAngle(fourtTestDirectionPairs.rDir_2,trdTestDirectionPairs.rDir_2) - mv4ChainAngles[3] ) < EPS_DIR_TEST && fourtTestDirectionPairs.rDist < MAX_DIST_TEST[multiScale_level] &&  fourtTestDirectionPairs.rDist > MIN_DIST_TEST[multiScale_level]) // fourth edge link
                {

                  //numberOfCompletedChains++;
                  // Check with codebook insdie vThirdPassCodeBook;

                  /*mvCandidateObjects.clear();
                  // this was added to check how many chains actually match BEFORE the distance transform error is calculated
                  int counterMatches = 0;
                  int counterHomMatches = 0;
                  for( unsigned int trd_pass = 0; trd_pass < vThirdPassCodeBook.size(); ++trd_pass )
                  {
                    if( fabs( vThirdPassCodeBook[trd_pass]->rAngle_3 - fourtTestDirectionPairs.rAngle ) < EPS_ANGLE 
                       && fabs( vThirdPassCodeBook[trd_pass]->rRelativeDist_3 - fourtTestDirectionPairs.rDist/trdTestDirectionPairs.rDist ) < EPS_DIST )
                    {
                       counterMatches++; 
                       CodeBookElement & cd = *(vThirdPassCodeBook[trd_pass]);
                       View & view = mvObjDescriptorClasses[no_class].mObjectsTemplate.vObjectElements[ cd.iObject_ID ].vViews[ cd.iView_ID ];
                       std::vector<Edgelet> & vEdgelets = view.vEdgelets;
                       Edges & allEdges = view.allEdges;
                       v5MatchedPairs[0].v2FstView = vEdgelets[cd.iIdx_0].v2Pose; v5MatchedPairs[0].v2SndView = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose;
                       v5MatchedPairs[1].v2FstView = vEdgelets[cd.iIdx_1].v2Pose; v5MatchedPairs[1].v2SndView = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose;
                       v5MatchedPairs[2].v2FstView = vEdgelets[cd.iIdx_2].v2Pose; v5MatchedPairs[2].v2SndView = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose;
                       v5MatchedPairs[3].v2FstView = vEdgelets[cd.iIdx_3].v2Pose; v5MatchedPairs[3].v2SndView = mvEdgelets[trdTestDirectionPairs.iIdx_1].v2Pose;
                       v5MatchedPairs[4].v2FstView = vEdgelets[cd.iIdx_4].v2Pose; v5MatchedPairs[4].v2SndView = mvEdgelets[fourtTestDirectionPairs.iIdx_1].v2Pose;
                       TooN::Matrix<3> h;
                       if ( mHomography.estimateHomography_group_of_5(h, v5MatchedPairs) )
                       {
                         counterHomMatches++;
                         std::vector<Edgelet> vProjectedEdgelets;
                         project_edgelets( vProjectedEdgelets, vEdgelets, h );
                         std::vector<TooN::Vector<2> > vv2Pnts;
                         Edges vProjectedEdges;
                         project_edges (vProjectedEdges, vv2Pnts, allEdges, h, false);
                         mvCandidateObjects.resize( mvCandidateObjects.size() + 1);
                         int obj_correct_ID = static_cast<int>(cd.obj_correct_ID);
                         mvCandidateObjects.back().iObject_ID = obj_correct_ID;

                         for (int ipt = 0; ipt < vv2Pnts.size(); ipt++){      
                            mvCandidateObjects.back().vv2Pnts.push_back (vv2Pnts[ipt]);
                         } 
                       }

                    }
                  } 
                  if (counterMatches > 0) 
                     std::cout << "number of matches with codebook " << counterMatches << " with homographies " << counterHomMatches << std::endl;                
                  */
                //  std::cerr << "Potential codebook matches " << vThirdPassCodeBook.size() << std::endl;
                  for( unsigned int trd_pass = 0; trd_pass < vThirdPassCodeBook.size(); ++trd_pass )
                  {
                    delta_time = CVD::timer.get_time() - start_time; 
                    if( delta_time > MAX_TIME ) goto END;
#ifdef TEST_RESTART                    
                    if( delta_time > step_time + RESTART_TIME )
                    {
                      //std::cout << delta_time <<" " << step_time << std::endl;
                      step_time += RESTART_TIME;
                      goto RESTART;
                    }
#endif                 
                    if( fabs( vThirdPassCodeBook[trd_pass]->rAngle_3 - fourtTestDirectionPairs.rAngle ) < EPS_ANGLE 
                       && fabs( vThirdPassCodeBook[trd_pass]->rRelativeDist_3 - fourtTestDirectionPairs.rDist/trdTestDirectionPairs.rDist ) < EPS_DIST )
                    { 
                     // std::cerr << "passed check in codebook for vThirdPassCodeBook" << std::endl;
                      CodeBookElement & cd = *(vThirdPassCodeBook[trd_pass]);
                      View & view = mvObjDescriptorClasses[no_class].mObjectsTemplate.vObjectElements[ cd.iObject_ID ].vViews[ cd.iView_ID ];
                      std::vector<Edgelet> & vEdgelets = view.vEdgelets;
                      Edges & allEdges = view.allEdges;
                      if (cd.iIdx_0 >= vEdgelets.size() || cd.iIdx_1 >= vEdgelets.size() || cd.iIdx_2 >= vEdgelets.size() || cd.iIdx_3 >= vEdgelets.size() || cd.iIdx_4 >= vEdgelets.size())
                      {
                          continue;
                      }
                      v5MatchedPairs[0].v2FstView = vEdgelets.at(cd.iIdx_0).v2Pose; v5MatchedPairs[0].v2SndView = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose;
                      v5MatchedPairs[1].v2FstView = vEdgelets.at(cd.iIdx_1).v2Pose; v5MatchedPairs[1].v2SndView = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose;
                      v5MatchedPairs[2].v2FstView = vEdgelets.at(cd.iIdx_2).v2Pose; v5MatchedPairs[2].v2SndView = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose;
                      v5MatchedPairs[3].v2FstView = vEdgelets.at(cd.iIdx_3).v2Pose; v5MatchedPairs[3].v2SndView = mvEdgelets[trdTestDirectionPairs.iIdx_1].v2Pose;
                      v5MatchedPairs[4].v2FstView = vEdgelets.at(cd.iIdx_4).v2Pose; v5MatchedPairs[4].v2SndView = mvEdgelets[fourtTestDirectionPairs.iIdx_1].v2Pose;
                      // add the chain
                      TooN::Vector<2> v1, v2, v3, v4, v5;
                      v1[0] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[0]*imgSizeRatio;
                      v1[1] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[1]*imgSizeRatio;
                      v2[0] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
                      v2[1] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
                      v3[0] = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
                      v3[1] = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
                      v4[0] = mvEdgelets[trdTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
                      v4[1] = mvEdgelets[trdTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
                      v5[0] = mvEdgelets[fourtTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
                      v5[1] = mvEdgelets[fourtTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
                      CheckedChains ch;
                      ch.chains.reserve (5);
                      ch.chains.push_back (v1);
                      ch.chains.push_back (v2);
                      ch.chains.push_back (v3);
                      ch.chains.push_back (v4);
                      ch.chains.push_back (v5);
                      mvCheckedChains.push_back (ch);
                      TooN::Matrix<3> h;
                      if ( mHomography.estimateHomography_group_of_5(h, v5MatchedPairs) )
                      {
                        std::vector<Edgelet> vProjectedEdgelets;
                        project_edgelets( vProjectedEdgelets, vEdgelets, h );
                        std::vector<REAL_TYPE> vrDistanceSum( vProjectedEdgelets.size() );

                        if (vProjectedEdgelets.size() == 0)
                           continue;
                        for( unsigned int no_edgelet = 0; no_edgelet < vProjectedEdgelets.size(); ++ no_edgelet )
                        {
                            if (isnan(vProjectedEdgelets[no_edgelet].rSlope))
                              continue;
                            int binNo = static_cast<int>( (atan( vProjectedEdgelets[no_edgelet].rSlope )/M_PI + 0.5f ) * 10.f );
                            vrDistanceSum[no_edgelet] = mvEdgeOrientationDistTransforms[binNo][vProjectedEdgelets[no_edgelet].v2Pose];
                        }
                        std::sort( vrDistanceSum.begin(), vrDistanceSum.end() );
                        REAL_TYPE rMedian = vrDistanceSum[vrDistanceSum.size()/2];
                       // std::cerr << "Before checking MAX_DT " << rMedian << " " << MAX_DT << std::endl;
                        if( rMedian < MAX_DT )
                        {
                          project_edgelets( vProjectedEdgelets, vEdgelets, h, false );

                          TooN::Matrix<3> newH;
                          std::vector<int> viFoundEdgeIndexes;
                          std::vector<int> viFullIndexes;
                          std::vector<MatchedPair> vMatchedPairs;
                          double angle;
                          REAL_TYPE rErr1 = iterativeClosestEdges2_NP( newH, viFullIndexes, viFoundEdgeIndexes, vProjectedEdgelets, mvEdgelets, vMatchedPairs, 25.f, 30.f, 0.2f);
                          std::vector<MatchedPair> matchedCorrespondences;
                          for (int i = 0; i < viFullIndexes.size(); i++)
                          {
                             if (viFullIndexes[i] == -1)
                                continue;
                             MatchedPair mp;
                             mp.v2FstView[0] = vEdgelets[i].v2Pose[0];
                             mp.v2FstView[1] = vEdgelets[i].v2Pose[1];
                             mp.v2SndView[0] = mvEdgelets[viFullIndexes[i]].v2Pose[0];
                             mp.v2SndView[1] = mvEdgelets[viFullIndexes[i]].v2Pose[1];
                             matchedCorrespondences.push_back (mp);
                          }
                          double stdDev;
                          angle = mHomography.estimateRotationAngle (stdDev, matchedCorrespondences);
                          //std::cerr << "Standard Deviation of Rotation Angle " << stdDev << " " << viFoundEdgeIndexes.size() << std::endl;
                          if (stdDev > 0.2)
                              continue;
                          std::vector<bool> vbBins( mvEdgelets.size(), true );
                          if( viFoundEdgeIndexes.size() > 5 )
                          {
                            int no_unique = 0;
                            std::vector<int>::iterator iIter = viFoundEdgeIndexes.begin();
                            std::vector<int>::iterator iEndIter = viFoundEdgeIndexes.end();
                            while( iIter != iEndIter )
                            {
                              if( vbBins[*iIter] )
                              {
                                vbBins[*iIter] = false;
                                ++no_unique;
                              }
                              ++iIter;
                            }
                              
                            REAL_TYPE rErr2 = rErr1*static_cast<REAL_TYPE>(viFoundEdgeIndexes.size())/(SAMPLING_FACTOR*static_cast<REAL_TYPE>(vMatchedPairs.size()));
                  //          std::cerr << "reached here with " << cd.obj_correct_ID << " " << rErr1 << " " << rErr2 << " " << (static_cast<REAL_TYPE>(vMatchedPairs.size())/static_cast<REAL_TYPE>(vEdgelets.size())) << std::endl;
                            
                            if( rErr1 < 0.07 && rErr2 < 0.04 &&  static_cast<REAL_TYPE>(vMatchedPairs.size())/static_cast<REAL_TYPE>(vEdgelets.size()) > 0.55)
                            {//FOUND
                            //  
                              double lastRecordedTime = CVD::timer.get_time(); 
                              current_time = CVD::timer.get_time();
                              std::string objName;
                              obj_correct_ID = static_cast<int>(cd.obj_correct_ID);
                              switch( obj_correct_ID )
                              {
                                case 3:
                                  objName = "Plier";
                                  break;
                                case 5:
                                  objName = "Box";
                                  break;
                                case 10:
                                  objName = "Charger";
                                  break;
                                case 4:
                                  objName = "Driver";
                                  break;
                                case 2:
                                  objName = "Hammer";
                                  break;
                                case 1:
                                  objName = "Screw Driver";
                                  break;
                                case 7:
                                  objName = "Stapler";
                                  break;
                                case 6:
                                  objName = "Wood";
                                  break;
                                case 8:
                                  objName = "Wrench";
                                  break;
                                case 9:
                                  objName = "Yellow";
                                  break;
                                default:
                                  objName = "Other";
                              }
                              h = newH*h;
                              project_edgelets( vProjectedEdgelets, vEdgelets, h, false );
                              mvDetectedObjects.resize( mvDetectedObjects.size() + 1);
                              mvDetectedObjects.back().iObject_ID = obj_correct_ID;
                              mvDetectedObjects.back().iView_ID = cd.iView_ID + 1;
                              // viFoundEdgeIndexes
                              std::vector<int> posRegionNo;
                              std::vector<int> posRegionNoCount;
                              bool found = false;
                              for (int vi = 0; vi < viFoundEdgeIndexes.size(); vi++)
                              {
                                 found = false;
                                 for (int ri = 0; ri < posRegionNo.size(); ri++)
                                    if (posRegionNo[ri] == mvEdgelets[viFoundEdgeIndexes[vi]].regionNo)
                                    {
                                        posRegionNoCount[ri]++;
                                        found = true;
                                        break;
                                    }
                                 if (!found)
                                 {
                                    posRegionNo.push_back(mvEdgelets[viFoundEdgeIndexes[vi]].regionNo);
                                    posRegionNoCount.push_back(1);
                                 } 
                              }
                              int max = 0;
                              int chosenRegionNo = -1;
                              for (int ri = 0; ri < posRegionNo.size(); ri++)
                                if (posRegionNoCount[ri] > max)
                                {
                                    max = posRegionNoCount[ri];
                                    chosenRegionNo = ri;
                                }
                              mvDetectedObjects.back().iRegion_No = chosenRegionNo;
                              Edges vProjectedEdges;
                              std::vector<TooN::Vector<2> > vv2Pnts;
                              project_edges (vProjectedEdges, vv2Pnts, allEdges, h, false);
                              for (int ipt = 0; ipt < vv2Pnts.size(); ipt++){      
                                  mvDetectedObjects.back().vv2Pnts.push_back (vv2Pnts[ipt]);
                              }
                              // change to the bigger size

                              TooN::Vector<2> y1, y2, y3, y4, y5;
           	              y1[0] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[0]*imgSizeRatio;
                    	      y1[1] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[1]*imgSizeRatio;
	                      y2[0] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
        	              y2[1] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
                	      y3[0] = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
	                      y3[1] = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
        	              y4[0] = mvEdgelets[trdTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
                	      y4[1] = mvEdgelets[trdTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
           	              y5[0] = mvEdgelets[fourtTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
                	      y5[1] = mvEdgelets[fourtTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
                              mvDetectedObjects.back().foundChainPnts.push_back(y1);
                              mvDetectedObjects.back().foundChainPnts.push_back(y2);
                              mvDetectedObjects.back().foundChainPnts.push_back(y3);
                              mvDetectedObjects.back().foundChainPnts.push_back(y4);
                              mvDetectedObjects.back().foundChainPnts.push_back(y5);
                              mvDetectedObjects.back().detection_time = lastRecordedTime - start_time;
                              std::cout << "Detection time : " << mvDetectedObjects.back().detection_time << " " << no_class << std::endl;
                              mvDetectedObjects.back().error = rErr2;

                              // homography
                              /*double trace = 1 + h(0,0) + h(1,1) + h(2,2);
                              if (trace > 1.000001)
                              {
                                 double s = 1/(sqrt (trace) * 2);
                                 mvDetectedObjects.back().pose_x = (h(1,2) - h(2,1)) / s;
                                 mvDetectedObjects.back().pose_y = (h(2,0) - h(0,2)) / s;
                                 mvDetectedObjects.back().pose_z = (h(1,0) - h(1,0)) / s;
                                 mvDetectedObjects.back().pose_w = trace * s;
                              }
                              else
                              {
                                 std::cout << "below trace = " << trace << std::endl;
                              }*/
                              mvDetectedObjects.back().pose_x = 0;
                              mvDetectedObjects.back().pose_y = 0;
                              mvDetectedObjects.back().pose_z = sin(angle/2);
                              mvDetectedObjects.back().pose_w = cos(angle/2);

                              // record the object ID in case a later view is added
                              latestObjID = cd.iObject_ID;
                              latest_correctObjID = obj_correct_ID;
                              
                              std::vector<Edgelet>::iterator tedIter = vProjectedEdgelets.begin();
                              CVD::ImageRef irTopLeft(img_h,img_w), irButtomRight(-1,-1);
                              while( tedIter != vProjectedEdgelets.end() )
                              {
                                if( (*tedIter).v2Pose[0] < irTopLeft.x )
                                  irTopLeft.x = (*tedIter).v2Pose[0];
                                else if( (*tedIter).v2Pose[0] > irButtomRight.x )
                                  irButtomRight.x = (*tedIter).v2Pose[0];
                                if( (*tedIter).v2Pose[1] < irTopLeft.y )
                                  irTopLeft.y = (*tedIter).v2Pose[1];
                                else if( (*tedIter).v2Pose[1] > irButtomRight.y )
                                  irButtomRight.y = (*tedIter).v2Pose[1];
                                ++tedIter;
                              }
                              mvDetectedObjects.back().irTopLeft = irTopLeft;
                              mvDetectedObjects.back().irBottomRight = irButtomRight;
                              count++;
                              if( current_time - subtract_time - start_time > MAX_TIME )
                                goto END;
                              // mark detected edges
                              #if 1
                              std::vector<TestDirection>::iterator bTestIter = mvTestDirections.begin();
                              std::vector<TestDirection>::iterator eTestIter = mvTestDirections.end();
                              std::vector<bool>::iterator bBinIter = vBinTestDirections.begin();
                              while( bTestIter != eTestIter )
                              {
                                if( !vbBins[bTestIter->iIdx_0] )
                                {
                                  *bBinIter = true;
                                  mvEdgelets[bTestIter->iIdx_0].v2Pose = TooN::Zeros;
                                }
                                else if ( !vbBins[bTestIter->iIdx_1] )
                                {
                                  *bBinIter = true;
                                  mvEdgelets[bTestIter->iIdx_1].v2Pose = TooN::Zeros;
                                }
                                ++bTestIter;
                                ++bBinIter;
                              }
                              #endif

                           //   goto NEXT_ROUND;
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      NEXT_ROUND:;
    }
    }
  }
  END:

  // time analysis stages

/*  for (int i = 1; i < NumberOfSteps; i++)
  {
     std::cout << "Counter of step " << i << " is " << counterSteps[i] << ", average = " << timestampsteps[i]/counterSteps[i] << ", total = " << timestampsteps[i] << std::endl;
  }*/

//  std::cout << "number of pairs tested " << numberOfStartingEdges << " out of " << numberOfPairs << std::endl;
//  std::cout << "number of chains of length 2 " << numberOf2LongChains << std::endl;
//  std::cout << "number of chains of length 3 " << numberOf3LongChains << std::endl;
//  std::cout << "number of chains of length 4 " << numberOf4LongChains << std::endl;
//  std::cout << "number of full length chains " << numberOfCompletedChains << std::endl;
//  std::cout << "average time for chains " << chainDuration/numberOfPairs << std::endl;
  
  logFile.close();
 // std::cout << "exiting with ********** " << mvDetectedObjects.size() << " and checked chains " << mvCheckedChains.size () << std::endl;
 // for (int i = 0; i < mvDetectedObjects.size(); i++)
 //    std::cout << "Object " << i << " took " << mvDetectedObjects[i].detection_time << " seconds" << std::endl;
  return ( mvDetectedObjects.size() != 0) ? true : false;
}

int EdgeObjectDetection::detect_star(int subtract_time)
{
  std::ofstream logFile;
  logFile.open( outputFileName, std::ios::app );
  double current_time;
  int count=0;
  int count_2 = 0;
  int beginIdx = 0;
  int endIdx;
  int numberOfStartingEdges = 0;
  int numberOf2LongChains = 0;
  int numberOf3LongChains = 0;
  int numberOf4LongChains = 0;
  int numberOfCompletedChains = 0;
  double chainDuration = 0;
  int numberOfPairs;
  for (int i = 0; i < NumberOfSteps; i++)
  {
     counterSteps[i] = 0;
     timestampsteps[i] = 0;
  }

  std::vector<MatchedPair> v5MatchedPairs(5);
  for( unsigned int no_class = 0; no_class < mvObjDescriptorClasses.size(); ++no_class )
  {
    TooN::Vector<4> & mv4ChainAngles = mvObjDescriptorClasses[no_class].mv4ChainAngles;
    std::vector<CodeBookElement> & mvCodeBook = mvObjDescriptorClasses[no_class].mvCodeBook;
    std::vector<CVD::ImageRef> & mvirFirstAngleIndexes = mvObjDescriptorClasses[no_class].mvirFirstAngleIndexes;
    std::vector< std::vector<int> > & mvirSecondIndexes = mvObjDescriptorClasses[no_class].mvirSecondIndexes;

#if 0
#else
  std::vector<bool> vBinTestDirections(mvTestDirections.size(),false);
  std::vector<int> vRandomSamples(mvTestDirections.size());
  for( unsigned int i = 0; i < mvTestDirections.size(); ++ i ) vRandomSamples[i]=i;
  endIdx = mvTestDirections.size() - 1;
  numberOfPairs = mvTestDirections.size();

  mvCheckedChains.clear();
  double delta_time;
  double step_time = 0;

  RESTART:
 // std::cout <<"Shuffle ..." << std::endl;
  std::random_shuffle( vRandomSamples.begin(), vRandomSamples.end() );
  for (int multiScale_level = 0; multiScale_level <= MAX_MULTISCALE_LEVEL; multiScale_level++)
  {
  //std::cerr << "Entering Multi-Scale level = " << multiScale_level << std::endl;
  beginIdx = 0;
  for( ; beginIdx <= endIdx; ++ beginIdx )
  {
    chainDuration += CVD::timer.get_time() - delta_time;
    numberOfStartingEdges++;
    delta_time = CVD::timer.get_time() - start_time; 
    if( delta_time > MAX_TIME ) goto END;
#ifdef TEST_RESTART    
    if( delta_time > step_time + RESTART_TIME )
    {
      //std::cout << delta_time <<" " << step_time << std::endl;
      step_time += RESTART_TIME;
      goto RESTART;
    }
#endif    
    if( vBinTestDirections[vRandomSamples[beginIdx]] )
      continue;
    TestDirection & firstTestDirectionPairs = mvTestDirections[vRandomSamples[beginIdx]];
    count_2 ++;


    if(!( (fabs( firstTestDirectionPairs.rDir_1-mv4ChainAngles[0] ) < EPS_DIR_TEST || fabs( firstTestDirectionPairs.rDir_1+mv4ChainAngles[0] ) < EPS_DIR_TEST ) && firstTestDirectionPairs.rDist < MAX_DIST_TEST[multiScale_level] &&  firstTestDirectionPairs.rDist > MIN_DIST_TEST[multiScale_level]) )
      continue;
    CheckedChains ch;
    ch.chains.reserve (2);
    TooN::Vector<2> v1, v2;
    v1[0] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[0]*imgSizeRatio;
    v1[1] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[1]*imgSizeRatio;
    v2[0] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
    v2[1] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;

    ch.chains.push_back (v1);
    ch.chains.push_back (v2);
    mvCheckedChains.push_back (ch);

     // try to get rid of parallel edges?
     if ((firstTestDirectionPairs.rAngle < 0.05 && firstTestDirectionPairs.rAngle > -0.05) ||  (fabs(firstTestDirectionPairs.rAngle - M_PI) < 0.05))
     {
       // std::cerr << "dropping pair with angle " << firstTestDirectionPairs.rAngle << std::endl;
        continue;
     }
#endif
      // looking for a second pair that fst.e2 == snd.e1
      
      int beginSndIndex = mvirDirectionIndexes[firstTestDirectionPairs.iIdx_0].x;
      if( beginSndIndex == -1 ) continue;
      int endSndIndex = mvirDirectionIndexes[firstTestDirectionPairs.iIdx_0].y;
      // searching for a first pair inside codebook.
      REAL_TYPE rLower_Angle = firstTestDirectionPairs.rAngle - EPS_ANGLE;
      REAL_TYPE rUpper_Angle = firstTestDirectionPairs.rAngle + EPS_ANGLE;
      returnedEdgelets.clear();
      returnedCB.clear();

      for( ; beginSndIndex <= endSndIndex; ++beginSndIndex ) // search through second pair
      {
        delta_time = CVD::timer.get_time() - start_time; 

        if( delta_time > MAX_TIME ) goto END;
#ifdef TEST_RESTART        
        if( delta_time > step_time + RESTART_TIME )
        {
          //std::cout << delta_time <<" " << step_time << std::endl;
          step_time += RESTART_TIME;
          goto RESTART;
        }
#endif        
        if( vBinTestDirections[beginSndIndex] )
          continue;
        const TestDirection & sndTestDirectionPairs = mvTestDirections[beginSndIndex];
        if( fabs( calculateAngle(sndTestDirectionPairs.rDir_2,firstTestDirectionPairs.rDir_2) - mv4ChainAngles[1] ) < EPS_DIR_TEST && sndTestDirectionPairs.rDist < MAX_DIST_TEST[multiScale_level] &&  sndTestDirectionPairs.rDist > MIN_DIST_TEST[multiScale_level]) // second edge lin k
        {
                /*CheckedChains ch;
		ch.chains.reserve (3);
		 TooN::Vector<2> v1, v2, v3;
		 v1[0] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[0]*imgSizeRatio;
		 v1[1] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[1]*imgSizeRatio;
		 v2[0] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
		 v2[1] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
		 v3[0] = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
		 v3[1] = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;

		ch.chains.push_back (v1);
		ch.chains.push_back (v2);
		ch.chains.push_back (v3);
		mvCheckedChains.push_back (ch);*/
        // Check with codebook betweeen lowerCodeBook and upperCodeBook 
          std::vector< std::vector<CodeBookElement>::iterator > vSecondPassCodeBook;
      /*    vSecondPassCodeBook.reserve( upperCodeBook - lowerCodeBook );
          lowerCodeBook = minCodeBook;
          for( ; lowerCodeBook != upperCodeBook; ++lowerCodeBook )
          {
            thisCB.push_back(lowerCodeBook);
            if( fabs( (*lowerCodeBook).rAngle_1 - sndTestDirectionPairs.rAngle ) < EPS_ANGLE 
               && fabs( (*lowerCodeBook).rRelativeDist_1 - sndTestDirectionPairs.rDist/firstTestDirectionPairs.rDist ) < EPS_DIST )
              vSecondPassCodeBook.push_back( lowerCodeBook );
          }*/
          // trying to read from the codebook

          int indexNo = 64.f*(firstTestDirectionPairs.rAngle/M_PI + 0.5f); if( indexNo < 0 ) indexNo = 0;
          int indexNo2 = 64.f*(sndTestDirectionPairs.rAngle/M_PI + 0.5f); if( indexNo2 < 0 ) indexNo2 = 0;
          int indexNo3;
          REAL_TYPE ratio = sndTestDirectionPairs.rDist/firstTestDirectionPairs.rDist;
          if (ratio < 1)
              indexNo3 = int(ratio/0.1);
          else if (ratio > 1 && ratio < 29.5)
              indexNo3 = int (ratio)+9;
          else
              indexNo3 = 38;
          std::vector<int> & thisList = mvirSecondIndexes[(indexNo)*2496+(indexNo2)*39+indexNo3];
          vSecondPassCodeBook.reserve( thisList.size() );
         // int distBin [25] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
          for(int k = 0; k < thisList.size(); k++)
          {
            std::vector<CodeBookElement>::iterator lx = mvCodeBook.begin()+thisList[k];
            //thisCB.push_back(lx);
            if( fabs( (*lx).rAngle_1 - sndTestDirectionPairs.rAngle ) < EPS_ANGLE 
               && fabs( (*lx).rRelativeDist_1 - sndTestDirectionPairs.rDist/firstTestDirectionPairs.rDist ) < EPS_DIST )
             {
                 vSecondPassCodeBook.push_back( lx );
                 REAL_TYPE lookingforDist = (*lx).rRelativeDist_2 * sndTestDirectionPairs.rDist;
              //   int binNo = (int)(lookingforDist/10);
              //   distBin[binNo]++;
             }
          }   

          if( vSecondPassCodeBook.size() == 0 )
          {
             continue;
          }
          //numberOf3LongChains++;
          int beginTrdIndex = mvirDirectionIndexes[sndTestDirectionPairs.iIdx_0].x;
          if( beginTrdIndex == -1 ) continue;
          int endTrdIndex = mvirDirectionIndexes[sndTestDirectionPairs.iIdx_0].y;
          
          for( ; beginTrdIndex <= endTrdIndex; ++beginTrdIndex ) // search through third pair
          {
            delta_time = CVD::timer.get_time() - start_time;
          /*  counterSteps[5]++;
           timestampsteps[5] += CVD::timer.get_time() - lastRecordedTime; 
           lastRecordedTime = CVD::timer.get_time();*/

            if( delta_time > MAX_TIME ) goto END;
#ifdef TEST_RESTART            
            if( delta_time > step_time + RESTART_TIME )
            {
              step_time += RESTART_TIME;
              goto RESTART;
            }
#endif            
            if( vBinTestDirections[beginTrdIndex] )
              continue;
            TestDirection & trdTestDirectionPairs = mvTestDirections[beginTrdIndex];
            if( fabs( calculateAngle(trdTestDirectionPairs.rDir_2, sndTestDirectionPairs.rDir_2) - mv4ChainAngles[2] ) < EPS_DIR_TEST && trdTestDirectionPairs.rDist < MAX_DIST_TEST[multiScale_level] &&  trdTestDirectionPairs.rDist > MIN_DIST_TEST[multiScale_level]) // third edge link
            {
/*              int binNo = (int)(trdTestDirectionPairs.rDist/10);
              int countBins = distBin[binNo];
              if (binNo > 0)
                  binNo += distBin[binNo-1];
              if (binNo < 24)
                  binNo += distBin[binNo+1];*/
              
            //  std::cerr << ((binNo > 0) ? distBin[binNo-1] : "") << " " << distBin[binNo] << " " << ((binNo < 24) ? distBin[binNo+1] : "") << std::endl;
              std::vector< std::vector<CodeBookElement>::iterator > vThirdPassCodeBook;
              vThirdPassCodeBook.reserve( vSecondPassCodeBook.size() );
              for( unsigned int snd_pass = 0; snd_pass < vSecondPassCodeBook.size(); ++snd_pass )
              {
                if( fabs( vSecondPassCodeBook[snd_pass]->rAngle_2 - trdTestDirectionPairs.rAngle ) < EPS_ANGLE
                   && fabs( vSecondPassCodeBook[snd_pass]->rRelativeDist_2 - trdTestDirectionPairs.rDist/sndTestDirectionPairs.rDist ) < EPS_DIST)
                  vThirdPassCodeBook.push_back( vSecondPassCodeBook[snd_pass] );
                  // here we can identify which distances and relative orientations are we looking at?
                  //std::cerr << "Should be looking for distance: " << vSecondPassCodeBook[snd_pass]->rRelativeDist_2 << " " << sndTestDirectionPairs.rDist << " " << std::endl;
              }
              if( vThirdPassCodeBook.size() == 0 ) continue;

              //numberOf4LongChains++;
              int beginFrtIndex = mvirDirectionIndexes[trdTestDirectionPairs.iIdx_0].x;
              if( beginFrtIndex == -1 ) continue;
              int endFrtIndex = mvirDirectionIndexes[trdTestDirectionPairs.iIdx_0].y;

              for( ; beginFrtIndex <= endFrtIndex; ++beginFrtIndex ) // search through fourth pair
              {
                delta_time = CVD::timer.get_time() - start_time;
               /* counterSteps[7]++;
                timestampsteps[7] += CVD::timer.get_time() - lastRecordedTime; 
                lastRecordedTime = CVD::timer.get_time();*/
 
                if( delta_time > MAX_TIME ) goto END;
#ifdef TEST_RESTART                
                if( delta_time > step_time + RESTART_TIME )
                {
                  //std::cout << delta_time <<" " << step_time << std::endl;
                  step_time += RESTART_TIME;
                  goto RESTART;
                }
#endif                
                if( vBinTestDirections[beginFrtIndex] )
                  continue;
                TestDirection & fourtTestDirectionPairs = mvTestDirections[beginFrtIndex];
                if( fabs( calculateAngle(fourtTestDirectionPairs.rDir_2,trdTestDirectionPairs.rDir_2) - mv4ChainAngles[3] ) < EPS_DIR_TEST && fourtTestDirectionPairs.rDist < MAX_DIST_TEST[multiScale_level] &&  fourtTestDirectionPairs.rDist > MIN_DIST_TEST[multiScale_level]) // fourth edge link
                {

                  //numberOfCompletedChains++;
                  // Check with codebook insdie vThirdPassCodeBook;

                  /*mvCandidateObjects.clear();
                  // this was added to check how many chains actually match BEFORE the distance transform error is calculated
                  int counterMatches = 0;
                  int counterHomMatches = 0;
                  for( unsigned int trd_pass = 0; trd_pass < vThirdPassCodeBook.size(); ++trd_pass )
                  {
                    if( fabs( vThirdPassCodeBook[trd_pass]->rAngle_3 - fourtTestDirectionPairs.rAngle ) < EPS_ANGLE 
                       && fabs( vThirdPassCodeBook[trd_pass]->rRelativeDist_3 - fourtTestDirectionPairs.rDist/trdTestDirectionPairs.rDist ) < EPS_DIST )
                    {
                       counterMatches++; 
                       CodeBookElement & cd = *(vThirdPassCodeBook[trd_pass]);
                       View & view = mvObjDescriptorClasses[no_class].mObjectsTemplate.vObjectElements[ cd.iObject_ID ].vViews[ cd.iView_ID ];
                       std::vector<Edgelet> & vEdgelets = view.vEdgelets;
                       Edges & allEdges = view.allEdges;
                       v5MatchedPairs[0].v2FstView = vEdgelets[cd.iIdx_0].v2Pose; v5MatchedPairs[0].v2SndView = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose;
                       v5MatchedPairs[1].v2FstView = vEdgelets[cd.iIdx_1].v2Pose; v5MatchedPairs[1].v2SndView = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose;
                       v5MatchedPairs[2].v2FstView = vEdgelets[cd.iIdx_2].v2Pose; v5MatchedPairs[2].v2SndView = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose;
                       v5MatchedPairs[3].v2FstView = vEdgelets[cd.iIdx_3].v2Pose; v5MatchedPairs[3].v2SndView = mvEdgelets[trdTestDirectionPairs.iIdx_1].v2Pose;
                       v5MatchedPairs[4].v2FstView = vEdgelets[cd.iIdx_4].v2Pose; v5MatchedPairs[4].v2SndView = mvEdgelets[fourtTestDirectionPairs.iIdx_1].v2Pose;
                       TooN::Matrix<3> h;
                       if ( mHomography.estimateHomography_group_of_5(h, v5MatchedPairs) )
                       {
                         counterHomMatches++;
                         std::vector<Edgelet> vProjectedEdgelets;
                         project_edgelets( vProjectedEdgelets, vEdgelets, h );
                         std::vector<TooN::Vector<2> > vv2Pnts;
                         Edges vProjectedEdges;
                         project_edges (vProjectedEdges, vv2Pnts, allEdges, h, false);
                         mvCandidateObjects.resize( mvCandidateObjects.size() + 1);
                         int obj_correct_ID = static_cast<int>(cd.obj_correct_ID);
                         mvCandidateObjects.back().iObject_ID = obj_correct_ID;

                         for (int ipt = 0; ipt < vv2Pnts.size(); ipt++){      
                            mvCandidateObjects.back().vv2Pnts.push_back (vv2Pnts[ipt]);
                         } 
                       }

                    }
                  } 
                  if (counterMatches > 0) 
                     std::cout << "number of matches with codebook " << counterMatches << " with homographies " << counterHomMatches << std::endl;                
                  */
                //  std::cerr << "Potential codebook matches " << vThirdPassCodeBook.size() << std::endl;
                  for( unsigned int trd_pass = 0; trd_pass < vThirdPassCodeBook.size(); ++trd_pass )
                  {
                    delta_time = CVD::timer.get_time() - start_time; 
                    if( delta_time > MAX_TIME ) goto END;
#ifdef TEST_RESTART                    
                    if( delta_time > step_time + RESTART_TIME )
                    {
                      //std::cout << delta_time <<" " << step_time << std::endl;
                      step_time += RESTART_TIME;
                      goto RESTART;
                    }
#endif                 
                    if( fabs( vThirdPassCodeBook[trd_pass]->rAngle_3 - fourtTestDirectionPairs.rAngle ) < EPS_ANGLE 
                       && fabs( vThirdPassCodeBook[trd_pass]->rRelativeDist_3 - fourtTestDirectionPairs.rDist/trdTestDirectionPairs.rDist ) < EPS_DIST )
                    { 
                     // std::cerr << "passed check in codebook for vThirdPassCodeBook" << std::endl;
                      CodeBookElement & cd = *(vThirdPassCodeBook[trd_pass]);
                      View & view = mvObjDescriptorClasses[no_class].mObjectsTemplate.vObjectElements[ cd.iObject_ID ].vViews[ cd.iView_ID ];
                      std::vector<Edgelet> & vEdgelets = view.vEdgelets;
                      Edges & allEdges = view.allEdges;
                      if (cd.iIdx_0 >= vEdgelets.size() || cd.iIdx_1 >= vEdgelets.size() || cd.iIdx_2 >= vEdgelets.size() || cd.iIdx_3 >= vEdgelets.size() || cd.iIdx_4 >= vEdgelets.size())
                      {
                          continue;
                      }
                      v5MatchedPairs[0].v2FstView = vEdgelets.at(cd.iIdx_0).v2Pose; v5MatchedPairs[0].v2SndView = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose;
                      v5MatchedPairs[1].v2FstView = vEdgelets.at(cd.iIdx_1).v2Pose; v5MatchedPairs[1].v2SndView = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose;
                      v5MatchedPairs[2].v2FstView = vEdgelets.at(cd.iIdx_2).v2Pose; v5MatchedPairs[2].v2SndView = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose;
                      v5MatchedPairs[3].v2FstView = vEdgelets.at(cd.iIdx_3).v2Pose; v5MatchedPairs[3].v2SndView = mvEdgelets[trdTestDirectionPairs.iIdx_1].v2Pose;
                      v5MatchedPairs[4].v2FstView = vEdgelets.at(cd.iIdx_4).v2Pose; v5MatchedPairs[4].v2SndView = mvEdgelets[fourtTestDirectionPairs.iIdx_1].v2Pose;
                      // add the chain
                      TooN::Vector<2> v1, v2, v3, v4, v5;
                      v1[0] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[0]*imgSizeRatio;
                      v1[1] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[1]*imgSizeRatio;
                      v2[0] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
                      v2[1] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
                      v3[0] = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
                      v3[1] = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
                      v4[0] = mvEdgelets[trdTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
                      v4[1] = mvEdgelets[trdTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
                      v5[0] = mvEdgelets[fourtTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
                      v5[1] = mvEdgelets[fourtTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
                      CheckedChains ch;
                      ch.chains.reserve (5);
                      ch.chains.push_back (v1);
                      ch.chains.push_back (v2);
                      ch.chains.push_back (v3);
                      ch.chains.push_back (v4);
                      ch.chains.push_back (v5);
                      mvCheckedChains.push_back (ch);
                      TooN::Matrix<3> h;
                      if ( mHomography.estimateHomography_group_of_5(h, v5MatchedPairs) )
                      {
                        std::vector<Edgelet> vProjectedEdgelets;
                        project_edgelets( vProjectedEdgelets, vEdgelets, h );
                        std::vector<REAL_TYPE> vrDistanceSum( vProjectedEdgelets.size() );

                        if (vProjectedEdgelets.size() == 0)
                           continue;
                        for( unsigned int no_edgelet = 0; no_edgelet < vProjectedEdgelets.size(); ++ no_edgelet )
                        {
                            if (isnan(vProjectedEdgelets[no_edgelet].rSlope))
                              continue;
                            int binNo = static_cast<int>( (atan( vProjectedEdgelets[no_edgelet].rSlope )/M_PI + 0.5f ) * 10.f );
                            vrDistanceSum[no_edgelet] = mvEdgeOrientationDistTransforms[binNo][vProjectedEdgelets[no_edgelet].v2Pose];
                        }
                        std::sort( vrDistanceSum.begin(), vrDistanceSum.end() );
                        REAL_TYPE rMedian = vrDistanceSum[vrDistanceSum.size()/2];
                       // std::cerr << "Before checking MAX_DT " << rMedian << " " << MAX_DT << std::endl;
                        if( rMedian < MAX_DT )
                        {
                          project_edgelets( vProjectedEdgelets, vEdgelets, h, false );

                          TooN::Matrix<3> newH;
                          std::vector<int> viFoundEdgeIndexes;
                          std::vector<int> viFullIndexes;
                          std::vector<MatchedPair> vMatchedPairs;
                          double angle;
                          REAL_TYPE rErr1 = iterativeClosestEdges2_NP( newH, viFullIndexes, viFoundEdgeIndexes, vProjectedEdgelets, mvEdgelets, vMatchedPairs, 25.f, 30.f, 0.2f);
                          std::vector<MatchedPair> matchedCorrespondences;
                          for (int i = 0; i < viFullIndexes.size(); i++)
                          {
                             if (viFullIndexes[i] == -1)
                                continue;
                             MatchedPair mp;
                             mp.v2FstView[0] = vEdgelets[i].v2Pose[0];
                             mp.v2FstView[1] = vEdgelets[i].v2Pose[1];
                             mp.v2SndView[0] = mvEdgelets[viFullIndexes[i]].v2Pose[0];
                             mp.v2SndView[1] = mvEdgelets[viFullIndexes[i]].v2Pose[1];
                             matchedCorrespondences.push_back (mp);
                          }
                          double stdDev;
                          angle = mHomography.estimateRotationAngle (stdDev, matchedCorrespondences);
                          //std::cerr << "Standard Deviation of Rotation Angle " << stdDev << " " << viFoundEdgeIndexes.size() << std::endl;
                          if (stdDev > 0.2)
                              continue;
                          std::vector<bool> vbBins( mvEdgelets.size(), true );
                          if( viFoundEdgeIndexes.size() > 5 )
                          {
                            int no_unique = 0;
                            std::vector<int>::iterator iIter = viFoundEdgeIndexes.begin();
                            std::vector<int>::iterator iEndIter = viFoundEdgeIndexes.end();
                            while( iIter != iEndIter )
                            {
                              if( vbBins[*iIter] )
                              {
                                vbBins[*iIter] = false;
                                ++no_unique;
                              }
                              ++iIter;
                            }
                              
                            REAL_TYPE rErr2 = rErr1*static_cast<REAL_TYPE>(viFoundEdgeIndexes.size())/(SAMPLING_FACTOR*static_cast<REAL_TYPE>(vMatchedPairs.size()));
                  //          std::cerr << "reached here with " << cd.obj_correct_ID << " " << rErr1 << " " << rErr2 << " " << (static_cast<REAL_TYPE>(vMatchedPairs.size())/static_cast<REAL_TYPE>(vEdgelets.size())) << std::endl;
                            
                            if( rErr1 < 0.12 && rErr2 < 0.05 &&  static_cast<REAL_TYPE>(vMatchedPairs.size())/static_cast<REAL_TYPE>(vEdgelets.size()) > 0.55)
                            {//FOUND
                            //  
                              double lastRecordedTime = CVD::timer.get_time(); 
                              current_time = CVD::timer.get_time();
                              std::string objName;
                              obj_correct_ID = static_cast<int>(cd.obj_correct_ID);
                              switch( obj_correct_ID )
                              {
                                case 3:
                                  objName = "Plier";
                                  break;
                                case 5:
                                  objName = "Box";
                                  break;
                                case 10:
                                  objName = "Charger";
                                  break;
                                case 4:
                                  objName = "Driver";
                                  break;
                                case 2:
                                  objName = "Hammer";
                                  break;
                                case 1:
                                  objName = "Screw Driver";
                                  break;
                                case 7:
                                  objName = "Stapler";
                                  break;
                                case 6:
                                  objName = "Wood";
                                  break;
                                case 8:
                                  objName = "Wrench";
                                  break;
                                case 9:
                                  objName = "Yellow";
                                  break;
                                default:
                                  objName = "Other";
                              }
                              h = newH*h;
                              project_edgelets( vProjectedEdgelets, vEdgelets, h, false );
                              mvDetectedObjects.resize( mvDetectedObjects.size() + 1);
                              mvDetectedObjects.back().iObject_ID = obj_correct_ID;
                              mvDetectedObjects.back().iView_ID = cd.iView_ID + 1;
                              // viFoundEdgeIndexes
                              std::vector<int> posRegionNo;
                              std::vector<int> posRegionNoCount;
                              bool found = false;
                              for (int vi = 0; vi < viFoundEdgeIndexes.size(); vi++)
                              {
                                 found = false;
                                 for (int ri = 0; ri < posRegionNo.size(); ri++)
                                    if (posRegionNo[ri] == mvEdgelets[viFoundEdgeIndexes[vi]].regionNo)
                                    {
                                        posRegionNoCount[ri]++;
                                        found = true;
                                        break;
                                    }
                                 if (!found)
                                 {
                                    posRegionNo.push_back(mvEdgelets[viFoundEdgeIndexes[vi]].regionNo);
                                    posRegionNoCount.push_back(1);
                                 } 
                              }
                              int max = 0;
                              int chosenRegionNo = -1;
                              for (int ri = 0; ri < posRegionNo.size(); ri++)
                                if (posRegionNoCount[ri] > max)
                                {
                                    max = posRegionNoCount[ri];
                                    chosenRegionNo = ri;
                                }
                              mvDetectedObjects.back().iRegion_No = chosenRegionNo;
                              Edges vProjectedEdges;
                              std::vector<TooN::Vector<2> > vv2Pnts;
                              project_edges (vProjectedEdges, vv2Pnts, allEdges, h, false);
                              for (int ipt = 0; ipt < vv2Pnts.size(); ipt++){      
                                  mvDetectedObjects.back().vv2Pnts.push_back (vv2Pnts[ipt]);
                              }
                              // change to the bigger size

                              TooN::Vector<2> y1, y2, y3, y4, y5;
           	              y1[0] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[0]*imgSizeRatio;
                    	      y1[1] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[1]*imgSizeRatio;
	                      y2[0] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
        	              y2[1] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
                	      y3[0] = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
	                      y3[1] = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
        	              y4[0] = mvEdgelets[trdTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
                	      y4[1] = mvEdgelets[trdTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
           	              y5[0] = mvEdgelets[fourtTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
                	      y5[1] = mvEdgelets[fourtTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
                              mvDetectedObjects.back().foundChainPnts.push_back(y1);
                              mvDetectedObjects.back().foundChainPnts.push_back(y2);
                              mvDetectedObjects.back().foundChainPnts.push_back(y3);
                              mvDetectedObjects.back().foundChainPnts.push_back(y4);
                              mvDetectedObjects.back().foundChainPnts.push_back(y5);
                              mvDetectedObjects.back().detection_time = lastRecordedTime - start_time;
                              std::cout << "Detection time : " << mvDetectedObjects.back().detection_time << " " << no_class << std::endl;
                              mvDetectedObjects.back().error = rErr2;

                              // homography
                              /*double trace = 1 + h(0,0) + h(1,1) + h(2,2);
                              if (trace > 1.000001)
                              {
                                 double s = 1/(sqrt (trace) * 2);
                                 mvDetectedObjects.back().pose_x = (h(1,2) - h(2,1)) / s;
                                 mvDetectedObjects.back().pose_y = (h(2,0) - h(0,2)) / s;
                                 mvDetectedObjects.back().pose_z = (h(1,0) - h(1,0)) / s;
                                 mvDetectedObjects.back().pose_w = trace * s;
                              }
                              else
                              {
                                 std::cout << "below trace = " << trace << std::endl;
                              }*/
                              mvDetectedObjects.back().pose_x = 0;
                              mvDetectedObjects.back().pose_y = 0;
                              mvDetectedObjects.back().pose_z = sin(angle/2);
                              mvDetectedObjects.back().pose_w = cos(angle/2);

                              // record the object ID in case a later view is added
                              latestObjID = cd.iObject_ID;
                              latest_correctObjID = obj_correct_ID;
                              
                              std::vector<Edgelet>::iterator tedIter = vProjectedEdgelets.begin();
                              CVD::ImageRef irTopLeft(img_h,img_w), irButtomRight(-1,-1);
                              while( tedIter != vProjectedEdgelets.end() )
                              {
                                if( (*tedIter).v2Pose[0] < irTopLeft.x )
                                  irTopLeft.x = (*tedIter).v2Pose[0];
                                else if( (*tedIter).v2Pose[0] > irButtomRight.x )
                                  irButtomRight.x = (*tedIter).v2Pose[0];
                                if( (*tedIter).v2Pose[1] < irTopLeft.y )
                                  irTopLeft.y = (*tedIter).v2Pose[1];
                                else if( (*tedIter).v2Pose[1] > irButtomRight.y )
                                  irButtomRight.y = (*tedIter).v2Pose[1];
                                ++tedIter;
                              }
                              mvDetectedObjects.back().irTopLeft = irTopLeft;
                              mvDetectedObjects.back().irBottomRight = irButtomRight;
                              count++;
                              if( current_time - subtract_time - start_time > MAX_TIME )
                                goto END;
                              // mark detected edges
                              #if 1
                              std::vector<TestDirection>::iterator bTestIter = mvTestDirections.begin();
                              std::vector<TestDirection>::iterator eTestIter = mvTestDirections.end();
                              std::vector<bool>::iterator bBinIter = vBinTestDirections.begin();
                              while( bTestIter != eTestIter )
                              {
                                if( !vbBins[bTestIter->iIdx_0] )
                                {
                                  *bBinIter = true;
                                  mvEdgelets[bTestIter->iIdx_0].v2Pose = TooN::Zeros;
                                }
                                else if ( !vbBins[bTestIter->iIdx_1] )
                                {
                                  *bBinIter = true;
                                  mvEdgelets[bTestIter->iIdx_1].v2Pose = TooN::Zeros;
                                }
                                ++bTestIter;
                                ++bBinIter;
                              }
                              #endif

                           //   goto NEXT_ROUND;
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      NEXT_ROUND:;
    }
    }
  }
  END:

  // time analysis stages

/*  for (int i = 1; i < NumberOfSteps; i++)
  {
     std::cout << "Counter of step " << i << " is " << counterSteps[i] << ", average = " << timestampsteps[i]/counterSteps[i] << ", total = " << timestampsteps[i] << std::endl;
  }*/

//  std::cout << "number of pairs tested " << numberOfStartingEdges << " out of " << numberOfPairs << std::endl;
//  std::cout << "number of chains of length 2 " << numberOf2LongChains << std::endl;
//  std::cout << "number of chains of length 3 " << numberOf3LongChains << std::endl;
//  std::cout << "number of chains of length 4 " << numberOf4LongChains << std::endl;
//  std::cout << "number of full length chains " << numberOfCompletedChains << std::endl;
//  std::cout << "average time for chains " << chainDuration/numberOfPairs << std::endl;
  
  logFile.close();
 // std::cout << "exiting with ********** " << mvDetectedObjects.size() << " and checked chains " << mvCheckedChains.size () << std::endl;
 // for (int i = 0; i < mvDetectedObjects.size(); i++)
 //    std::cout << "Object " << i << " took " << mvDetectedObjects[i].detection_time << " seconds" << std::endl;
  return ( mvDetectedObjects.size() != 0) ? true : false;
}

int EdgeObjectDetection::detect_constellation(int subtract_time)
{
  std::ofstream logFile;
  logFile.open( outputFileName, std::ios::app );
  double current_time;
  int count=0;
  int count_2 = 0;
  int beginIdx = 0;
  int endIdx;
  int numberOfStartingEdges = 0;
  int numberOf2LongChains = 0;
  int numberOf3LongChains = 0;
  int numberOf4LongChains = 0;
  int numberOfCompletedChains = 0;
  double chainDuration = 0;
  int numberOfPairs;
  for (int i = 0; i < NumberOfSteps; i++)
  {
     counterSteps[i] = 0;
     timestampsteps[i] = 0;
  }

  std::vector<MatchedPair> v5MatchedPairs(5);
  for( unsigned int no_class = 0; no_class < mvObjDescriptorClasses.size(); ++no_class )
  {
    TooN::Vector<4> & mv4ChainAngles = mvObjDescriptorClasses[no_class].mv4ChainAngles;
    std::vector<CodeBookElement> & mvCodeBook = mvObjDescriptorClasses[no_class].mvCodeBook;
    std::vector<CVD::ImageRef> & mvirFirstAngleIndexes = mvObjDescriptorClasses[no_class].mvirFirstAngleIndexes;
    std::vector< std::vector<int> > & mvirSecondIndexes = mvObjDescriptorClasses[no_class].mvirSecondIndexes;

#if 0
// This block may help speed up the program if there are more than one class ( fixed angle chains ).
    std::vector<std::pair<REAL_TYPE, int> > mvpairDir_1_Indexes( mvTestDirections.size() ); // having mvpairDir_1_Indexes will help if there are more classes ( more v4ChainAngles ) to be used.
    for( unsigned int no_testdirection = 0; no_testdirection < mvTestDirections.size(); ++ no_testdirection )
      mvpairDir_1_Indexes[no_testdirection] = std::make_pair( mvTestDirections[no_testdirection].rDir_1, no_testdirection );
    std::sort( mvpairDir_1_Indexes.begin(), mvpairDir_1_Indexes.end(), comp_pair_r_i );
    // | dir - mv4ChainAngles[0] | < EPS_DIR_TEST;
    // first part
    REAL_TYPE rUpper_bound = mv4ChainAngles[0] + EPS_DIR_TEST; // first edge link
    REAL_TYPE rLower_bound = mv4ChainAngles[0] - EPS_DIR_TEST;

    std::vector<std::pair<REAL_TYPE, int> >::iterator fstBeginPairIter = std::lower_bound( mvpairDir_1_Indexes.begin(), mvpairDir_1_Indexes.end(), rLower_bound,   comp_lower_pair );
    std::vector<std::pair<REAL_TYPE, int> >::iterator fstEndPairIter = std::upper_bound( fstBeginPairIter, mvpairDir_1_Indexes.end(), rUpper_bound, comp_upper_pair );
    for( ; fstBeginPairIter != fstEndPairIter; ++fstBeginPairIter ) // after passing the first check from mv4ChainAngles[0], 
    {
//      std::cout << fstBeginPairIter->first << " " << fstBeginPairIter->second << std::endl;
    // then check for angle between the pair and elements inside codebooks.
      TestDirection & firstTestDirectionPairs = mvTestDirections[(*fstBeginPairIter).second];
      if( firstTestDirectionPairs.rDist > FIRST_DIST_TRAIN )
        continue;
      count_2 ++;
#else
  #ifdef ECM
  int sum2 = 0, posTD, edgeNo2;

 // double howLong1 = CVD::timer.get_time();
  // here I need to search through the edges only not the pairs
  for (int i = 0; i < mvCellDirections.size(); i++)
  {     
     CellDirection thisCD = mvCellDirections[i];
     if (fabs(thisCD.rDir_1 - mv4ChainAngles[0]) < (EPS_DIR_TEST) || fabs(thisCD.rDir_1 + mv4ChainAngles[0]) < (EPS_DIR_TEST))
     {
        for (int thisCellEdgelets = 0; thisCellEdgelets < ecm[thisCD.cellNo].cellEdgelets.size(); thisCellEdgelets++)
        {
           int edgeNo2 = ecm[thisCD.cellNo].cellEdgelets[thisCellEdgelets];
           if (thisCD.iIdx == edgeNo2)
              continue;
           if (edgeNo2 < thisCD.iIdx)
              posTD = thisCD.iIdx*(mvEdgelets.size()-1) + edgeNo2; 
           else
              posTD = thisCD.iIdx*(mvEdgelets.size()-1) + edgeNo2 - 1;
           vCellTestDirections.push_back (mvTestDirections[posTD]);
           sum2++;
        }
     }
  }
  std::vector<bool> vBinTestDirections(vCellTestDirections.size(),false);
  std::vector<int> vRandomSamples(vCellTestDirections.size());
  for( unsigned int i = 0; i < vCellTestDirections.size(); ++ i ) vRandomSamples[i]=i;  
  endIdx = sum2 - 1;
  numberOfPairs = sum2;
  #else
  std::vector<bool> vBinTestDirections(mvTestDirections.size(),false);
  std::vector<int> vRandomSamples(mvTestDirections.size());
  for( unsigned int i = 0; i < mvTestDirections.size(); ++ i ) vRandomSamples[i]=i;
  endIdx = mvTestDirections.size() - 1;
  numberOfPairs = mvTestDirections.size();
  #endif

  mvCheckedChains.clear();
  double delta_time;
  double step_time = 0;
 // double lastRecordedTime = start_time;

  RESTART:
 // std::cout <<"Shuffle ..." << std::endl;
  std::random_shuffle( vRandomSamples.begin(), vRandomSamples.end() );
  for (int multiScale_level = 0; multiScale_level <= MAX_MULTISCALE_LEVEL; multiScale_level++)
  {
  //std::cerr << "Entering Multi-Scale level = " << multiScale_level << std::endl;
  beginIdx = 0;
  for( ; beginIdx <= endIdx; ++ beginIdx )
  {
    chainDuration += CVD::timer.get_time() - delta_time;
    numberOfStartingEdges++;
    delta_time = CVD::timer.get_time() - start_time; 
  /*  counterSteps[1]++;
    timestampsteps[1] += CVD::timer.get_time() - lastRecordedTime; 
    lastRecordedTime = CVD::timer.get_time();*/
    if( delta_time > MAX_TIME ) goto END;
#ifdef TEST_RESTART    
    if( delta_time > step_time + RESTART_TIME )
    {
      //std::cout << delta_time <<" " << step_time << std::endl;
      step_time += RESTART_TIME;
      goto RESTART;
    }
#endif    
    if( vBinTestDirections[vRandomSamples[beginIdx]] )
      continue;
    #ifdef ECM
    TestDirection & firstTestDirectionPairs = vCellTestDirections[vRandomSamples[beginIdx]];
    #else
    TestDirection & firstTestDirectionPairs = mvTestDirections[vRandomSamples[beginIdx]];
    #endif
    count_2 ++;


    if(!( (fabs( firstTestDirectionPairs.rDir_1-mv4ChainAngles[0] ) < EPS_DIR_TEST || fabs( firstTestDirectionPairs.rDir_1+mv4ChainAngles[0] ) < EPS_DIR_TEST ) && firstTestDirectionPairs.rDist < MAX_DIST_TEST[multiScale_level] &&  firstTestDirectionPairs.rDist > MIN_DIST_TEST[multiScale_level]) )
//    if( !((fabs( firstTestDirectionPairs.rDir_1 - mv4ChainAngles[0] ) < EPS_DIR_TEST) && firstTestDirectionPairs.rDist < FIRST_DIST_TRAIN ))
      continue;
   /* CheckedChains ch;
    ch.chains.reserve (2);
    TooN::Vector<2> v1, v2;
    v1[0] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[0]*imgSizeRatio;
    v1[1] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[1]*imgSizeRatio;
    v2[0] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
    v2[1] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;

    ch.chains.push_back (v1);
    ch.chains.push_back (v2);
    mvCheckedChains.push_back (ch);*/

     // try to get rid of parallel edges?
     if ((firstTestDirectionPairs.rAngle < 0.05 && firstTestDirectionPairs.rAngle > -0.05) ||  (fabs(firstTestDirectionPairs.rAngle - M_PI) < 0.05))
     {
       // std::cerr << "dropping pair with angle " << firstTestDirectionPairs.rAngle << std::endl;
        continue;
     }
#endif
      // looking for a second pair that fst.e2 == snd.e1
      
      int beginSndIndex = mvirDirectionIndexes[firstTestDirectionPairs.iIdx_1].x;
      if( beginSndIndex == -1 ) continue;
      int endSndIndex = mvirDirectionIndexes[firstTestDirectionPairs.iIdx_1].y;
      // searching for a first pair inside codebook.
    /*  counterSteps[2]++;
      timestampsteps[2] += CVD::timer.get_time() - lastRecordedTime; 
      lastRecordedTime = CVD::timer.get_time();*/
      REAL_TYPE rLower_Angle = firstTestDirectionPairs.rAngle - EPS_ANGLE;
      REAL_TYPE rUpper_Angle = firstTestDirectionPairs.rAngle + EPS_ANGLE;
      returnedEdgelets.clear();
      returnedCB.clear();

      for( ; beginSndIndex <= endSndIndex; ++beginSndIndex ) // search through second pair
      {
        delta_time = CVD::timer.get_time() - start_time; 
      /*  counterSteps[3]++;
        timestampsteps[3] += CVD::timer.get_time() - lastRecordedTime; 
        lastRecordedTime = CVD::timer.get_time();*/

        if( delta_time > MAX_TIME ) goto END;
#ifdef TEST_RESTART        
        if( delta_time > step_time + RESTART_TIME )
        {
          //std::cout << delta_time <<" " << step_time << std::endl;
          step_time += RESTART_TIME;
          goto RESTART;
        }
#endif        
        if( vBinTestDirections[beginSndIndex] )
          continue;
        const TestDirection & sndTestDirectionPairs = mvTestDirections[beginSndIndex];
        if( fabs( calculateAngle(sndTestDirectionPairs.rDir_2,firstTestDirectionPairs.rDir_2) - mv4ChainAngles[1] ) < EPS_DIR_TEST && sndTestDirectionPairs.rDist < MAX_DIST_TEST[multiScale_level] &&  sndTestDirectionPairs.rDist > MIN_DIST_TEST[multiScale_level]) // second edge lin k
        {
/*CheckedChains ch;
		ch.chains.reserve (3);
		 TooN::Vector<2> v1, v2, v3;
		 v1[0] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[0]*imgSizeRatio;
		 v1[1] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[1]*imgSizeRatio;
		 v2[0] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
		 v2[1] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
		 v3[0] = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
		 v3[1] = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;

		ch.chains.push_back (v1);
		ch.chains.push_back (v2);
		ch.chains.push_back (v3);
		mvCheckedChains.push_back (ch);*/
        // Check with codebook betweeen lowerCodeBook and upperCodeBook 
          std::vector< std::vector<CodeBookElement>::iterator > vSecondPassCodeBook;
      /*    vSecondPassCodeBook.reserve( upperCodeBook - lowerCodeBook );
          lowerCodeBook = minCodeBook;
          for( ; lowerCodeBook != upperCodeBook; ++lowerCodeBook )
          {
            thisCB.push_back(lowerCodeBook);
            if( fabs( (*lowerCodeBook).rAngle_1 - sndTestDirectionPairs.rAngle ) < EPS_ANGLE 
               && fabs( (*lowerCodeBook).rRelativeDist_1 - sndTestDirectionPairs.rDist/firstTestDirectionPairs.rDist ) < EPS_DIST )
              vSecondPassCodeBook.push_back( lowerCodeBook );
          }*/
          // trying to read from the codebook

          int indexNo = 64.f*(firstTestDirectionPairs.rAngle/M_PI + 0.5f); if( indexNo < 0 ) indexNo = 0;
          int indexNo2 = 64.f*(sndTestDirectionPairs.rAngle/M_PI + 0.5f); if( indexNo2 < 0 ) indexNo2 = 0;
          int indexNo3;
          REAL_TYPE ratio = sndTestDirectionPairs.rDist/firstTestDirectionPairs.rDist;
          if (ratio < 1)
              indexNo3 = int(ratio/0.1);
          else if (ratio > 1 && ratio < 29.5)
              indexNo3 = int (ratio)+9;
          else
              indexNo3 = 38;
          std::vector<int> & thisList = mvirSecondIndexes[(indexNo)*2496+(indexNo2)*39+indexNo3];
          vSecondPassCodeBook.reserve( thisList.size() );
         // int distBin [25] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
          for(int k = 0; k < thisList.size(); k++)
          {
            std::vector<CodeBookElement>::iterator lx = mvCodeBook.begin()+thisList[k];
            //thisCB.push_back(lx);
            if( fabs( (*lx).rAngle_1 - sndTestDirectionPairs.rAngle ) < EPS_ANGLE 
               && fabs( (*lx).rRelativeDist_1 - sndTestDirectionPairs.rDist/firstTestDirectionPairs.rDist ) < EPS_DIST )
             {
                 vSecondPassCodeBook.push_back( lx );
                 REAL_TYPE lookingforDist = (*lx).rRelativeDist_2 * sndTestDirectionPairs.rDist;
              //   int binNo = (int)(lookingforDist/10);
              //   distBin[binNo]++;
             }
          }   

          if( vSecondPassCodeBook.size() == 0 )
          {
             continue;
          }
          //numberOf3LongChains++;
          int beginTrdIndex = mvirDirectionIndexes[sndTestDirectionPairs.iIdx_1].x;
          if( beginTrdIndex == -1 ) continue;
          int endTrdIndex = mvirDirectionIndexes[sndTestDirectionPairs.iIdx_1].y;
          
          for( ; beginTrdIndex <= endTrdIndex; ++beginTrdIndex ) // search through third pair
          {
            delta_time = CVD::timer.get_time() - start_time;
          /*  counterSteps[5]++;
           timestampsteps[5] += CVD::timer.get_time() - lastRecordedTime; 
           lastRecordedTime = CVD::timer.get_time();*/

            if( delta_time > MAX_TIME ) goto END;
#ifdef TEST_RESTART            
            if( delta_time > step_time + RESTART_TIME )
            {
              step_time += RESTART_TIME;
              goto RESTART;
            }
#endif            
            if( vBinTestDirections[beginTrdIndex] )
              continue;
            TestDirection & trdTestDirectionPairs = mvTestDirections[beginTrdIndex];
            if( fabs( calculateAngle(trdTestDirectionPairs.rDir_2, sndTestDirectionPairs.rDir_2) - mv4ChainAngles[2] ) < EPS_DIR_TEST && trdTestDirectionPairs.rDist < MAX_DIST_TEST[multiScale_level] &&  trdTestDirectionPairs.rDist > MIN_DIST_TEST[multiScale_level]) // third edge link
            {
/*              int binNo = (int)(trdTestDirectionPairs.rDist/10);
              int countBins = distBin[binNo];
              if (binNo > 0)
                  binNo += distBin[binNo-1];
              if (binNo < 24)
                  binNo += distBin[binNo+1];*/
              
            //  std::cerr << ((binNo > 0) ? distBin[binNo-1] : "") << " " << distBin[binNo] << " " << ((binNo < 24) ? distBin[binNo+1] : "") << std::endl;
              std::vector< std::vector<CodeBookElement>::iterator > vThirdPassCodeBook;
              vThirdPassCodeBook.reserve( vSecondPassCodeBook.size() );
              for( unsigned int snd_pass = 0; snd_pass < vSecondPassCodeBook.size(); ++snd_pass )
              {
                if( fabs( vSecondPassCodeBook[snd_pass]->rAngle_2 - trdTestDirectionPairs.rAngle ) < EPS_ANGLE
                   && fabs( vSecondPassCodeBook[snd_pass]->rRelativeDist_2 - trdTestDirectionPairs.rDist/sndTestDirectionPairs.rDist ) < EPS_DIST)
                  vThirdPassCodeBook.push_back( vSecondPassCodeBook[snd_pass] );
                  // here we can identify which distances and relative orientations are we looking at?
                  //std::cerr << "Should be looking for distance: " << vSecondPassCodeBook[snd_pass]->rRelativeDist_2 << " " << sndTestDirectionPairs.rDist << " " << std::endl;
              }
              if( vThirdPassCodeBook.size() == 0 ) continue;

              //numberOf4LongChains++;
              int beginFrtIndex = mvirDirectionIndexes[trdTestDirectionPairs.iIdx_1].x;
              if( beginFrtIndex == -1 ) continue;
              int endFrtIndex = mvirDirectionIndexes[trdTestDirectionPairs.iIdx_1].y;

              for( ; beginFrtIndex <= endFrtIndex; ++beginFrtIndex ) // search through fourth pair
              {
                delta_time = CVD::timer.get_time() - start_time;
               /* counterSteps[7]++;
                timestampsteps[7] += CVD::timer.get_time() - lastRecordedTime; 
                lastRecordedTime = CVD::timer.get_time();*/
 
                if( delta_time > MAX_TIME ) goto END;
#ifdef TEST_RESTART                
                if( delta_time > step_time + RESTART_TIME )
                {
                  //std::cout << delta_time <<" " << step_time << std::endl;
                  step_time += RESTART_TIME;
                  goto RESTART;
                }
#endif                
                if( vBinTestDirections[beginFrtIndex] )
                  continue;
                TestDirection & fourtTestDirectionPairs = mvTestDirections[beginFrtIndex];
                if( fabs( calculateAngle(fourtTestDirectionPairs.rDir_2,trdTestDirectionPairs.rDir_2) - mv4ChainAngles[3] ) < EPS_DIR_TEST && fourtTestDirectionPairs.rDist < MAX_DIST_TEST[multiScale_level] &&  fourtTestDirectionPairs.rDist > MIN_DIST_TEST[multiScale_level]) // fourth edge link
                {

                  //numberOfCompletedChains++;
                  // Check with codebook insdie vThirdPassCodeBook;

                  /*mvCandidateObjects.clear();
                  // this was added to check how many chains actually match BEFORE the distance transform error is calculated
                  int counterMatches = 0;
                  int counterHomMatches = 0;
                  for( unsigned int trd_pass = 0; trd_pass < vThirdPassCodeBook.size(); ++trd_pass )
                  {
                    if( fabs( vThirdPassCodeBook[trd_pass]->rAngle_3 - fourtTestDirectionPairs.rAngle ) < EPS_ANGLE 
                       && fabs( vThirdPassCodeBook[trd_pass]->rRelativeDist_3 - fourtTestDirectionPairs.rDist/trdTestDirectionPairs.rDist ) < EPS_DIST )
                    {
                       counterMatches++; 
                       CodeBookElement & cd = *(vThirdPassCodeBook[trd_pass]);
                       View & view = mvObjDescriptorClasses[no_class].mObjectsTemplate.vObjectElements[ cd.iObject_ID ].vViews[ cd.iView_ID ];
                       std::vector<Edgelet> & vEdgelets = view.vEdgelets;
                       Edges & allEdges = view.allEdges;
                       v5MatchedPairs[0].v2FstView = vEdgelets[cd.iIdx_0].v2Pose; v5MatchedPairs[0].v2SndView = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose;
                       v5MatchedPairs[1].v2FstView = vEdgelets[cd.iIdx_1].v2Pose; v5MatchedPairs[1].v2SndView = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose;
                       v5MatchedPairs[2].v2FstView = vEdgelets[cd.iIdx_2].v2Pose; v5MatchedPairs[2].v2SndView = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose;
                       v5MatchedPairs[3].v2FstView = vEdgelets[cd.iIdx_3].v2Pose; v5MatchedPairs[3].v2SndView = mvEdgelets[trdTestDirectionPairs.iIdx_1].v2Pose;
                       v5MatchedPairs[4].v2FstView = vEdgelets[cd.iIdx_4].v2Pose; v5MatchedPairs[4].v2SndView = mvEdgelets[fourtTestDirectionPairs.iIdx_1].v2Pose;
                       TooN::Matrix<3> h;
                       if ( mHomography.estimateHomography_group_of_5(h, v5MatchedPairs) )
                       {
                         counterHomMatches++;
                         std::vector<Edgelet> vProjectedEdgelets;
                         project_edgelets( vProjectedEdgelets, vEdgelets, h );
                         std::vector<TooN::Vector<2> > vv2Pnts;
                         Edges vProjectedEdges;
                         project_edges (vProjectedEdges, vv2Pnts, allEdges, h, false);
                         mvCandidateObjects.resize( mvCandidateObjects.size() + 1);
                         int obj_correct_ID = static_cast<int>(cd.obj_correct_ID);
                         mvCandidateObjects.back().iObject_ID = obj_correct_ID;

                         for (int ipt = 0; ipt < vv2Pnts.size(); ipt++){      
                            mvCandidateObjects.back().vv2Pnts.push_back (vv2Pnts[ipt]);
                         } 
                       }

                    }
                  } 
                  if (counterMatches > 0) 
                     std::cout << "number of matches with codebook " << counterMatches << " with homographies " << counterHomMatches << std::endl;                
                  */
                //  std::cerr << "Potential codebook matches " << vThirdPassCodeBook.size() << std::endl;
                  for( unsigned int trd_pass = 0; trd_pass < vThirdPassCodeBook.size(); ++trd_pass )
                  {
                    delta_time = CVD::timer.get_time() - start_time; 
                    if( delta_time > MAX_TIME ) goto END;
#ifdef TEST_RESTART                    
                    if( delta_time > step_time + RESTART_TIME )
                    {
                      //std::cout << delta_time <<" " << step_time << std::endl;
                      step_time += RESTART_TIME;
                      goto RESTART;
                    }
#endif                 
                    if( fabs( vThirdPassCodeBook[trd_pass]->rAngle_3 - fourtTestDirectionPairs.rAngle ) < EPS_ANGLE 
                       && fabs( vThirdPassCodeBook[trd_pass]->rRelativeDist_3 - fourtTestDirectionPairs.rDist/trdTestDirectionPairs.rDist ) < EPS_DIST )
                    { 
                     // std::cerr << "passed check in codebook for vThirdPassCodeBook" << std::endl;
                      CodeBookElement & cd = *(vThirdPassCodeBook[trd_pass]);
                      View & view = mvObjDescriptorClasses[no_class].mObjectsTemplate.vObjectElements[ cd.iObject_ID ].vViews[ cd.iView_ID ];
                      std::vector<Edgelet> & vEdgelets = view.vEdgelets;
                      Edges & allEdges = view.allEdges;
                      if (cd.iIdx_0 >= vEdgelets.size() || cd.iIdx_1 >= vEdgelets.size() || cd.iIdx_2 >= vEdgelets.size() || cd.iIdx_3 >= vEdgelets.size() || cd.iIdx_4 >= vEdgelets.size())
                      {
                          continue;
                      }
                      v5MatchedPairs[0].v2FstView = vEdgelets.at(cd.iIdx_0).v2Pose; v5MatchedPairs[0].v2SndView = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose;
                      v5MatchedPairs[1].v2FstView = vEdgelets.at(cd.iIdx_1).v2Pose; v5MatchedPairs[1].v2SndView = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose;
                      v5MatchedPairs[2].v2FstView = vEdgelets.at(cd.iIdx_2).v2Pose; v5MatchedPairs[2].v2SndView = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose;
                      v5MatchedPairs[3].v2FstView = vEdgelets.at(cd.iIdx_3).v2Pose; v5MatchedPairs[3].v2SndView = mvEdgelets[trdTestDirectionPairs.iIdx_1].v2Pose;
                      v5MatchedPairs[4].v2FstView = vEdgelets.at(cd.iIdx_4).v2Pose; v5MatchedPairs[4].v2SndView = mvEdgelets[fourtTestDirectionPairs.iIdx_1].v2Pose;
                      // add the chain
                      TooN::Vector<2> v1, v2, v3, v4, v5;
                      v1[0] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[0]*imgSizeRatio;
                      v1[1] = mvEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose[1]*imgSizeRatio;
                      v2[0] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
                      v2[1] = mvEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
                      v3[0] = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
                      v3[1] = mvEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
                      v4[0] = mvEdgelets[trdTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
                      v4[1] = mvEdgelets[trdTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
                      v5[0] = mvEdgelets[fourtTestDirectionPairs.iIdx_1].v2Pose[0]*imgSizeRatio;
                      v5[1] = mvEdgelets[fourtTestDirectionPairs.iIdx_1].v2Pose[1]*imgSizeRatio;
                      CheckedChains ch;
                      ch.chains.reserve (5);
                      ch.chains.push_back (v1);
                      ch.chains.push_back (v2);
                      ch.chains.push_back (v3);
                      ch.chains.push_back (v4);
                      ch.chains.push_back (v5);
                      mvCheckedChains.push_back (ch);
                    }
                  }
                }
              }
            }
          }
        }
      }
      NEXT_ROUND:;
    }
    }
  }
  END:

  // time analysis stages

/*  for (int i = 1; i < NumberOfSteps; i++)
  {
     std::cout << "Counter of step " << i << " is " << counterSteps[i] << ", average = " << timestampsteps[i]/counterSteps[i] << ", total = " << timestampsteps[i] << std::endl;
  }*/

//  std::cout << "number of pairs tested " << numberOfStartingEdges << " out of " << numberOfPairs << std::endl;
//  std::cout << "number of chains of length 2 " << numberOf2LongChains << std::endl;
//  std::cout << "number of chains of length 3 " << numberOf3LongChains << std::endl;
//  std::cout << "number of chains of length 4 " << numberOf4LongChains << std::endl;
//  std::cout << "number of full length chains " << numberOfCompletedChains << std::endl;
//  std::cout << "average time for chains " << chainDuration/numberOfPairs << std::endl;
  
  logFile.close();
 // std::cout << "exiting with ********** " << mvDetectedObjects.size() << " and checked chains " << mvCheckedChains.size () << std::endl;
 // for (int i = 0; i < mvDetectedObjects.size(); i++)
 //    std::cout << "Object " << i << " took " << mvDetectedObjects[i].detection_time << " seconds" << std::endl;
  return ( mvDetectedObjects.size() != 0) ? true : false;
}

int EdgeObjectDetection::max (int a, int b)
{
   if (a > b) return a;
   return b;
}

int EdgeObjectDetection::min (int a, int b)
{
   if (a < b) return a;
   return b;
}

void EdgeObjectDetection::preCalculateEdgeLnksWithRegions( LowLevelImageData<CVD::byte> & llImage, std::vector<SegmentedRegion> imageRegions, std::vector<SegmentedMask> imageMasks, int mfi)
{
  std::srand( time(NULL) );
  start_time = CVD::timer.get_time();
  std::vector<Line> vLines;
  master_file_id = mfi;
  double line_time;
  thisLlImage = llImage;
  // copy detected objects onto previous frame's detected objects
  mvDetectedObjects_previousFrame.clear();
  for (int i = 0; i < mvDetectedObjects.size(); i++)
  {
     mvDetectedObjects[i].missingCount = 0;
     mvDetectedObjects[i].irTopLeft.x -= 5;
     mvDetectedObjects[i].irTopLeft.y -= 5;
     mvDetectedObjects[i].irBottomRight.x += 5;
     mvDetectedObjects[i].irBottomRight.y += 5;
     mvDetectedObjects_previousFrame.push_back(mvDetectedObjects[i]);
  }
  for (int i = 0; i < mvDetectedObjects_undetected.size(); i++)
  {
/*    bool overlap = false;
    for (int j = 0; j < mvDetectedObjects.size(); i++)
    {
       int a_x1 = min (mvDetectedObjects[j].irTopLeft.x, mvDetectedObjects_undetected[i].irTopLeft.x);
       int a_y1 = min (mvDetectedObjects[j].irTopLeft.y, mvDetectedObjects_undetected[i].irTopLeft.y);
       int b_x1 = max (mvDetectedObjects[j].irTopLeft.x, mvDetectedObjects_undetected[i].irTopLeft.x);
       int b_y1 = max (mvDetectedObjects[j].irTopLeft.y, mvDetectedObjects_undetected[i].irTopLeft.y);
       int a_x2 = max (mvDetectedObjects[j].irBottomRight.x, mvDetectedObjects_undetected[i].irBottomRight.x);
       int a_y2 = max (mvDetectedObjects[j].irBottomRight.y, mvDetectedObjects_undetected[i].irBottomRight.y);
       int b_x2 = min (mvDetectedObjects[j].irBottomRight.x, mvDetectedObjects_undetected[i].irBottomRight.x);
       int b_y2 = min (mvDetectedObjects[j].irBottomRight.y, mvDetectedObjects_undetected[i].irBottomRight.y);
       double a1 = (a_y2-a_y1)*(a_x2-a_x1);
       double a2 = (b_y2-b_y1)*(b_x2-b_x1);
       if (a1 < 0 || a2 < 0)
         continue;
       double overlap = a1/a2;
       if (overlap > 0.6)
       {
          overlap = true;
          break;
       }
    }
    if (!overlap)
    {*/
       mvDetectedObjects_undetected[i].missingCount++;
       mvDetectedObjects_previousFrame.push_back (mvDetectedObjects_undetected[i]);
    //}
  }
  mvDetectedObjects_undetected.clear();
  mvDetectedObjects.clear();
  if (!isLSD)
  {
    //Edges edges;
    mCannyEdgeDetector.compute( mEdgeImage, edges, llImage, 10, 0.5, 0.2 );
    /* DIMA - tests to adjust the edgelet size and line length automatically */
    int edgeletNumber = (int)edges.virEdges.size();
    std::cerr << "Edgelet Number: " << edgeletNumber << std::endl;
    
    if (edgeletNumber > 3000){
	    lineMaxLen = (float)LINE_MIN_LEN * (float)2;
	    EDGELET_LENGTH_TEST = (float)EDGELET_MIN_LEN * 8;
    }
    if( edgeletNumber > 1500){
	    lineMaxLen = (float)LINE_MIN_LEN * (float)2;
	    EDGELET_LENGTH_TEST = (float)EDGELET_MIN_LEN * 4;
    }
    else
    {
	    lineMaxLen = (float)LINE_MIN_LEN;
	    EDGELET_LENGTH_TEST = (float)EDGELET_MIN_LEN;
    }
    std::cerr << "New EDGELET_LENGTH_TEST: " << EDGELET_LENGTH_TEST << std::endl;
    get_lines( vLines, edges );
    line_time = CVD::timer.get_time();
  /*  char filename[100];
    sprintf(filename,"%s/master_%d.xvy", folderName, master_file_id);
    std::ofstream out_xvy(filename);
    for( unsigned int i = 0; i < edges.virEdges.size(); ++i )
    {
      mEdgeImage[edges.virEdges[i]] = 255;
      out_xvy << CVD::vec(edges.virEdges[i]) << std::endl;
    }
    out_xvy.close();
    sprintf(filename,"%s/master_%d.png",folderName, master_file_id);
    std::cout << "LSD " << filename << std::endl;
    //++master_file_id;
    CVD::img_save(mEdgeImage, filename);*/
  }
  else
  {
    int no_data = 0;
    double * img_ptr = lsd_image->data;
    CVD::byte * int_img_ptr = llImage.mImage.data();
    while( no_data != data_size )
    {
      *img_ptr = static_cast<double>(*int_img_ptr);
      ++ no_data;
      ++ img_ptr;
      ++ int_img_ptr;
    }
    ntuple_list out;
    out = lsd(lsd_image);
    lsd_to_lines( vLines, out );
    lines_to_edges (edges, vLines);
    free_ntuple_list(out);
    line_time = CVD::timer.get_time();
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
      std::cout << "CANNY " << filename << std::endl;
      CVD::img_save(mEdgeImage, filename);
    #endif
  }
  // calculate the overlap between image regions
  //std::cerr << "Reached here with " << EDGELET_LENGTH_TEST << std::endl;
  if (isMask)
     gen_Edgelets_with_Directional_Masks( mvEdgelets, mviLineIndexes, vLines, imageMasks, EDGELET_LENGTH_TEST);     
  else
     gen_Edgelets_with_Directional_Regions( mvEdgelets, mviLineIndexes, vLines, imageRegions, 0, EDGELET_LENGTH_TEST);
  double edgelet_time = CVD::timer.get_time();
  gen_TestDirections( mvEdgelets, mviLineIndexes );
  #ifdef ECM
  gen_EdgeletCells (mvEdgelets);
  double beforeECM = CVD::timer.get_time ();
  EdgeCellMatrix* ecm =  getEdgeletImageMatrix (mvEdgelets, img_h, img_w, 10, 10);
  double afterECM = CVD::timer.get_time ();
  #endif
  double td_time = CVD::timer.get_time();
}

void EdgeObjectDetection::getDistanceTransformImage (CVD::Image<CVD::Rgb<CVD::byte> >& transformImage)
{
    int width = mvEdgeDistanceTransformImage.size().x;
    int height = mvEdgeDistanceTransformImage.size().y;
    for( int row = 0; row < height; ++row )
    {
      for( int col = 0; col < width; ++col )
      {
         transformImage[row*(int)(SAMPLING_FACTOR)][col*(int)(SAMPLING_FACTOR)].red = (CVD::byte)(mvEdgeDistanceTransformImage[row][col]*255.f);
         transformImage[row*(int)(SAMPLING_FACTOR)][col*(int)(SAMPLING_FACTOR)].green = (CVD::byte)(mvEdgeDistanceTransformImage[row][col]*255.f);
         transformImage[row*(int)(SAMPLING_FACTOR)][col*(int)(SAMPLING_FACTOR)].blue = (CVD::byte)(mvEdgeDistanceTransformImage[row][col]*255.f);
         transformImage[row*(int)(SAMPLING_FACTOR)+1][col*(int)(SAMPLING_FACTOR)].red = (CVD::byte)(mvEdgeDistanceTransformImage[row][col]*255.f);
         transformImage[row*(int)(SAMPLING_FACTOR)+1][col*(int)(SAMPLING_FACTOR)].green = (CVD::byte)(mvEdgeDistanceTransformImage[row][col]*255.f);
         transformImage[row*(int)(SAMPLING_FACTOR)+1][col*(int)(SAMPLING_FACTOR)].blue = (CVD::byte)(mvEdgeDistanceTransformImage[row][col]*255.f);
         transformImage[row*(int)(SAMPLING_FACTOR)][col*(int)(SAMPLING_FACTOR)+1].red = (CVD::byte)(mvEdgeDistanceTransformImage[row][col]*255.f);
         transformImage[row*(int)(SAMPLING_FACTOR)][col*(int)(SAMPLING_FACTOR)+1].green = (CVD::byte)(mvEdgeDistanceTransformImage[row][col]*255.f);
         transformImage[row*(int)(SAMPLING_FACTOR)][col*(int)(SAMPLING_FACTOR)+1].blue = (CVD::byte)(mvEdgeDistanceTransformImage[row][col]*255.f);
         transformImage[row*(int)(SAMPLING_FACTOR)+1][col*(int)(SAMPLING_FACTOR)+1].red = (CVD::byte)(mvEdgeDistanceTransformImage[row][col]*255.f);
         transformImage[row*(int)(SAMPLING_FACTOR)+1][col*(int)(SAMPLING_FACTOR)+1].green = (CVD::byte)(mvEdgeDistanceTransformImage[row][col]*255.f);
         transformImage[row*(int)(SAMPLING_FACTOR)+1][col*(int)(SAMPLING_FACTOR)+1].blue = (CVD::byte)(mvEdgeDistanceTransformImage[row][col]*255.f);
      }
    }
}

void EdgeObjectDetection::getEdgeImage (CVD::Image<CVD::Rgb<CVD::byte> >& edgeImage)
{
    int width = mEdgeImage.size().x;
    int height = mEdgeImage.size().y;
    for( int row = 0; row < height; ++row )
    {
      for( int col = 0; col < width; ++col )
      {
         edgeImage[row*(int)(SAMPLING_FACTOR)][col*(int)(SAMPLING_FACTOR)].red = mEdgeImage[row][col];
         edgeImage[row*(int)(SAMPLING_FACTOR)][col*(int)(SAMPLING_FACTOR)].green = mEdgeImage[row][col];
         edgeImage[row*(int)(SAMPLING_FACTOR)][col*(int)(SAMPLING_FACTOR)].blue = mEdgeImage[row][col];
         edgeImage[row*(int)(SAMPLING_FACTOR)+1][col*(int)(SAMPLING_FACTOR)].red = mEdgeImage[row][col];
         edgeImage[row*(int)(SAMPLING_FACTOR)+1][col*(int)(SAMPLING_FACTOR)].green = mEdgeImage[row][col];
         edgeImage[row*(int)(SAMPLING_FACTOR)+1][col*(int)(SAMPLING_FACTOR)].blue = mEdgeImage[row][col];
         edgeImage[row*(int)(SAMPLING_FACTOR)][col*(int)(SAMPLING_FACTOR)+1].red = mEdgeImage[row][col];
         edgeImage[row*(int)(SAMPLING_FACTOR)][col*(int)(SAMPLING_FACTOR)+1].green = mEdgeImage[row][col];
         edgeImage[row*(int)(SAMPLING_FACTOR)][col*(int)(SAMPLING_FACTOR)+1].blue = mEdgeImage[row][col];
         edgeImage[row*(int)(SAMPLING_FACTOR)+1][col*(int)(SAMPLING_FACTOR)+1].red = mEdgeImage[row][col];
         edgeImage[row*(int)(SAMPLING_FACTOR)+1][col*(int)(SAMPLING_FACTOR)+1].green = mEdgeImage[row][col];
         edgeImage[row*(int)(SAMPLING_FACTOR)+1][col*(int)(SAMPLING_FACTOR)+1].blue = mEdgeImage[row][col];
      }
    }
}

void EdgeObjectDetection::preCalculateEdgeLnksWithRegions_MultiScale( LowLevelImageData<CVD::byte> & llImage, std::vector<SegmentedRegion> imageRegions, std::vector<SegmentedMask> imageMasks, int mfi)
{
  std::srand( time(NULL) );
  start_time = CVD::timer.get_time();
  std::vector<Line> vLines;
  master_file_id = mfi;
  double line_time;
  thisLlImage = llImage;
  // copy detected objects onto previous frame's detected objects
  mvDetectedObjects_previousFrame.clear();
  for (int i = 0; i < mvDetectedObjects.size(); i++)
  {
     mvDetectedObjects[i].missingCount = 0;
     mvDetectedObjects[i].irTopLeft.x -= 5;
     mvDetectedObjects[i].irTopLeft.y -= 5;
     mvDetectedObjects[i].irBottomRight.x += 5;
     mvDetectedObjects[i].irBottomRight.y += 5;
     mvDetectedObjects_previousFrame.push_back(mvDetectedObjects[i]);
  }
  for (int i = 0; i < mvDetectedObjects_undetected.size(); i++)
  {
       mvDetectedObjects_undetected[i].missingCount++;
       mvDetectedObjects_previousFrame.push_back (mvDetectedObjects_undetected[i]);
  }
  mvDetectedObjects_undetected.clear();
  mvDetectedObjects.clear();
  if (!isLSD)
  {
    //Edges edges;
    mCannyEdgeDetector.compute( mEdgeImage, edges, llImage, 8, 0.5, 0.3 ); //10, 0.5, 0.2
  }
  else
  {
    int no_data = 0;
    double * img_ptr = lsd_image->data;
    CVD::byte * int_img_ptr = llImage.mImage.data();
    while( no_data != data_size )
    {
      *img_ptr = static_cast<double>(*int_img_ptr);
      ++ no_data;
      ++ img_ptr;
      ++ int_img_ptr;
    }
    ntuple_list out;
    out = lsd(lsd_image);
    lsd_to_lines( vLines, out );
    lines_to_edges (edges, vLines);
    free_ntuple_list(out);
    line_time = CVD::timer.get_time();
  }

 // int maxValue = GetTransformDist::compute( mvEdgeDistanceTransformImage, mEdgeImage); 
  int maxValue = 1;
 // GetTransformDist::integral(mvEdgeDistanceTransformImageIntegral , mvEdgeDistanceTransformImage);
  GetTransformDist::integral(mvEdgeDistanceTransformImageIntegral , mEdgeImage);

  int edgeletNumber = (int)edges.virEdges.size();
  int edgelet_length = (ceil)((float)5 * ((float)(edgeletNumber) / (float)1000)) + 1;
//  int edgelet_length = EDGELET_LENGTH_TEST;
 // lineMaxLen = (ceil)((float)5 * ((float)(edgeletNumber) / (float)1000)) + 1;
  lineMaxLen = 4;

  if (!isLSD)
  {
    //get_lines( vLines, edges );
    get_lines_DT( vLines, edges, mvEdgeDistanceTransformImageIntegral );
    line_time = CVD::timer.get_time();
  }

 // EPS_DIST = 0.1 * (ceil)((float)edgelet_length/(float)5.f); //0.06 ICCV 0.15 without depth
  EPS_ANGLE  = 0.1 * (ceil)((float)edgelet_length/(float)5.f); //0.06 ICCV 0.15 without depth
 // std::cerr << "To be: " << EPS_DIST << " AND " << EPS_ANGLE << std::endl;
 for (int level = 0; level >= 0; level--)
  {
     gen_Edgelets_with_Directional_Regions( mvEdgelets, mviLineIndexes, vLines, imageRegions, maxValue, edgelet_length);
  //   std::cerr << "Number of edgelets found at level " << edgelet_length << " = " << mvEdgelets.size() << ", lines = " << lineMaxLen << std::endl;
     edgelet_length /= 2;
  }
  double edgelet_time = CVD::timer.get_time();
  gen_TestDirections( mvEdgelets, mviLineIndexes );
  #ifdef ECM
  gen_EdgeletCells (mvEdgelets);
  double beforeECM = CVD::timer.get_time ();
  EdgeCellMatrix* ecm =  getEdgeletImageMatrix (mvEdgelets, img_h, img_w, 10, 10);
  double afterECM = CVD::timer.get_time ();
  #endif
  double td_time = CVD::timer.get_time();
}

void EdgeObjectDetection::preCalculateEdgeLnksWithRegions_depth( LowLevelImageData<CVD::byte> & llImage, CVD::Image<CVD::byte> & depth_image, std::vector<SegmentedRegion> imageRegions, int mfi)
{
  std::srand( time(NULL) );
  start_time = CVD::timer.get_time();
  std::vector<Line> vLines;
  master_file_id = mfi;
  double line_time;
  thisLlImage = llImage;
  // copy detected objects onto previous frame's detected objects
  mvDetectedObjects_previousFrame.clear();
  for (int i = 0; i < mvDetectedObjects.size(); i++)
  {
     mvDetectedObjects[i].missingCount = 0;
     mvDetectedObjects[i].irTopLeft.x -= 5;
     mvDetectedObjects[i].irTopLeft.y -= 5;
     mvDetectedObjects[i].irBottomRight.x += 5;
     mvDetectedObjects[i].irBottomRight.y += 5;
     mvDetectedObjects_previousFrame.push_back(mvDetectedObjects[i]);
  }
  for (int i = 0; i < mvDetectedObjects_undetected.size(); i++)
  {
       mvDetectedObjects_undetected[i].missingCount++;
       mvDetectedObjects_previousFrame.push_back (mvDetectedObjects_undetected[i]);
    //}
  }
  mvDetectedObjects_undetected.clear();
  mvDetectedObjects.clear();
  if (!isLSD)
  {
    //Edges edges;
       mCannyEdgeDetector.compute( mEdgeImage, edges, llImage, 10, 0.5, 0.2 );
    get_lines( vLines, edges );
//    std::cout << edges.size() << " " << vLines.size() << std::endl;
    line_time = CVD::timer.get_time();
    char filename[100];
  //  sprintf(filename,"%s/master_%d.xvy", folderName, master_file_id);
  //  std::ofstream out_xvy(filename);
  //  for( unsigned int i = 0; i < edges.virEdges.size(); ++i )
  //  {
  //    mEdgeImage[edges.virEdges[i]] = 255;
  //    out_xvy << CVD::vec(edges.virEdges[i]) << std::endl;
  //  }
  //  out_xvy.close();
  //  sprintf(filename,"%s/master_%d.png",folderName, master_file_id);
  //  //++master_file_id;
  //  CVD::img_save(mEdgeImage, filename);
  }
  else
  {
    int no_data = 0;
    double * img_ptr = lsd_image->data;
    CVD::byte * int_img_ptr = llImage.mImage.data();
    while( no_data != data_size )
    {
      *img_ptr = static_cast<double>(*int_img_ptr);
      ++ no_data;
      ++ img_ptr;
      ++ int_img_ptr;
    }
    ntuple_list out;
    out = lsd(lsd_image);
    lsd_to_lines( vLines, out );
    lines_to_edges (edges, vLines);
    free_ntuple_list(out);
    line_time = CVD::timer.get_time();
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
      std::cout << "CANNY " << filename << std::endl;
      CVD::img_save(mEdgeImage, filename);
    #endif
  }
  // calculate the overlap between image regions
  
  gen_Edgelets_with_Directional_Regions( mvEdgelets, mviLineIndexes, vLines, imageRegions, 0, EDGELET_LENGTH_TEST);
  double edgelet_time = CVD::timer.get_time();
  gen_TestDirections( mvEdgelets, mviLineIndexes );
  #ifdef ECM
  gen_EdgeletCells (mvEdgelets);
  double beforeECM = CVD::timer.get_time ();
  EdgeCellMatrix* ecm =  getEdgeletImageMatrix (mvEdgelets, img_h, img_w, 10, 10);
  double afterECM = CVD::timer.get_time ();
  #endif
  double td_time = CVD::timer.get_time();
}

void EdgeObjectDetection::preCalculateEdgeLnks( LowLevelImageData<CVD::byte> & llImage, int mfi)
{
  std::srand( time(NULL) );
  start_time = CVD::timer.get_time();
  std::vector<Line> vLines;
  master_file_id = mfi;
  double line_time;
  // copy detected objects onto previous frame's detected objects
  mvDetectedObjects_previousFrame.clear();
  for (int i = 0; i < mvDetectedObjects_undetected.size(); i++)
  {
    mvDetectedObjects_undetected[i].missingCount++;
    mvDetectedObjects_previousFrame.push_back (mvDetectedObjects_undetected[i]);
  }
  mvDetectedObjects_undetected.clear();
  for (int i = 0; i < mvDetectedObjects.size(); i++)
  {
     mvDetectedObjects[i].missingCount = 0;
     mvDetectedObjects_previousFrame.push_back(mvDetectedObjects[i]);
  }
  mvDetectedObjects.clear();
//  std::cout << "NUMBER OF PREVIOUS FRAMES = " << mvDetectedObjects_previousFrame.size() << std::endl;
  //std::cout << " rewriting master_file_id to " << master_file_id << std::endl;
if (!isLSD)
{
  //Edges edges;
  mCannyEdgeDetector.compute( mEdgeImage, edges, llImage, 10, 0.5, 0.2 );
//  for( int i = 0; i < edges.virEdges.size(); ++i )
//    std::cout <<  CVD::vec(edges.virEdges[i]) << std::endl;
  get_lines( vLines, edges );
  line_time = CVD::timer.get_time();
  mEdgeImage.fill(0);
//  char filename[100];
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
  CVD::img_save(mEdgeImage, filename);*/
}
else
{
  int no_data = 0;
  double * img_ptr = lsd_image->data;
  CVD::byte * int_img_ptr = llImage.mImage.data();
  while( no_data != data_size )
  {
    *img_ptr = static_cast<double>(*int_img_ptr);
    ++ no_data;
    ++ img_ptr;
    ++ int_img_ptr;
  }
  ntuple_list out;
  out = lsd(lsd_image);
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
  gen_Edgelets_with_Directional( mvEdgelets, mviLineIndexes, vLines, EDGELET_LENGTH_TEST );
//  gen_TestDirections_for_detection( mvEdgelets, mviLineIndexes ); // switching between these two lines need to change mvTestDirections and mlTestDirections in detect().
  double edgelet_time = CVD::timer.get_time();
  gen_TestDirections( mvEdgelets, mviLineIndexes );
  #ifdef ECM
  gen_EdgeletCells (mvEdgelets);
  double beforeECM = CVD::timer.get_time ();
  EdgeCellMatrix* ecm =  getEdgeletImageMatrix (mvEdgelets, img_h, img_w, 10, 10);
  double afterECM = CVD::timer.get_time ();
  #endif
  double td_time = CVD::timer.get_time();
//  std::cout << "Number of edgelets = " << mvTestDirections.size() << std::endl;

//  std::cout << "PreCalculate time = Line : " << line_time - start_time<< ", Edgelet : " << edgelet_time - line_time << ", TestDirection : " << td_time - edgelet_time << std::endl;
}

void EdgeObjectDetection::preCalculateEdgeLnks_depth( LowLevelImageData<CVD::byte> & llImage, CVD::Image<CVD::byte> & depth_image, int mfi)
{
  std::srand( time(NULL) );
  start_time = CVD::timer.get_time();
  std::vector<Line> vLines;
  master_file_id = mfi;
  double line_time;
  // copy detected objects onto previous frame's detected objects
  mvDetectedObjects_previousFrame.clear();
  for (int i = 0; i < mvDetectedObjects_undetected.size(); i++)
  {
    mvDetectedObjects_undetected[i].missingCount++;
    mvDetectedObjects_previousFrame.push_back (mvDetectedObjects_undetected[i]);
  }
  mvDetectedObjects_undetected.clear();
  for (int i = 0; i < mvDetectedObjects.size(); i++)
  {
     mvDetectedObjects[i].missingCount = 0;
     mvDetectedObjects_previousFrame.push_back(mvDetectedObjects[i]);
  }
  mvDetectedObjects.clear();
//  std::cout << "NUMBER OF PREVIOUS FRAMES = " << mvDetectedObjects_previousFrame.size() << std::endl;
  //std::cout << " rewriting master_file_id to " << master_file_id << std::endl;
if (!isLSD)
{
  //Edges edges;
     mCannyEdgeDetector.compute( mEdgeImage, edges, llImage, 10, 0.5, 0.2 );
//  for( int i = 0; i < edges.virEdges.size(); ++i )
//    std::cout <<  CVD::vec(edges.virEdges[i]) << std::endl;
  get_lines( vLines, edges );
  line_time = CVD::timer.get_time();
  mEdgeImage.fill(0);
//  char filename[100];
// SAVING MASTER FILE
//  sprintf(filename,"%s/master_%d.xvy", folderName, master_file_id);
//  std::ofstream out_xvy(filename);
//  for( unsigned int i = 0; i < edges.virEdges.size(); ++i )
//  {
//    mEdgeImage[edges.virEdges[i]] = 255;
//    out_xvy << CVD::vec(edges.virEdges[i]) << std::endl;
//  }
//  out_xvy.close();
//  sprintf(filename,"%s/master_%d.png",folderName, master_file_id);
  //++master_file_id;
//  CVD::img_save(mEdgeImage, filename);
}
else
{
  int no_data = 0;
  double * img_ptr = lsd_image->data;
  CVD::byte * int_img_ptr = llImage.mImage.data();
  while( no_data != data_size )
  {
    *img_ptr = static_cast<double>(*int_img_ptr);
    ++ no_data;
    ++ img_ptr;
    ++ int_img_ptr;
  }
  ntuple_list out;
  out = lsd(lsd_image);
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
  gen_Edgelets_with_Directional( mvEdgelets, mviLineIndexes, vLines, EDGELET_LENGTH_TEST );
//  gen_TestDirections_for_detection( mvEdgelets, mviLineIndexes ); // switching between these two lines need to change mvTestDirections and mlTestDirections in detect().
  double edgelet_time = CVD::timer.get_time();
  gen_TestDirections( mvEdgelets, mviLineIndexes );
  #ifdef ECM
  gen_EdgeletCells (mvEdgelets);
  double beforeECM = CVD::timer.get_time ();
  EdgeCellMatrix* ecm =  getEdgeletImageMatrix (mvEdgelets, img_h, img_w, 10, 10);
  double afterECM = CVD::timer.get_time ();
  #endif
  double td_time = CVD::timer.get_time();
//  std::cout << "Number of edgelets = " << mvTestDirections.size() << std::endl;

//  std::cout << "PreCalculate time = Line : " << line_time - start_time<< ", Edgelet : " << edgelet_time - line_time << ", TestDirection : " << td_time - edgelet_time << std::endl;
}

void EdgeObjectDetection::getLSDEdges (Edges & thisEdges)
{
   for (int i = 0; i < edges.virEdges.size(); i++){
      edges.virEdges[i].x *= SAMPLING_FACTOR;
      edges.virEdges[i].y *= SAMPLING_FACTOR;
      thisEdges.virEdges.push_back (edges.virEdges[i]);
   }
}

void EdgeObjectDetection::getEdges (Edges & thisEdges)
{
   for (int i = 0; i < edges.virEdges.size(); i++){
      edges.virEdges[i].x *= SAMPLING_FACTOR;
      edges.virEdges[i].y *= SAMPLING_FACTOR;
      thisEdges.virEdges.push_back (edges.virEdges[i]);
   }
}

void EdgeObjectDetection::getLSDEdgesInRegions (Edges & thisEdges, std::vector<SegmentedRegion> imageRegions)
{
   for (int i = 0; i < edges.virEdges.size(); i++)
   {
      edges.virEdges[i].x *= 4;
      edges.virEdges[i].y *= 4;
      for (int r = 0; r < imageRegions.size(); r++)
          if (edges.virEdges[i].x >= imageRegions[r].irTopLeft.x && edges.virEdges[i].x <= imageRegions[r].irBottomRight.x && edges.virEdges[i].y >= imageRegions[r].irTopLeft.y && edges.virEdges[i].y <= imageRegions[r].irBottomRight.y)
          {
             thisEdges.virEdges.push_back (edges.virEdges[i]);
             break;
          }
   }
}

void EdgeObjectDetection::getFoundChains (std::vector<CheckedChains> & thisChecked)
{
   for (int i = 0; i < mvCheckedChains.size(); i++){
      thisChecked.push_back (mvCheckedChains[i]);
   }
}

void EdgeObjectDetection::getMatchedConstellations (std::vector< TooN::Vector<2> > & thisChecked)
{
    for (int i = 0; i < mvCheckedChains.size(); i++){
      for (int j = 0; j < mvCheckedChains[i].chains.size(); j++)
         thisChecked.push_back (mvCheckedChains[i].chains[j]);
   }
}

void EdgeObjectDetection::get_lines( std::vector<Line> & vLines, const Edges & edges )
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
      if( line.no_pixel >  (float)(lineMaxLen))
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

void EdgeObjectDetection::get_lines_DT( std::vector<Line> & vLines, const Edges & edges , const CVD::Image<REAL_TYPE> & mvEdgeDistanceTransformImage )
{
  CVD::ImageRef window_dt1, window_dt2, window_dt3;
  const int window_size = 15;
  window_dt1.x = window_dt2.x = window_size;
  window_dt1.y = window_dt3.y = window_size;
  window_dt3.x = -1*window_size;
  window_dt2.y = -1*window_size;
  CVD::ImageRef imageSize = mvEdgeDistanceTransformImage.size();
  CVD::ImageRef top_left, top_right, bottom_left, bottom_right;
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
      REAL_TYPE dt_sum = 0.f;
      int dt_count = 0;
      
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
          //std::cerr << "Can read distance transform at the point to be " << (REAL_TYPE)(mvEdgeDistanceTransformImage[irPoint]) << std::endl;
          // here we need to calculate the dimensions of the DT around this point
          top_left = irPoint - window_dt1;
          top_left.x = ( top_left.x < 0 ) ? 0 : top_left.x;
          top_left.y = ( top_left.y < 0 ) ? 0 : top_left.y;
          top_right = irPoint + window_dt2;
          top_right.x = ( top_right.x >=  imageSize.x) ? imageSize.x-1 : top_right.x;
          top_right.y = ( top_right.y < 0 ) ? 0 : top_right.y;
          bottom_left = irPoint + window_dt3;
          bottom_left.x = ( bottom_left.x < 0 ) ? 0 : bottom_left.x;
          bottom_left.y = ( bottom_left.y >=  imageSize.y) ? imageSize.y-1 : bottom_left.y;
          bottom_right = irPoint + window_dt1;
          bottom_right.x = ( bottom_right.x >=  imageSize.x) ? imageSize.x-1 : bottom_right.x;
          bottom_right.y = ( bottom_right.y >=  imageSize.y) ? imageSize.y-1 : bottom_right.y;

          // calculate and average distance transform using integral image
          REAL_TYPE integral_cal = mvEdgeDistanceTransformImage[bottom_right]-mvEdgeDistanceTransformImage[bottom_left]-mvEdgeDistanceTransformImage[top_right]+mvEdgeDistanceTransformImage[top_left];
          integral_cal /= (REAL_TYPE)(window_size*window_size*4*128);
          dt_sum += integral_cal;
         // if (integral_cal > 1)
         //    std::cerr << "***: " << integral_cal << " " << top_left << " " << bottom_left << " " << top_right << " " << bottom_right << std::endl;
          dt_count++;
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
      if( line.no_pixel >  (float)(lineMaxLen))
      {
        line.distanceTransformAverage = dt_sum/((REAL_TYPE) dt_count);
//        std::cerr << "Dist Trans Avg = " << line.distanceTransformAverage << std::endl;
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

void EdgeObjectDetection::gen_Edgelets( std::vector<Edgelet> & vEdgelets, std::vector<int> & viLineIndexes, CVD::ImageRef irTopLeft, CVD::ImageRef irButtomRight, const std::vector<Line> & vLines, int iMaxLength )
{
  vEdgelets.clear();
  REAL_TYPE halfMaxLength = static_cast<REAL_TYPE>(iMaxLength)*0.5f;
  for( unsigned int i = 0; i < vLines.size(); ++i )
  {
    const Line & line = vLines[i];
    int no_segment = line.no_pixel/iMaxLength;// iMaxLength pixel per edgelet
    viLineIndexes.push_back( vEdgelets.size() );
    if( no_segment == 0 )
    {
      if( line.v2Start[0] < 0 || line.v2Start[1] < 0 || line.v2End[0] < 0 || line.v2End[1] < 0 )
        continue;      
      Edgelet edgeLet;
      edgeLet.v2Pose = (line.v2Start + line.v2End)*0.5f;
      edgeLet.v2Pose[0] = round(edgeLet.v2Pose[0] ); edgeLet.v2Pose[1] = round(edgeLet.v2Pose[1]);
      if( edgeLet.v2Pose[0] < irTopLeft.x || edgeLet.v2Pose[1] < irTopLeft.y || edgeLet.v2Pose[0] > irButtomRight.x || edgeLet.v2Pose[1] > irButtomRight.y )
        continue;
      edgeLet.previousObjNo = -1;
//      edgeLet.rAngle = atan(line.rSlope);
      edgeLet.rSlope = line.rSlope;
      vEdgelets.push_back(edgeLet);
    }
    else
    {
//      REAL_TYPE rAngle = atan(line.rSlope);
      for( int j = 0; j < no_segment; ++ j )
      {
        Edgelet edgeLet;
        edgeLet.v2Pose = line.v2Start + (iMaxLength*j+halfMaxLength)*line.v2Direction;
        edgeLet.v2Pose[0] = round(edgeLet.v2Pose[0] ); edgeLet.v2Pose[1] = round(edgeLet.v2Pose[1]);
        if( edgeLet.v2Pose[0] < irTopLeft.x || edgeLet.v2Pose[1] < irTopLeft.y || edgeLet.v2Pose[0] > irButtomRight.x || edgeLet.v2Pose[1] > irButtomRight.y )
          continue;
//        edgeLet.rAngle = rAngle;
        edgeLet.rSlope = line.rSlope;
        edgeLet.previousObjNo = -1;
        vEdgelets.push_back(edgeLet);
      }
      if( vEdgelets.size() > 0 && norm(vEdgelets.back().v2Pose - line.v2End) > 1 )
      {
        if( line.v2End[0] < 0 || line.v2End[1] < 0 )
          continue;
        Edgelet edgeLet;
        edgeLet.v2Pose = ( line.v2Start + iMaxLength*no_segment*line.v2Direction + line.v2End )*0.5f;
        edgeLet.v2Pose[0] = round(edgeLet.v2Pose[0] ); edgeLet.v2Pose[1] = round(edgeLet.v2Pose[1]);
        if( edgeLet.v2Pose[0] < irTopLeft.x || edgeLet.v2Pose[1] < irTopLeft.y || edgeLet.v2Pose[0] > irButtomRight.x || edgeLet.v2Pose[1] > irButtomRight.y )
          continue;
//        edgeLet.rAngle = rAngle;
        edgeLet.rSlope = line.rSlope;
        edgeLet.previousObjNo = -1;
        vEdgelets.push_back(edgeLet);
      }
    }
  }
  viLineIndexes.push_back( vEdgelets.size() );
}

void EdgeObjectDetection::gen_Edgelets( std::vector<Edgelet> & vEdgelets, std::vector<int> & viLineIndexes, const std::vector<Line> & vLines, int iMaxLength )
{
  vEdgelets.clear();
  REAL_TYPE halfMaxLength = static_cast<REAL_TYPE>(iMaxLength)*0.5f;
  for( unsigned int i = 0; i < vLines.size(); ++i )
  {
    const Line & line = vLines[i];
    int no_segment = line.no_pixel/iMaxLength;// iMaxLength pixel per edgelet
    viLineIndexes.push_back( vEdgelets.size() );
    if( no_segment == 0 )
    {
      if( line.v2Start[0] < 0 || line.v2Start[1] < 0 || line.v2End[0] < 0 || line.v2End[1] < 0 )
        continue;      
      Edgelet edgeLet;
      edgeLet.v2Pose = (line.v2Start + line.v2End)*0.5f;
      edgeLet.v2Pose[0] = round(edgeLet.v2Pose[0] ); edgeLet.v2Pose[1] = round(edgeLet.v2Pose[1]);
//      edgeLet.rAngle = atan(line.rSlope);
      edgeLet.rSlope = line.rSlope;
      edgeLet.previousObjNo = -1;
      vEdgelets.push_back(edgeLet);
    }
    else
    {
//      REAL_TYPE rAngle = atan(line.rSlope);
      for( int j = 0; j < no_segment; ++ j )
      {
        Edgelet edgeLet;
        edgeLet.v2Pose = line.v2Start + (iMaxLength*j+halfMaxLength)*line.v2Direction;
        edgeLet.v2Pose[0] = round(edgeLet.v2Pose[0] ); edgeLet.v2Pose[1] = round(edgeLet.v2Pose[1]);
//        edgeLet.rAngle = rAngle;
        edgeLet.rSlope = line.rSlope;
        edgeLet.previousObjNo = -1;
        vEdgelets.push_back(edgeLet);
      }
      if( norm(vEdgelets.back().v2Pose - line.v2End) > 1 )
      {
        if( line.v2End[0] < 0 || line.v2End[1] < 0 )
          continue;
        Edgelet edgeLet;
        edgeLet.v2Pose = ( line.v2Start + iMaxLength*no_segment*line.v2Direction + line.v2End )*0.5f;
        edgeLet.v2Pose[0] = round(edgeLet.v2Pose[0] ); edgeLet.v2Pose[1] = round(edgeLet.v2Pose[1]);
//        edgeLet.rAngle = rAngle;
        edgeLet.rSlope = line.rSlope;
        edgeLet.previousObjNo = -1;
        vEdgelets.push_back(edgeLet);
      }
    }
  }
  viLineIndexes.push_back( vEdgelets.size() );
}

void EdgeObjectDetection::gen_Edgelets_with_Directional_Masks( std::vector<Edgelet> & vEdgelets, std::vector<int> & viLineIndexes, const std::vector<Line> & vLines, std::vector<SegmentedMask> imageMasks, long oriSize, int iMaxLength)
{
  double local_start_time = CVD::timer.get_time();
  bool pushed = false;
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
      pushed = false;
     /* for (int p = 0; p < mvDetectedObjects_previousFrame.size(); p++)
         if (edgeLet.v2Pose[0] >= (mvDetectedObjects_previousFrame[p].irTopLeft.x) && edgeLet.v2Pose[0] <= (mvDetectedObjects_previousFrame[p].irBottomRight.x) && edgeLet.v2Pose[1] >= (mvDetectedObjects_previousFrame[p].irTopLeft.y) && edgeLet.v2Pose[1] <= (mvDetectedObjects_previousFrame[p].irBottomRight.y))
         {
            edgeLet.previousObjNo = p;
            vEdgelets.push_back(edgeLet);
            pushed = true;
            break;
         }  
      if (!pushed)
         edgeLet.previousObjNo = -1;       */
      for (int r = 0; r < imageMasks.size(); r++) 
         if (imageMasks[r].image.data[2*edgeLet.v2Pose[1]*imageMasks[r].image.width+2*edgeLet.v2Pose[0]] > 0)
      //   if (edgeLet.v2Pose[0] >= (imageMasks[r].irTopLeft.x/2) && edgeLet.v2Pose[0] <= (imageRegions[r].irBottomRight.x/2) && edgeLet.v2Pose[1] >= (imageRegions[r].irTopLeft.y/2) && edgeLet.v2Pose[1] <= (imageRegions[r].irBottomRight.y/2))
         {
            edgeLet.regionNo = r;
            if (!pushed)
               vEdgelets.push_back(edgeLet);
            break;
         }
    }
    else
    {
      bool addedOne = false;
      for( int j = 0; j < no_segment; ++ j )
      {
        Edgelet edgeLet;
        edgeLet.v2Pose = line.v2Start + (iMaxLength*j+halfMaxLength)*line.v2Direction;
        edgeLet.v2Pose[0] = round(edgeLet.v2Pose[0] ); edgeLet.v2Pose[1] = round(edgeLet.v2Pose[1]);

        edgeLet.rSlope = line.rSlope;
        pushed = false;
  /*      for (int p = 0; p < mvDetectedObjects_previousFrame.size(); p++)
           if (edgeLet.v2Pose[0] >= (mvDetectedObjects_previousFrame[p].irTopLeft.x) && edgeLet.v2Pose[0] <= (mvDetectedObjects_previousFrame[p].irBottomRight.x) && edgeLet.v2Pose[1] >= (mvDetectedObjects_previousFrame[p].irTopLeft.y) && edgeLet.v2Pose[1] <= (mvDetectedObjects_previousFrame[p].irBottomRight.y))
           {
              edgeLet.previousObjNo = p;
              vEdgelets.push_back(edgeLet);
              pushed = true;
              break;
           }         
        if (!pushed)
           edgeLet.previousObjNo = -1;       */
        for (int r = 0; r < imageMasks.size(); r++)
        {
           if (imageMasks[r].image.data[2*edgeLet.v2Pose[1]*imageMasks[r].image.width+2*edgeLet.v2Pose[0]] > 0)
//           if (edgeLet.v2Pose[0] >= (imageRegions[r].irTopLeft.x/2) && edgeLet.v2Pose[0] <= (imageRegions[r].irBottomRight.x/2) && edgeLet.v2Pose[1] >= (imageRegions[r].irTopLeft.y/2) && edgeLet.v2Pose[1] <= (imageRegions[r].irBottomRight.y/2))
           {
              edgeLet.regionNo = r;
              if (!pushed)
                 vEdgelets.push_back(edgeLet);
              addedOne = true;
              break;
           }
         }
      }
      if(addedOne && norm(vEdgelets.back().v2Pose - line.v2End) > 1 )
      {
        if( line.v2End[0] < 0 || line.v2End[1] < 0 )
          continue;
        Edgelet edgeLet;
        edgeLet.v2Pose = ( line.v2Start + iMaxLength*no_segment*line.v2Direction + line.v2End )*0.5f;
        edgeLet.v2Pose[0] = round(edgeLet.v2Pose[0] ); edgeLet.v2Pose[1] = round(edgeLet.v2Pose[1]);
        pushed = false;
/*        for (int p = 0; p < mvDetectedObjects_previousFrame.size(); p++)
           if (edgeLet.v2Pose[0] >= (mvDetectedObjects_previousFrame[p].irTopLeft.x) && edgeLet.v2Pose[0] <= (mvDetectedObjects_previousFrame[p].irBottomRight.x) && edgeLet.v2Pose[1] >= (mvDetectedObjects_previousFrame[p].irTopLeft.y) && edgeLet.v2Pose[1] <= (mvDetectedObjects_previousFrame[p].irBottomRight.y))
           {
              edgeLet.previousObjNo = p;
              vEdgelets.push_back(edgeLet);
              pushed = true;
              break;
           }         
        if (!pushed)
           edgeLet.previousObjNo = -1;       */

        if (edgeLet.v2Pose[0] > 0 && edgeLet.v2Pose[0] <= img_h && edgeLet.v2Pose[1] > 0 && edgeLet.v2Pose[1] <= img_w)
        {
           edgeLet.rSlope = line.rSlope;
           for (int r = 0; r < imageMasks.size(); r++)
              if (imageMasks[r].image.data[2*edgeLet.v2Pose[1]*imageMasks[r].image.width+2*edgeLet.v2Pose[0]] > 0)
   //           if (edgeLet.v2Pose[0] >= (imageRegions[r].irTopLeft.x/2) && edgeLet.v2Pose[0] <= (imageRegions[r].irBottomRight.x/2) && edgeLet.v2Pose[1] >= (imageRegions[r].irTopLeft.y/2) && edgeLet.v2Pose[1] <= (imageRegions[r].irBottomRight.y/2))
              {
                 edgeLet.regionNo = r;
                 if (!pushed)
                    vEdgelets.push_back(edgeLet);
                break;
              }
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
  double edge_ori_time = CVD::timer.get_time();
}

void EdgeObjectDetection::gen_Edgelets_with_Directional_Regions( std::vector<Edgelet> & vEdgelets, std::vector<int> & viLineIndexes, const std::vector<Line> & vLines, std::vector<SegmentedRegion> imageRegions, int maxValue, int iMaxLength, long oriSize)
{
  
  double local_start_time = CVD::timer.get_time();
  bool pushed = false;
  vEdgelets.clear();
  viLineIndexes.clear();
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
    
    // here is the part where we can replace this by the information we have from the distance Transform
    // iMaxLength pixel per edgelet
    // if we wish to calculate from distance transform
  //  int edited_max_length = (int)(log((line.distanceTransformAverage/(REAL_TYPE)maxValue))*0.0953f+30);
    int edited_max_length = (int)((line.distanceTransformAverage*100)*(line.distanceTransformAverage*100)-60);
    if (edited_max_length < 3)
       edited_max_length = 3;
    if (edited_max_length > 40)
       edited_max_length = 40;
    if (edited_max_length > line.no_pixel)
       edited_max_length = line.no_pixel;
   // iMaxLength = edited_max_length;
    REAL_TYPE halfMaxLength = static_cast<REAL_TYPE>(iMaxLength)*0.5f;

    int no_segment = line.no_pixel/iMaxLength;

    //std::cerr << line.distanceTransformAverage << " " << edited_max_length << std::endl;
    viLineIndexes.push_back( vEdgelets.size() );
    if( no_segment == 0 )
    {
      if( line.v2Start[0] < 0 || line.v2Start[1] < 0 || line.v2End[0] < 0 || line.v2End[1] < 0 )
        continue;
      Edgelet edgeLet;
      edgeLet.v2Pose = (line.v2Start + line.v2End)*0.5;
      edgeLet.v2Pose[0] = round(edgeLet.v2Pose[0] ); edgeLet.v2Pose[1] = round(edgeLet.v2Pose[1]);
      edgeLet.rSlope = line.rSlope;
      pushed = false;
      for (int p = 0; p < mvDetectedObjects_previousFrame.size(); p++)
         if (edgeLet.v2Pose[0] >= (mvDetectedObjects_previousFrame[p].irTopLeft.x) && edgeLet.v2Pose[0] <= (mvDetectedObjects_previousFrame[p].irBottomRight.x) && edgeLet.v2Pose[1] >= (mvDetectedObjects_previousFrame[p].irTopLeft.y) && edgeLet.v2Pose[1] <= (mvDetectedObjects_previousFrame[p].irBottomRight.y))
         {
            edgeLet.previousObjNo = p;
            vEdgelets.push_back(edgeLet);
            pushed = true;
            break;
         }  
      if (!pushed)
         edgeLet.previousObjNo = -1;       
      for (int r = 0; r < imageRegions.size(); r++)
         if (edgeLet.v2Pose[0] >= (imageRegions[r].irTopLeft.x/2) && edgeLet.v2Pose[0] <= (imageRegions[r].irBottomRight.x/2) && edgeLet.v2Pose[1] >= (imageRegions[r].irTopLeft.y/2) && edgeLet.v2Pose[1] <= (imageRegions[r].irBottomRight.y/2))
         {
            thisIrTopLeft = imageRegions[r].irTopLeft;
            thisIrBottomRight = imageRegions[r].irBottomRight;
            edgeLet.regionNo = r;
            if (!pushed)
               vEdgelets.push_back(edgeLet);
            break;
         }
    }
    else
    {
      bool addedOne = false;
      for( int j = 0; j < no_segment; ++ j )
      {
        Edgelet edgeLet;
        edgeLet.v2Pose = line.v2Start + (iMaxLength*j+halfMaxLength)*line.v2Direction;
        edgeLet.v2Pose[0] = round(edgeLet.v2Pose[0] ); edgeLet.v2Pose[1] = round(edgeLet.v2Pose[1]);

        edgeLet.rSlope = line.rSlope;
        pushed = false;
        for (int p = 0; p < mvDetectedObjects_previousFrame.size(); p++)
           if (edgeLet.v2Pose[0] >= (mvDetectedObjects_previousFrame[p].irTopLeft.x) && edgeLet.v2Pose[0] <= (mvDetectedObjects_previousFrame[p].irBottomRight.x) && edgeLet.v2Pose[1] >= (mvDetectedObjects_previousFrame[p].irTopLeft.y) && edgeLet.v2Pose[1] <= (mvDetectedObjects_previousFrame[p].irBottomRight.y))
           {
              edgeLet.previousObjNo = p;
              vEdgelets.push_back(edgeLet);
              pushed = true;
              break;
           }         
        if (!pushed)
           edgeLet.previousObjNo = -1;       
        for (int r = 0; r < imageRegions.size(); r++)
        {
           if (edgeLet.v2Pose[0] >= (imageRegions[r].irTopLeft.x/2) && edgeLet.v2Pose[0] <= (imageRegions[r].irBottomRight.x/2) && edgeLet.v2Pose[1] >= (imageRegions[r].irTopLeft.y/2) && edgeLet.v2Pose[1] <= (imageRegions[r].irBottomRight.y/2))
           {
              edgeLet.regionNo = r;
              thisIrTopLeft = imageRegions[r].irTopLeft;
              thisIrBottomRight = imageRegions[r].irBottomRight;
              if (!pushed)
                 vEdgelets.push_back(edgeLet);
              addedOne = true;
              break;
           }
         }
      }
      if(addedOne && norm(vEdgelets.back().v2Pose - line.v2End) > 1 )
      {
        if( line.v2End[0] < 0 || line.v2End[1] < 0 )
          continue;
        Edgelet edgeLet;
        edgeLet.v2Pose = ( line.v2Start + iMaxLength*no_segment*line.v2Direction + line.v2End )*0.5f;
        edgeLet.v2Pose[0] = round(edgeLet.v2Pose[0] ); edgeLet.v2Pose[1] = round(edgeLet.v2Pose[1]);
        pushed = false;
        for (int p = 0; p < mvDetectedObjects_previousFrame.size(); p++)
           if (edgeLet.v2Pose[0] >= (mvDetectedObjects_previousFrame[p].irTopLeft.x) && edgeLet.v2Pose[0] <= (mvDetectedObjects_previousFrame[p].irBottomRight.x) && edgeLet.v2Pose[1] >= (mvDetectedObjects_previousFrame[p].irTopLeft.y) && edgeLet.v2Pose[1] <= (mvDetectedObjects_previousFrame[p].irBottomRight.y))
           {
              edgeLet.previousObjNo = p;
              vEdgelets.push_back(edgeLet);
              pushed = true;
              break;
           }         
        if (!pushed)
           edgeLet.previousObjNo = -1;       

        if (edgeLet.v2Pose[0] > 0 && edgeLet.v2Pose[0] <= img_h && edgeLet.v2Pose[1] > 0 && edgeLet.v2Pose[1] <= img_w)
        {
           edgeLet.rSlope = line.rSlope;
           for (int r = 0; r < imageRegions.size(); r++)
              if (edgeLet.v2Pose[0] >= (imageRegions[r].irTopLeft.x/2) && edgeLet.v2Pose[0] <= (imageRegions[r].irBottomRight.x/2) && edgeLet.v2Pose[1] >= (imageRegions[r].irTopLeft.y/2) && edgeLet.v2Pose[1] <= (imageRegions[r].irBottomRight.y/2))
              {
                 edgeLet.regionNo = r;
                 if (!pushed)
                    vEdgelets.push_back(edgeLet);
                 break;
              }
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
  double edge_ori_time = CVD::timer.get_time();
}

void EdgeObjectDetection::gen_Edgelets_with_Directional( std::vector<Edgelet> & vEdgelets, std::vector<int> & viLineIndexes, const std::vector<Line> & vLines, long oriSize, int iMaxLength)
{
  double local_start_time = CVD::timer.get_time();
  vEdgelets.clear();
  viLineIndexes.clear();
  bool pushed = false;
  for( int i = 0; i < 11; ++ i )
  {
//    mvEdgeOrientationImages[i].fill(1);
    mvEdgeOrientations[i].reset();
  }
  REAL_TYPE halfMaxLength = static_cast<REAL_TYPE>(iMaxLength)*0.5f;
//  std::cout << "Lines ..." << std::endl;
  for( unsigned int i = 0; i < vLines.size(); ++i )
  {
    const Line & line = vLines[i];
//    std::cout << line.v2Start << " " << line.v2End << std::endl;
    int binNo = static_cast<int>( (atan( line.rSlope )/M_PI + 0.5f )*10.f );
//    mvEdgeOrientationImages[binNo][CVD::ir(line.v2End)] = 0;
    mvEdgeOrientations[binNo][line.v2End] = 1;    
    for( int no_pixel = 0; no_pixel <= line.no_pixel; ++no_pixel )
    {
      TooN::Vector<2> v2Pose = line.v2Start + no_pixel*line.v2Direction;
//      mvEdgeOrientationImages[binNo][CVD::ir(v2Pose)] = 0;
        mvEdgeOrientations[binNo][v2Pose] = 1;      
//      if( binNo == 0 )
//        std::cout << CVD::ir(v2Pose) << " " << v2Pose << " " << line.v2Start << " " << line.v2End << " " << line.no_pixel << std::endl;
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

//      edgeLet.rAngle = atan(line.rSlope);
      edgeLet.rSlope = line.rSlope;
      edgeLet.regionNo = 0;
      pushed = false;
      for (int p = 0; p < mvDetectedObjects_previousFrame.size(); p++)
        if (edgeLet.v2Pose[0] >= (mvDetectedObjects_previousFrame[p].irTopLeft.x) && edgeLet.v2Pose[0] <= (mvDetectedObjects_previousFrame[p].irBottomRight.x) && edgeLet.v2Pose[1] >= (mvDetectedObjects_previousFrame[p].irTopLeft.y) && edgeLet.v2Pose[1] <= (mvDetectedObjects_previousFrame[p].irBottomRight.y))
        {
           edgeLet.previousObjNo = p;
           pushed = true;
           break;
        }
      if (!pushed)
         edgeLet.previousObjNo = -1;
      vEdgelets.push_back(edgeLet);
//      std::cout << "Check : " << vEdgelets.back().v2Pose << " " << line.v2Start << " " << line.v2End << std::endl;
    }
    else
    {
//      REAL_TYPE rAngle = atan(line.rSlope);
      for( int j = 0; j < no_segment; ++ j )
      {
        Edgelet edgeLet;
        edgeLet.v2Pose = line.v2Start + (iMaxLength*j+halfMaxLength)*line.v2Direction;
        edgeLet.v2Pose[0] = round(edgeLet.v2Pose[0] ); edgeLet.v2Pose[1] = round(edgeLet.v2Pose[1]);

//        edgeLet.rAngle = rAngle;
        edgeLet.rSlope = line.rSlope;
        pushed = false;
        for (int p = 0; p < mvDetectedObjects_previousFrame.size(); p++)
          if (edgeLet.v2Pose[0] >= (mvDetectedObjects_previousFrame[p].irTopLeft.x) && edgeLet.v2Pose[0] <= (mvDetectedObjects_previousFrame[p].irBottomRight.x) && edgeLet.v2Pose[1] >= (mvDetectedObjects_previousFrame[p].irTopLeft.y) && edgeLet.v2Pose[1] <= (mvDetectedObjects_previousFrame[p].irBottomRight.y))
          {
             edgeLet.previousObjNo = p;
             pushed = true;
             break;
          }
        if (!pushed)
          edgeLet.previousObjNo = -1;
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
        if (edgeLet.v2Pose[0] > 0 && edgeLet.v2Pose[0] <= img_h && edgeLet.v2Pose[1] > 0 && edgeLet.v2Pose[1] <= img_w)
        {
//        edgeLet.rAngle = rAngle;
           edgeLet.rSlope = line.rSlope;
           edgeLet.regionNo = 0;
           pushed = false;
           for (int p = 0; p < mvDetectedObjects_previousFrame.size(); p++)
             if (edgeLet.v2Pose[0] >= (mvDetectedObjects_previousFrame[p].irTopLeft.x) && edgeLet.v2Pose[0] <= (mvDetectedObjects_previousFrame[p].irBottomRight.x) && edgeLet.v2Pose[1] >= (mvDetectedObjects_previousFrame[p].irTopLeft.y) && edgeLet.v2Pose[1] <= (mvDetectedObjects_previousFrame[p].irBottomRight.y))
             {
               edgeLet.previousObjNo = p;
               pushed = true;
               break;
             }
           if (!pushed)
             edgeLet.previousObjNo = -1;
           vEdgelets.push_back(edgeLet);
        }
      }
    }
  }
  viLineIndexes.push_back( vEdgelets.size() );
#if 0
  std::cout << "Edgelet ..." << std::endl;
  for( int i = 0; i < vEdgelets.size(); ++i )
    std::cout << vEdgelets[i].v2Pose << " " << vEdgelets[i].rSlope << " " << vEdgelets[i].rAngle << std::endl;
#endif
  double edgelet_time = CVD::timer.get_time();
  double distTransform1 = 0, distTransform2 = 0, distTransform3 = 0, distTransform4 = 0;
// Generate Edge Orientation Image 
  for( unsigned int i = 0; i < mvEdgeOrientationImages.size(); ++i )
  {
#if 0  
    char filename[100];
    sprintf(filename,"%s/EdgeOrientation%d.png",folderName,i);
    CVD::byte *dataDst = mEdgeImage.data();
//    CVD::byte *dataSrc = mvEdgeOrientationImages[i].data();
    REAL_TYPE *dataSrc = mvEdgeOrientations[i].data;    
    int no_data = mEdgeImage.size().x*mEdgeImage.size().y;
    while( no_data > 0 )
    {
      *dataDst = (*dataSrc)*255;
      ++dataSrc;
      ++dataDst;
      --no_data;
    }
    CVD::img_save(mEdgeImage, filename);
#endif
//    GetTransformDist::compute( mvEdgeOrientationDistTransformImages[i], mvEdgeOrientationImages[i] );
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
/*    
    CVD::BasicImage<REAL_TYPE>::iterator beginIter = mvEdgeOrientationDistTransformImages[i].begin();
    CVD::BasicImage<REAL_TYPE>::iterator endIter = mvEdgeOrientationDistTransformImages[i].end();
    while( beginIter != endIter )
    {
      if( max < *beginIter )
        max = *beginIter;
      ++beginIter;
    }
    beginIter = mvEdgeOrientationDistTransformImages[i].begin();
    while( beginIter != endIter )
    {
      *beginIter /= max; // normalize
      ++beginIter;
    }
*/  
  }
  double edge_ori_time = CVD::timer.get_time();
 // std::cout <<"Edgelet with Ori = Edgelet : " << edgelet_time - local_start_time << ", Ori : " << edge_ori_time - edgelet_time << std::endl;
 // std::cout<<"DistanceTransform takes: " << distTransform2 - distTransform1 << " then " << distTransform3 - distTransform2 << " then " << distTransform4 - distTransform3 << std::endl;
}

bool EdgeObjectDetection::calculateTwoEdgelets( TestDirection & testDirection, const std::vector<Edgelet> & vEdgelets, int iIdx_0, int iIdx_1 )
{
  const Edgelet & fstEdgelet = vEdgelets[iIdx_0];
  const Edgelet & sndEdgelet = vEdgelets[iIdx_1];
  const REAL_TYPE & m1 = fstEdgelet.rSlope; 
  const REAL_TYPE & m3 = sndEdgelet.rSlope;
  TooN::Vector<2>  v2Del = sndEdgelet.v2Pose - fstEdgelet.v2Pose;
  REAL_TYPE m2 = ( (v2Del[0]) != 0 ) ? (v2Del[1])/(v2Del[0]) : (v2Del[1])/0.01;

  testDirection.iIdx_0 = iIdx_0;
  testDirection.iIdx_1 = iIdx_1;
  testDirection.rDir_1 = atan( (m2-m1)/(1.f+m1*m2) ); // angle between fst edgelet and direction
/*  if( m1 == 10000000 )
  {
    testDirection.rDir_1 = atan(m2);
    if( testDirection.rDir_1 < 0 ) testDirection.rDir_1 = -M_PI*0.5 - testDirection.rDir_1;
    else testDirection.rDir_1 = M_PI*0.5 - testDirection.rDir_1;
  }
  else
    testDirection.rDir_1 = atan( (m2-m1)/1.f+m1*m2); */
  
  testDirection.rDist = norm(v2Del);
  testDirection.rAngle = atan( (m3 - m1)/(1.f+m1*m3) );

#if 1
  testDirection.rDir_2 = atan2(v2Del[1], v2Del[0]); 
#else  
  normalize(v2Del);
  REAL_TYPE test=0;
  if( v2Del[1] < -1 ) test = acos(-1);
  else if( v2Del[1] > 1 ) test = acos(1);
  else test = acos(v2Del[0]);
  if( asin(v2Del[1]) < 0 )
    test = -test;
  testDirection.rDir_2 = test;
#endif    
  return true;
}

bool EdgeObjectDetection::calculateEdgeletAndCell( CellDirection & cd, const std::vector<Edgelet> & vEdgelets, int edgeNo, int cellNo )
{
  const Edgelet & fstEdgelet = vEdgelets[edgeNo];
  const EdgeCellMatrix & sndEdgelet = ecm[cellNo];
  const REAL_TYPE & m1 = fstEdgelet.rSlope; 
  TooN::Vector<2>  v2Del = sndEdgelet.v2Pose - fstEdgelet.v2Pose;
  REAL_TYPE m2 = ( (v2Del[0]) != 0 ) ? (v2Del[1])/(v2Del[0]) : (v2Del[1])/0.01;

  cd.iIdx = edgeNo;
  cd.cellNo = cellNo;
  cd.rDir_1 = atan( (m2-m1)/(1.f+m1*m2) ); // angle between fst edgelet and direction
  cd.rDir_2 = atan2(v2Del[1], v2Del[0]); 
  return true;
}

int EdgeObjectDetection::calculateTwoCellsDirection (int cell1, int cell2, double cellInterval)
{
  TooN::Vector<2> v2Del = ecm[cell2].v2Pose - ecm[cell1].v2Pose;
  REAL_TYPE m2 = ( (v2Del[0]) != 0 ) ? (v2Del[1])/(v2Del[0]) : (v2Del[1])/0.01;
  int dir = (atan(m2)+pi/2)/cellInterval;
  //std::cout << "comparing (" << ecm[cell1].v2Pose[1] << ", " << ecm[cell1].v2Pose[0] << ") and (" << ecm[cell2].v2Pose[1] << ", " << ecm[cell1].v2Pose[0] << ") with " << m2 << " " << dir << std::endl;
  return dir;
}

void EdgeObjectDetection::gen_EdgeletCells( const std::vector<Edgelet> & vEdgelets )
{
  int no_edgelets = vEdgelets.size();
  mvCellDirections.clear();
  mvCellDirections.reserve( no_edgelets*cellCount );
  vCellTestDirections.clear();
  vCellTestDirections.reserve(no_edgelets*no_edgelets/2);

  CellDirection cd;
  for( unsigned int i = 0; i < no_edgelets; ++i )
  {
    for( unsigned int j = 0; j < cellCount; ++j )
    {
      if( calculateEdgeletAndCell( cd, vEdgelets, i, j ) )
        mvCellDirections.push_back(cd);
    }
  }
}

void EdgeObjectDetection::gen_TestDirections( const std::vector<Edgelet> & vEdgelets )
{
  std::cerr << "gen_TestDirections" << std::endl;
  int no_edgelets = vEdgelets.size();
  mvTestDirections.clear();
  mvTestDirections.reserve( no_edgelets*(no_edgelets-1) );
  mvirDirectionIndexes.clear(); 
  mvirDirectionIndexes.resize( vEdgelets.size() );
  if (vEdgelets.size() < 2)
    return;

  TestDirection td;
  for( unsigned int i = 0; i < no_edgelets; ++i )
  {
    mvirDirectionIndexes[i].x = mvTestDirections.size();
    for( unsigned int j = 0; j < no_edgelets; ++j )
    {
      if( i == j )
        continue;
      if( calculateTwoEdgelets( td, vEdgelets, i, j ) )
        mvTestDirections.push_back(td);
    }
  }
  if (mvirDirectionIndexes.size() == 0)
     return;
  int checkDirectionIndexes = 1;
  while( mvirDirectionIndexes[checkDirectionIndexes].x == 0 )
  {
    mvirDirectionIndexes[checkDirectionIndexes-1].x = -1;
    ++checkDirectionIndexes;
  }
  while( checkDirectionIndexes < (int)mvirDirectionIndexes.size() )
  {
    int idxValue = mvirDirectionIndexes[checkDirectionIndexes].x;
    mvirDirectionIndexes[checkDirectionIndexes-1].y = idxValue - 1;
    ++checkDirectionIndexes;
    while( idxValue == mvirDirectionIndexes[checkDirectionIndexes].x )
    {
      mvirDirectionIndexes[checkDirectionIndexes].x = -1;
      ++checkDirectionIndexes;
    }
  }
  mvirDirectionIndexes.back().y = mvTestDirections.size() - 1;
}

void EdgeObjectDetection::gen_TestDirections( const std::vector<Edgelet> & vEdgelets, const std::vector<int> & viLineIndexes )
{
  int no_edgelets = vEdgelets.size();
//  std::vector<int> vE1(no_edgelets), vE2(no_edgelets);
//  for( int i = 0; i < no_edgelets; ++ i )
//    vE1[i] = vE2[i] = i;
//  std::random_shuffle( vE1.begin(), vE1.end() );
//  std::random_shuffle( vE2.begin(), vE2.end() );
  mvTestDirections.clear();
  if (no_edgelets == 0)
      return;
  mvTestDirections.reserve( no_edgelets*(no_edgelets-1) );
  mvTestDirectionsPreviousObjects.clear();
  if (mvDetectedObjects_previousFrame.size() > 0)
     mvTestDirectionsPreviousObjects.reserve (mvDetectedObjects_previousFrame.size());
  for (int i = 0; i < mvDetectedObjects_previousFrame.size(); i++)
  {
     std::vector<TestDirection> thisObj;
     thisObj.reserve(no_edgelets*(no_edgelets-1));
     mvTestDirectionsPreviousObjects.push_back(thisObj);
  }
  mvirDirectionIndexes.clear(); 
  if (vEdgelets.size() == 0)
     return;
  mvirDirectionIndexes.resize( vEdgelets.size() );
  TestDirection td;
#if 1
  for( unsigned int i = 0; i < no_edgelets; ++i )
  {
    mvirDirectionIndexes[i].x = mvTestDirections.size();
    for( unsigned int j = 0; j < no_edgelets; ++j )
    {
      if (i == j || (vEdgelets[i].regionNo != vEdgelets[j].regionNo && (vEdgelets[i].previousObjNo != vEdgelets[j].previousObjNo || vEdgelets[i].previousObjNo < 0 || vEdgelets[j].previousObjNo < 0)))
        continue;
      if( calculateTwoEdgelets( td, vEdgelets, i, j ) )
      {
         mvTestDirections.push_back(td);
       /*  std::cerr << "B.3 " << vEdgelets[i].previousObjNo << std::endl;
         if (mvDetectedObjects_previousFrame.size() > 0)
         {
            if (vEdgelets[i].previousObjNo >= 0 && vEdgelets[j].previousObjNo >= 0 && vEdgelets[i].previousObjNo == vEdgelets[j].previousObjNo)
            {
               mvTestDirectionsPreviousObjects[vEdgelets[i].previousObjNo].push_back(td);
            }
         }*/
      }
    }
  }
#else

  for( unsigned int i = 1; i < viLineIndexes.size(); ++i )
  {
    int begin = viLineIndexes[i-1];
    int end = viLineIndexes[i];
    for( int fst = begin; fst < end; ++fst )
    {
      mvirDirectionIndexes[fst].x = mvTestDirections.size();
      int snd = 0;
      while( snd < begin )
      {
        if( calculateTwoEdgelets( td, vEdgelets, fst, snd ) )
          mvTestDirections.push_back(td);
        ++snd;
      }
      snd = end;
      while( snd < no_edgelets )
      {
        if( calculateTwoEdgelets( td, vEdgelets, fst, snd ) )
          mvTestDirections.push_back(td);
        ++snd;
      }
    }
  }
#endif
  if (mvirDirectionIndexes.size() == 0)
     return;
  int checkDirectionIndexes = 1;
  while( mvirDirectionIndexes[checkDirectionIndexes].x == 0 )
  {
    mvirDirectionIndexes[checkDirectionIndexes-1].x = -1;
    ++checkDirectionIndexes;
  }
  while( checkDirectionIndexes < (int)mvirDirectionIndexes.size() )
  {
    int idxValue = mvirDirectionIndexes[checkDirectionIndexes].x;
    mvirDirectionIndexes[checkDirectionIndexes-1].y = idxValue - 1;
    ++checkDirectionIndexes;
    while( idxValue == mvirDirectionIndexes[checkDirectionIndexes].x )
    {
      mvirDirectionIndexes[checkDirectionIndexes].x = -1;
      ++checkDirectionIndexes;
    }
  }
  mvirDirectionIndexes.back().y = mvTestDirections.size() - 1;
//  std::cout <<"No pairs " << no_edgelets*(no_edgelets-1) << " " << mvTestDirections.size() << std::endl;
 // for (int i = 0; i < mvTestDirectionsPreviousObjects.size(); i++)
 //    std::cout << "FOR PREVIOUS OBJECT " << i << " there are " << mvTestDirectionsPreviousObjects[i].size() << " pairs" << std::endl;
}

void EdgeObjectDetection::gen_CertainDirectionsChain_star( std::vector<CodeBookElement> & mvCodeBook, const std::vector<Edgelet> & vEdgelets, const TooN::Vector<4> & v4ChainAngles, int view_ID, int object_ID, int obj_correct_ID )
{
  int beginFstIndex = 0;
  int endFstIndex = mvTestDirections.size() - 1;
  int count=0;
  int pass_first = 0;
  for( ; beginFstIndex <= endFstIndex; ++ beginFstIndex )
  {
    const TestDirection & firstTestDirectionPairs = mvTestDirections[beginFstIndex];
    if( (fabs( firstTestDirectionPairs.rDir_1-v4ChainAngles[0] ) < EPS_DIR || fabs( firstTestDirectionPairs.rDir_1+v4ChainAngles[0] ) < EPS_DIR ) && firstTestDirectionPairs.rDist < FIRST_DIST_TRAIN )
    {
      ++pass_first;
      int beginSndIndex = mvirDirectionIndexes[firstTestDirectionPairs.iIdx_0].x;
      if( beginSndIndex == -1 ) continue;
      int endSndIndex = mvirDirectionIndexes[firstTestDirectionPairs.iIdx_0].y;

      for( ; beginSndIndex <= endSndIndex; ++beginSndIndex )
      {
        const TestDirection & sndTestDirectionPairs = mvTestDirections[beginSndIndex];
        if( fabs( calculateAngle(sndTestDirectionPairs.rDir_2,firstTestDirectionPairs.rDir_2) - v4ChainAngles[1] ) < EPS_DIR )
        {

          int beginTrdIndex = mvirDirectionIndexes[sndTestDirectionPairs.iIdx_0].x;
          if( beginTrdIndex == -1 ) continue;
          int endTrdIndex = mvirDirectionIndexes[sndTestDirectionPairs.iIdx_0].y;
          for( ; beginTrdIndex <= endTrdIndex; ++beginTrdIndex )
          {
            const TestDirection & trdTestDirectionPairs = mvTestDirections[beginTrdIndex];
            if( fabs( calculateAngle(trdTestDirectionPairs.rDir_2,sndTestDirectionPairs.rDir_2) - v4ChainAngles[2] ) < EPS_DIR ) 
            {
              
              if( isColinear( vEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose, vEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose, vEdgelets[trdTestDirectionPairs.iIdx_1].v2Pose ) ) continue;
              int beginFrtIndex = mvirDirectionIndexes[trdTestDirectionPairs.iIdx_0].x;
              if( beginFrtIndex == -1 ) continue;
              int endFrtIndex = mvirDirectionIndexes[trdTestDirectionPairs.iIdx_0].y;
              for( ; beginFrtIndex <= endFrtIndex; ++beginFrtIndex )
              {
                const TestDirection & fourtTestDirectionPairs = mvTestDirections[beginFrtIndex];
                if( fabs( calculateAngle(fourtTestDirectionPairs.rDir_2,trdTestDirectionPairs.rDir_2) - v4ChainAngles[3] ) < EPS_DIR )
                {
                  if( isColinear( vEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose, vEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose, vEdgelets[fourtTestDirectionPairs.iIdx_1].v2Pose ) || isColinear( vEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose, vEdgelets[trdTestDirectionPairs.iIdx_1].v2Pose, vEdgelets[fourtTestDirectionPairs.iIdx_1].v2Pose ) || isColinear( vEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose, vEdgelets[trdTestDirectionPairs.iIdx_1].v2Pose, vEdgelets[fourtTestDirectionPairs.iIdx_1].v2Pose ) )
                    continue;
                  ++count;
                  CodeBookElement cd;
                  cd.iIdx_0 = firstTestDirectionPairs.iIdx_0;
                  cd.iIdx_1 = firstTestDirectionPairs.iIdx_1;
                  cd.iIdx_2 = sndTestDirectionPairs.iIdx_1;
                  cd.iIdx_3 = trdTestDirectionPairs.iIdx_1;
                  cd.iIdx_4 = fourtTestDirectionPairs.iIdx_1;
                  cd.rAngle_0 = firstTestDirectionPairs.rAngle;
                
                  cd.rAngle_1 = sndTestDirectionPairs.rAngle;
                  cd.rRelativeDist_1 = sndTestDirectionPairs.rDist/firstTestDirectionPairs.rDist;

                  cd.rAngle_2 = trdTestDirectionPairs.rAngle;
                  cd.rRelativeDist_2 = trdTestDirectionPairs.rDist/sndTestDirectionPairs.rDist;

                  cd.rAngle_3 = fourtTestDirectionPairs.rAngle;
                  cd.rRelativeDist_3 = fourtTestDirectionPairs.rDist/trdTestDirectionPairs.rDist;

                  cd.iView_ID = view_ID;
                  cd.iObject_ID = object_ID;
                  cd.obj_correct_ID = obj_correct_ID;

                  mvCodeBook.push_back( cd );
                }
              }
            }
          }
        }
      }
    }
  }
  std::cerr <<"No code book for this view : " << count << " Pass : " << pass_first << std::endl;
}

void EdgeObjectDetection::gen_CertainDirectionsChain( std::vector<CodeBookElement> & mvCodeBook, const std::vector<Edgelet> & vEdgelets, const TooN::Vector<4> & v4ChainAngles, int view_ID, int object_ID, int obj_correct_ID )
{
  int beginFstIndex = 0;
  int endFstIndex = mvTestDirections.size() - 1;
  int count=0;
  int pass_first = 0;
  for( ; beginFstIndex <= endFstIndex; ++ beginFstIndex )
  {
    const TestDirection & firstTestDirectionPairs = mvTestDirections[beginFstIndex];
    if( (fabs( firstTestDirectionPairs.rDir_1-v4ChainAngles[0] ) < EPS_DIR || fabs( firstTestDirectionPairs.rDir_1+v4ChainAngles[0] ) < EPS_DIR ) && firstTestDirectionPairs.rDist < FIRST_DIST_TRAIN )
//    if( (fabs( firstTestDirectionPairs.rDir_1 - v4ChainAngles[0] ) < EPS_DIR) && firstTestDirectionPairs.rDist < FIRST_DIST_TRAIN )
    {
//      std::cout << firstTestDirectionPairs.rAngle << " " << firstTestDirectionPairs.rDir_1 << " " << firstTestDirectionPairs.iIdx_0 << " " << firstTestDirectionPairs.iIdx_1 << std::endl;
//      std::cout << firstTestDirectionPairs.iIdx_0 << " " << firstTestDirectionPairs.iIdx_1 << " " << vEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose << " " << vEdgelets[firstTestDirectionPairs.iIdx_1].v2Pose << " " << firstTestDirectionPairs.rDir_1 << std::endl;
      ++pass_first;
      int beginSndIndex = mvirDirectionIndexes[firstTestDirectionPairs.iIdx_1].x;
      if( beginSndIndex == -1 ) continue;
      int endSndIndex = mvirDirectionIndexes[firstTestDirectionPairs.iIdx_1].y;
//      std::cout <<"Test : " << endSndIndex << " " << beginSndIndex << std::endl;
      for( ; beginSndIndex <= endSndIndex; ++beginSndIndex )
      {
        const TestDirection & sndTestDirectionPairs = mvTestDirections[beginSndIndex];
//        std::cout << fabs( calculateTwoEdgelets(sndTestDirectionPairs.rDir_2,firstTestDirectionPairs.rDir_2) - v4ChainAngles[1] ) << std::endl;
        if( fabs( calculateAngle(sndTestDirectionPairs.rDir_2,firstTestDirectionPairs.rDir_2) - v4ChainAngles[1] ) < EPS_DIR )
        {
//          std::cout << vEdgelets[sndTestDirectionPairs.iIdx_0].v2Pose << " " << vEdgelets[sndTestDirectionPairs.iIdx_1].v2Pose << " " << sndTestDirectionPairs.rDir_2 << " " << firstTestDirectionPairs.rDir_2 << " " << calculateAngle(sndTestDirectionPairs.rDir_2, firstTestDirectionPairs.rDir_2) << std::endl;
//          std::cout << firstTestDirectionPairs.rAngle << " " << firstTestDirectionPairs.rDir_1 << " " << firstTestDirectionPairs.iIdx_0 << " " << firstTestDirectionPairs.iIdx_1 << " " << v4ChainAngles[1] << std::endl;
//          ++pass_first;

          int beginTrdIndex = mvirDirectionIndexes[sndTestDirectionPairs.iIdx_1].x;
          if( beginTrdIndex == -1 ) continue;
          int endTrdIndex = mvirDirectionIndexes[sndTestDirectionPairs.iIdx_1].y;
          for( ; beginTrdIndex <= endTrdIndex; ++beginTrdIndex )
          {
            const TestDirection & trdTestDirectionPairs = mvTestDirections[beginTrdIndex];
/*            if( firstTestDirectionPairs.iIdx_0 == 45 && firstTestDirectionPairs.iIdx_1 == 47 && sndTestDirectionPairs.iIdx_1 == 32 )
            {
              std::cout << "Check " << trdTestDirectionPairs.iIdx_1 << " " << trdTestDirectionPairs.rDir_2 << " " << sndTestDirectionPairs.rDir_2 << " " << calculateAngle(trdTestDirectionPairs.rDir_2, sndTestDirectionPairs.rDir_2) << std::endl;
            } */
            if( fabs( calculateAngle(trdTestDirectionPairs.rDir_2,sndTestDirectionPairs.rDir_2) - v4ChainAngles[2] ) < EPS_DIR ) 
            {
              if( isColinear( vEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose, vEdgelets[sndTestDirectionPairs.iIdx_0].v2Pose, vEdgelets[trdTestDirectionPairs.iIdx_0].v2Pose ) ) continue;
              int beginFrtIndex = mvirDirectionIndexes[trdTestDirectionPairs.iIdx_1].x;
              if( beginFrtIndex == -1 ) continue;
              int endFrtIndex = mvirDirectionIndexes[trdTestDirectionPairs.iIdx_1].y;
              for( ; beginFrtIndex <= endFrtIndex; ++beginFrtIndex )
              {
                const TestDirection & fourtTestDirectionPairs = mvTestDirections[beginFrtIndex];
                if( fabs( calculateAngle(fourtTestDirectionPairs.rDir_2,trdTestDirectionPairs.rDir_2) - v4ChainAngles[3] ) < EPS_DIR )
                {
                  if( isColinear( vEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose, vEdgelets[sndTestDirectionPairs.iIdx_0].v2Pose, vEdgelets[fourtTestDirectionPairs.iIdx_0].v2Pose ) || isColinear( vEdgelets[firstTestDirectionPairs.iIdx_0].v2Pose, vEdgelets[trdTestDirectionPairs.iIdx_0].v2Pose, vEdgelets[fourtTestDirectionPairs.iIdx_0].v2Pose ) || isColinear( vEdgelets[sndTestDirectionPairs.iIdx_0].v2Pose, vEdgelets[trdTestDirectionPairs.iIdx_0].v2Pose, vEdgelets[fourtTestDirectionPairs.iIdx_0].v2Pose ) )
                    continue;
//                  std::cout << firstTestDirectionPairs.rDir_1 << " " << vEdgelets[firstTestDirectionPairs.iIdx_0].rSlope << std::endl;
                  ++count;
                  CodeBookElement cd;
                  cd.iIdx_0 = firstTestDirectionPairs.iIdx_0;
                  cd.iIdx_1 = firstTestDirectionPairs.iIdx_1;
                  cd.iIdx_2 = sndTestDirectionPairs.iIdx_1;
                  cd.iIdx_3 = trdTestDirectionPairs.iIdx_1;
                  cd.iIdx_4 = fourtTestDirectionPairs.iIdx_1;
                  cd.rAngle_0 = firstTestDirectionPairs.rAngle;
                
                  cd.rAngle_1 = sndTestDirectionPairs.rAngle;
                  cd.rRelativeDist_1 = sndTestDirectionPairs.rDist/firstTestDirectionPairs.rDist;

                  cd.rAngle_2 = trdTestDirectionPairs.rAngle;
                  cd.rRelativeDist_2 = trdTestDirectionPairs.rDist/sndTestDirectionPairs.rDist;

                  cd.rAngle_3 = fourtTestDirectionPairs.rAngle;
                  cd.rRelativeDist_3 = fourtTestDirectionPairs.rDist/trdTestDirectionPairs.rDist;

                  cd.iView_ID = view_ID;
                  cd.iObject_ID = object_ID;
                  cd.obj_correct_ID = obj_correct_ID;

                  mvCodeBook.push_back( cd );
                }
              }
            }
          }
        }
      }
    }
  }
  std::cerr <<"No code book for this view : " << count << " Pass : " << pass_first << std::endl;
}

bool sort_codebook( const CodeBookElement & left, const CodeBookElement & right )
{
  return left.rAngle_0 < right.rAngle_0;
}

bool lower_codebook( const CodeBookElement & left, REAL_TYPE value )
{
  return left.rAngle_0 < value;
}

bool lower_codebook2 (const CodeBookElement & left, REAL_TYPE value)
{
  return left.rAngle_1 < value;
}

bool lower_codebook3 (const CodeBookElement & left, REAL_TYPE value )
{
  return left.rRelativeDist_1 < value;
}

void EdgeObjectDetection::sortCodeBooks()
{
  for( unsigned int no_class = 0; no_class < mvObjDescriptorClasses.size(); ++ no_class )
  {
    std::vector<CodeBookElement> & mvCodeBook = mvObjDescriptorClasses[no_class].mvCodeBook;

    std::cout << "Sort no. " << mvCodeBook.size() << " elements. " << std::endl; 
    std::sort( mvCodeBook.begin(), mvCodeBook.end(), sort_codebook );
    std::vector<CVD::ImageRef> & mvirFirstAngleIndexes = mvObjDescriptorClasses[no_class].mvirFirstAngleIndexes;
    mvirFirstAngleIndexes.resize(64);
    std::vector< std::vector<int> > & mvirSecondIndexes = mvObjDescriptorClasses[no_class].mvirSecondIndexes;
    mvirSecondIndexes.resize(163840);

    REAL_TYPE rStep = M_PI/64.f;
    REAL_TYPE rTheta = -0.5*M_PI;
    mvirFirstAngleIndexes[0].x = 0;
    std::vector<CodeBookElement>::iterator beginIter = mvCodeBook.begin();
    std::vector<CodeBookElement>::iterator endIter = mvCodeBook.end();
    REAL_TYPE rTheta2, rDist; 
    REAL_TYPE dStep1 = 0.1;
    REAL_TYPE dStep2 = 1;

    int counter = 0;
    int ind1 = 0, ind2 = 0, ind3 = 0;
    for( int i = 1; i < 64; ++i )
    {
      rTheta += rStep;
      std::vector<CodeBookElement>::iterator lowerIter = std::lower_bound( beginIter, endIter, rTheta, lower_codebook );
      mvirFirstAngleIndexes[i].x = int(lowerIter-mvCodeBook.begin());
      mvirFirstAngleIndexes[i-1].y = mvirFirstAngleIndexes[i].x - 1;

      rTheta2 = -0.5*M_PI;
      ind2 = 0;
      while (rTheta2 <= 0.5*M_PI)
      {
         rTheta2 += rStep;
         REAL_TYPE rDist = 0;
         ind3 = 0;
         while (rDist < 30)
         {
              REAL_TYPE dStep;
              
              if (rDist < 1)
                 dStep = dStep1;
              else
                 dStep = dStep2;
              rDist += dStep;
              std::vector<CodeBookElement>::iterator loop = mvirFirstAngleIndexes[i-1].x+mvCodeBook.begin();
              std::vector<int> & mvirVals = mvirSecondIndexes[counter];
              mvirVals.reserve (int(lowerIter-beginIter));
              while (loop < (mvirFirstAngleIndexes[i].x+mvCodeBook.begin()))
              {
                  if( (*loop).rAngle_1 > (rTheta2-rStep) && (*loop).rAngle_1 < (rTheta2) && (*loop).rRelativeDist_1 > (rDist-dStep) && (*loop).rRelativeDist_1 < rDist)
                  {
                     mvirVals.push_back (int(loop-mvCodeBook.begin()));
                  }
                  loop++;
              }   
              counter++;
              ind3++;
         }
         ind2++;
      }
      ind1++;
      beginIter = lowerIter;
    }
    mvirFirstAngleIndexes[63].y = mvCodeBook.size() - 1;

  } 
}

int EdgeObjectDetection::resetLatestObj ()
{
   latestObjID = mvObjDescriptorClasses[0].mObjectsTemplate.iNoObjects + 1;
   latest_correctObjID = latestObjID;
   return (latestObjID-1);
}

void EdgeObjectDetection::setLatestObjectID(int i)
{
   latest_correctObjID = i;
   if (i < 3)
      latestObjID = i - 1;
   else
      latestObjID = i - 3;
}

bool EdgeObjectDetection::addView (int objNo, SegmentedRegion region, bool andSave)
{
 // std::cout << "entering add view with object " << objNo << std::endl;
  for( unsigned int no_class = 0; no_class < mvObjDescriptorClasses.size(); ++no_class )
  {
     TooN::Vector<4> & mv4ChainAngles = mvObjDescriptorClasses[no_class].mv4ChainAngles;
     int sizeBefore = mvObjDescriptorClasses[no_class].mvCodeBook.size();
     std::vector<CodeBookElement> & mvCodeBook = mvObjDescriptorClasses[no_class].mvCodeBook;
     ObjectsTemplate & mObjectsTemplate = mvObjDescriptorClasses[no_class].mObjectsTemplate;
     if (mObjectsTemplate.vObjectElements.size() <= objNo)
     {
        mObjectsTemplate.vObjectElements.resize (mObjectsTemplate.iNoObjects+1);
        ObjectElement & newObject = mObjectsTemplate.vObjectElements[mObjectsTemplate.iNoObjects];
        newObject.iObject_ID = mObjectsTemplate.iNoObjects;
        newObject.obj_correct_ID = objNo;
        newObject.iNoViews = 0;
        latestObjID = mObjectsTemplate.iNoObjects;
        mObjectsTemplate.iNoObjects++;
     }
     learntObjNo = objNo + 1;
     ObjectElement & obj = mObjectsTemplate.vObjectElements[objNo];
     int view_ID = obj.iNoViews;
    // std::cerr << "adding view: " << objNo << " " << view_ID << " " << obj.iNoViews << std::endl;
     obj.iNoViews++;
     obj.vViews.resize(obj.iNoViews);
     View & view = obj.vViews[view_ID];
     view.iView_ID = view_ID;
     view.vEdgelets.reserve (mvEdgelets.size());
     for (int i = 0; i < mvEdgelets.size(); i++)
     {
        if (mvEdgelets[i].v2Pose[0] >= (region.irTopLeft.x/imgSizeRatio) && mvEdgelets[i].v2Pose[0] <= (region.irBottomRight.x/imgSizeRatio) && mvEdgelets[i].v2Pose[1] >= (region.irTopLeft.y/imgSizeRatio) && mvEdgelets[i].v2Pose[1] <= (region.irBottomRight.y/imgSizeRatio))
        {
           view.vEdgelets.push_back (mvEdgelets[i]);
        }
     }
     int counterIndexes = 1;
     view.allEdges.viEdgeIdxes.reserve (edges.virEdges.size());
     view.allEdges.virEdges.reserve (edges.virEdges.size());
        for (int i = 0; i < edges.virEdges.size(); i++)
        {
           if (edges.virEdges[i].x >= region.irTopLeft.x/imgSizeRatio && edges.virEdges[i].x <= region.irBottomRight.x/imgSizeRatio && edges.virEdges[i].y >= region.irTopLeft.y/imgSizeRatio && edges.virEdges[i].y <= region.irBottomRight.y/imgSizeRatio)
            {
               view.allEdges.viEdgeIdxes.push_back (counterIndexes);
               counterIndexes++;
               view.allEdges.virEdges.push_back (edges.virEdges[i]);
          //     view.allEdges.virEdges.back().x /= imgSizeRatio;
          //     view.allEdges.virEdges.back().y /= imgSizeRatio;
            }
        }
     view.iObject_ID = objNo;
     if (isStar)
        gen_CertainDirectionsChain_star(mvCodeBook, view.vEdgelets, mv4ChainAngles, view_ID, objNo, objNo);
     else
        gen_CertainDirectionsChain(mvCodeBook, view.vEdgelets, mv4ChainAngles, view_ID, objNo, objNo);
     
     for (int a = sizeBefore; a < mvCodeBook.size(); a++)
     {
          CodeBookElement h = mvCodeBook[a];
          int indexNo = 64.f*(h.rAngle_0/M_PI + 0.5f); if( indexNo < 0 ) indexNo = 0;
          int indexNo2 = 64.f*(h.rAngle_1/M_PI + 0.5f); if( indexNo2 < 0 ) indexNo2 = 0;
          int indexNo3;
          REAL_TYPE ratio = h.rRelativeDist_1;
          if (ratio < 1)
              indexNo3 = int(ratio/0.1);
          else if (ratio > 1 && ratio < 29.5)
              indexNo3 = int (ratio)+9;
          else
              indexNo3 = 38;
          std::vector< std::vector<int> > & mvirSecondIndexes = mvObjDescriptorClasses[no_class].mvirSecondIndexes;
          std::vector<int> & thisList = mvirSecondIndexes[(indexNo)*2496+(indexNo2)*39+indexNo3];
          thisList.reserve(thisList.size()+1);
          thisList.push_back(a);
     }     
  }
  viewAdded = true;
  learntViewNo++;
  return true;
}

bool EdgeObjectDetection::addView (int objNo, SegmentedRegion region, SegmentedMask mask, bool AndSave)
{
  std::cerr << "entering add view with object " << objNo << std::endl;
  int MAX_EDGELET_SIZE = 300; // to be changed later!! - BEWARE
  for( unsigned int no_class = 0; no_class < mvObjDescriptorClasses.size(); ++no_class )
  {
     TooN::Vector<4> & mv4ChainAngles = mvObjDescriptorClasses[no_class].mv4ChainAngles;
     int sizeBefore = mvObjDescriptorClasses[no_class].mvCodeBook.size();
     std::vector<CodeBookElement> & mvCodeBook = mvObjDescriptorClasses[no_class].mvCodeBook;
     ObjectsTemplate & mObjectsTemplate = mvObjDescriptorClasses[no_class].mObjectsTemplate;
     if (mObjectsTemplate.vObjectElements.size() <= objNo)
     {
        mObjectsTemplate.vObjectElements.resize (mObjectsTemplate.iNoObjects+1);
        ObjectElement & newObject = mObjectsTemplate.vObjectElements[mObjectsTemplate.iNoObjects];
        newObject.iObject_ID = mObjectsTemplate.iNoObjects;
        newObject.obj_correct_ID = objNo;
        newObject.iNoViews = 0;
        latestObjID = mObjectsTemplate.iNoObjects;
        mObjectsTemplate.iNoObjects++;
     }
     learntObjNo = objNo+1;
     ObjectElement & obj = mObjectsTemplate.vObjectElements[objNo];
     int view_ID = obj.iNoViews;
     std::cerr << "adding view: " << objNo << " " << view_ID << " " << obj.iNoViews << std::endl;
     obj.iNoViews++;
     obj.vViews.resize(obj.iNoViews);
     View & view = obj.vViews[view_ID];
     view.iView_ID = view_ID;
     view.vEdgelets.reserve (mvEdgelets.size());
     if (isMask)
     {
       for (int i = 0; i < mvEdgelets.size(); i++)
       {
          if (mask.image.data[imgSizeRatio*mvEdgelets[i].v2Pose[1]*mask.image.width+imgSizeRatio*mvEdgelets[i].v2Pose[0]] > 0)
              if (view.vEdgelets.size() < MAX_EDGELET_SIZE)
              {
                 view.vEdgelets.push_back (mvEdgelets[i]);
              }
              else
                break;
       }
     }
     else
     {
       for (int i = 0; i < mvEdgelets.size(); i++)
       {
           if (mvEdgelets[i].v2Pose[0] >= (region.irTopLeft.x/imgSizeRatio) && mvEdgelets[i].v2Pose[0] <= (region.irBottomRight.x/imgSizeRatio) && mvEdgelets[i].v2Pose[1] >= (region.irTopLeft.y/imgSizeRatio) && mvEdgelets[i].v2Pose[1] <= (region.irBottomRight.y/imgSizeRatio))
           {
              view.vEdgelets.push_back (mvEdgelets[i]);
           }
       }
     }
     std::cerr << "Number of edgelets within mask " << view.vEdgelets.size() << std::endl;
     int counterIndexes = 1;
     view.allEdges.viEdgeIdxes.reserve (edges.virEdges.size());
     view.allEdges.virEdges.reserve (edges.virEdges.size());
     if (isMask)
     {
        for (int i = 0; i < edges.virEdges.size(); i++)
        {
            if (mask.image.data[2*edges.virEdges[i].y*mask.image.width+2*edges.virEdges[i].x] > 0)
            {
               view.allEdges.viEdgeIdxes.push_back (counterIndexes);
               counterIndexes++;
               view.allEdges.virEdges.push_back (edges.virEdges[i]);
             //  view.allEdges.virEdges.back().x /= imgSizeRatio;
             //  view.allEdges.virEdges.back().y /= imgSizeRatio;
            }
        }

     }
     else
     {
        for (int i = 0; i < edges.virEdges.size(); i++)
        {
            if (edges.virEdges[i].x >= region.irTopLeft.x/imgSizeRatio && edges.virEdges[i].x <= region.irBottomRight.x/imgSizeRatio && edges.virEdges[i].y >= region.irTopLeft.y/imgSizeRatio && edges.virEdges[i].y <= region.irBottomRight.y/imgSizeRatio)
            {
               view.allEdges.viEdgeIdxes.push_back (counterIndexes);
               counterIndexes++;
               view.allEdges.virEdges.push_back (edges.virEdges[i]);
            //   view.allEdges.virEdges.back().x /= imgSizeRatio;
            //   view.allEdges.virEdges.back().y /= imgSizeRatio;
            }
        }
     }
     view.iObject_ID = objNo;
     if (!isStar)
        gen_CertainDirectionsChain(mvCodeBook, view.vEdgelets, mv4ChainAngles, view_ID, objNo, objNo);
     else
        gen_CertainDirectionsChain_star(mvCodeBook, view.vEdgelets, mv4ChainAngles, view_ID, objNo, objNo);
     for (int a = sizeBefore; a < mvCodeBook.size(); a++)
     {
          CodeBookElement h = mvCodeBook[a];
          int indexNo = 64.f*(h.rAngle_0/M_PI + 0.5f); if( indexNo < 0 ) indexNo = 0;
          int indexNo2 = 64.f*(h.rAngle_1/M_PI + 0.5f); if( indexNo2 < 0 ) indexNo2 = 0;
          int indexNo3;
          REAL_TYPE ratio = h.rRelativeDist_1;
          if (ratio < 1)
              indexNo3 = int(ratio/0.1);
          else if (ratio > 1 && ratio < 29.5)
              indexNo3 = int (ratio)+9;
          else
              indexNo3 = 38;
          std::vector< std::vector<int> > & mvirSecondIndexes = mvObjDescriptorClasses[no_class].mvirSecondIndexes;
          std::vector<int> & thisList = mvirSecondIndexes[(indexNo)*2496+(indexNo2)*39+indexNo3];
          thisList.reserve(thisList.size()+1);
          thisList.push_back(a);
     }     
  }
  viewAdded = true;
  learntViewNo++;
  return true;
}

bool EdgeObjectDetection::addView ()
{
  for( unsigned int no_class = 0; no_class < mvObjDescriptorClasses.size(); ++no_class )
  {
     TooN::Vector<4> & mv4ChainAngles = mvObjDescriptorClasses[no_class].mv4ChainAngles;
     int sizeBefore = mvObjDescriptorClasses[no_class].mvCodeBook.size();
     std::vector<CodeBookElement> & mvCodeBook = mvObjDescriptorClasses[no_class].mvCodeBook;
     ObjectsTemplate & mObjectsTemplate = mvObjDescriptorClasses[no_class].mObjectsTemplate;
     if (mObjectsTemplate.vObjectElements.size() <= latestObjID)
     {
        mObjectsTemplate.vObjectElements.resize (mObjectsTemplate.iNoObjects+1);
        ObjectElement & newObject = mObjectsTemplate.vObjectElements[mObjectsTemplate.iNoObjects];
        newObject.iObject_ID = mObjectsTemplate.iNoObjects;
        newObject.obj_correct_ID = 11;
        newObject.iNoViews = 0;
        latestObjID = mObjectsTemplate.iNoObjects;
        mObjectsTemplate.iNoObjects++;
     }
     learntObjNo = latestObjID+1;
     ObjectElement & obj = mObjectsTemplate.vObjectElements[latestObjID];
     int view_ID = obj.iNoViews;
     std::cout << "view id is " << latestObjID << " " << view_ID << " " << obj.iNoViews << std::endl;
     obj.iNoViews++;
     obj.vViews.resize(obj.iNoViews);
     View & view = obj.vViews[view_ID];
     view.iView_ID = view_ID;
     view.vEdgelets.reserve (mvEdgelets.size());
     for (int i = 0; i < mvEdgelets.size(); i++)
     {
         if (mvEdgelets[i].v2Pose[0] >= (thisIrTopLeft.x/imgSizeRatio) && mvEdgelets[i].v2Pose[0] <= (thisIrBottomRight.x/imgSizeRatio) && mvEdgelets[i].v2Pose[1] >= (thisIrTopLeft.y/imgSizeRatio) && mvEdgelets[i].v2Pose[1] <= (thisIrBottomRight.y/imgSizeRatio))
         {
            view.vEdgelets.push_back (mvEdgelets[i]);

         }
     }
     int counterIndexes = 1;
     view.allEdges.viEdgeIdxes.reserve (edges.virEdges.size());
     view.allEdges.virEdges.reserve (edges.virEdges.size());
     for (int i = 0; i < edges.virEdges.size(); i++)
     {
         if (edges.virEdges[i].x >= thisIrTopLeft.x && edges.virEdges[i].x <= thisIrBottomRight.x && edges.virEdges[i].y >= thisIrTopLeft.y && edges.virEdges[i].y <= thisIrBottomRight.y)
         {
            view.allEdges.viEdgeIdxes.push_back (counterIndexes);
            counterIndexes++;
            view.allEdges.virEdges.push_back (edges.virEdges[i]);
         //   view.allEdges.virEdges.back().x /= imgSizeRatio;
         //   view.allEdges.virEdges.back().y /= imgSizeRatio;
         }
     }
     view.iObject_ID = latestObjID;
     gen_CertainDirectionsChain(mvCodeBook, view.vEdgelets, mv4ChainAngles, view_ID, latestObjID, latest_correctObjID);
     for (int a = sizeBefore; a < mvCodeBook.size(); a++)
     {
          CodeBookElement h = mvCodeBook[a];
          int indexNo = 64.f*(h.rAngle_0/M_PI + 0.5f); if( indexNo < 0 ) indexNo = 0;
          int indexNo2 = 64.f*(h.rAngle_1/M_PI + 0.5f); if( indexNo2 < 0 ) indexNo2 = 0;
          int indexNo3;
          REAL_TYPE ratio = h.rRelativeDist_1;
          if (ratio < 1)
              indexNo3 = int(ratio/0.1);
          else if (ratio > 1 && ratio < 29.5)
              indexNo3 = int (ratio)+9;
          else
              indexNo3 = 38;
          std::vector< std::vector<int> > & mvirSecondIndexes = mvObjDescriptorClasses[no_class].mvirSecondIndexes;
          std::vector<int> & thisList = mvirSecondIndexes[(indexNo)*2496+(indexNo2)*39+indexNo3];
          thisList.reserve(thisList.size()+1);
          thisList.push_back(a);
     }     
  }
  viewAdded = true;
  learntViewNo++;
  return true;
}

int EdgeObjectDetection::addCodeBook( LowLevelImageData<CVD::byte> & llImage, CVD::ImageRef irTopLeft, CVD::ImageRef irButtomRight, const TooN::Vector<4> & v4ChainAngles, int object_ID, int obj_correct_ID )
{
  Edges edges, edgesInsideBox;
  std::vector<Line> vLines;
if (!isLSD) //#ifndef LSD 
{
    mCannyEdgeDetector.compute( mEdgeImage, edges, llImage, 8, 0.5, 0.3 ); //10, 0.8, 0.6
  int max_LinkedEdges = edges.viEdgeIdxes.size(), startIdx, endIdx;
//  for( int i = 0; i < edges.virEdges.size(); ++i )
//    std::cout << edges.virEdges[i].x << " " << edges.virEdges[i].y << std::endl;
  edgesInsideBox.viEdgeIdxes.push_back( 0 );
  for( int no_LinkedEdges = 1; no_LinkedEdges < max_LinkedEdges; ++ no_LinkedEdges )
  {
    startIdx = edges.viEdgeIdxes[ no_LinkedEdges - 1 ];
    endIdx = edges.viEdgeIdxes[ no_LinkedEdges ];
    for( int idx = startIdx; idx < endIdx; ++idx )
    {
      if( edges.virEdges[idx].x > irTopLeft.x && edges.virEdges[idx].y > irTopLeft.y && edges.virEdges[idx].x < irButtomRight.x && edges.virEdges[idx].y < irButtomRight.y )
      {
        edgesInsideBox.virEdges.push_back( edges.virEdges[idx] );
      }
      else
      {
        if( edgesInsideBox.viEdgeIdxes.back() != edgesInsideBox.virEdges.size() )
          edgesInsideBox.viEdgeIdxes.push_back( edgesInsideBox.virEdges.size() );
        do
        {
          ++idx;
        }while( idx < endIdx && (edges.virEdges[idx].x < irTopLeft.x || edges.virEdges[idx].y < irTopLeft.y || edges.virEdges[idx].x > irButtomRight.x || edges.virEdges[idx].y > irButtomRight.y) );
      }
    }
    if( edgesInsideBox.viEdgeIdxes.back() != edgesInsideBox.virEdges.size() )
      edgesInsideBox.viEdgeIdxes.push_back( edgesInsideBox.virEdges.size() );
  }
  get_lines( vLines, edgesInsideBox );
}

else
{
  int no_data = 0;
  double * img_ptr = lsd_image->data;
  CVD::byte * int_img_ptr = llImage.mImage.data();
  while( no_data != data_size )
  {
    *img_ptr = static_cast<double>(*int_img_ptr);
    ++ no_data;
    ++ img_ptr;
    ++ int_img_ptr;
  }
  ntuple_list out;
  out = lsd(lsd_image);
  lsd_to_lines( vLines, out );
  //lines_to_edges (vEdges, vLines);
  edgesInsideBox.viEdgeIdxes.push_back( 0 );
  for( unsigned int no_line = 0; no_line < vLines.size(); ++ no_line )
  {
    Line & line = vLines[no_line];
    for( int idx = 0; idx <= line.no_pixel; ++idx )
    {
      TooN::Vector<2> v2 = line.v2Start + static_cast<REAL_TYPE>(idx)*line.v2Direction;
      if( v2[0] > irTopLeft.x && v2[1] > irTopLeft.y && v2[0] < irButtomRight.x && v2[1] < irButtomRight.y )
      {
        edgesInsideBox.virEdges.push_back( CVD::ir(v2) );
      }
      else
      {
        if( edgesInsideBox.viEdgeIdxes.back() != (int)edgesInsideBox.virEdges.size() )
          edgesInsideBox.viEdgeIdxes.push_back( edgesInsideBox.virEdges.size() );
        do
        {
          ++idx;
          v2 = line.v2Start + static_cast<REAL_TYPE>(idx)*line.v2Direction;
        }while( idx < line.no_pixel && ( v2[0] < irTopLeft.x || v2[1] < irTopLeft.y || v2[0] > irButtomRight.x || v2[1] > irButtomRight.y) );
      }
    }
    if (edgesInsideBox.virEdges.size() > 0)
    {
      if( CVD::ir(line.v2End) != edgesInsideBox.virEdges.back() && line.v2End[0] > irTopLeft.x && line.v2End[1] > irTopLeft.y && 
        line.v2End[0] < irButtomRight.x && line.v2End[1] < irButtomRight.y )
      {
          edgesInsideBox.virEdges.push_back( CVD::ir(line.v2End) );
      }

      if( edgesInsideBox.viEdgeIdxes.back() != (int)edgesInsideBox.virEdges.size() )
      {
        edgesInsideBox.viEdgeIdxes.push_back( edgesInsideBox.virEdges.size() );
      }
    }
  }
  std::cout << "end of line entering? " << std::endl;
  free_ntuple_list(out);
/*  std::cout <<"LSD"<< std::endl;
  for( unsigned int i = 0; i < vLines.size(); ++i )
    std::cout << vLines[i].v2Start << " " << vLines[i].v2End << " " << vLines[i].no_pixel << " " << vLines[i].rSlope << std::endl; */
//  vLines.clear();
//  get_lines( vLines, edgesInsideBox );
}

// Have to remove comment from vLine.clear and get_lines RUT TEST
/*  std::cout <<"GetLines " << std::endl; 
  for( unsigned int i = 0; i < vLines.size(); ++i )
    std::cout << vLines[i].v2Start << " " << vLines[i].v2End << " " << vLines[i].no_pixel << " " << vLines[i].rSlope << std::endl; */
//  std::cout << edgesInsideBox.viEdgeIdxes.size() <<" " << edges.viEdgeIdxes.size() << std::endl;
//  exit(0);
#if 0
  mEdgeImage.fill(0);
  for( unsigned int i = 0; i < edgesInsideBox.virEdges.size(); ++i )
//    std::cout << vec(edgesInsideBox.virEdges[i]) << std::endl;
    mEdgeImage[edgesInsideBox.virEdges[i]] = 255;
  char filename[100];
  sprintf(filename,"template_%d.png",template_file_id);
  ++template_file_id;
  CVD::img_save(mEdgeImage, filename);
#endif  
  return addCodeBook( vLines, edgesInsideBox, irTopLeft, irButtomRight, v4ChainAngles, object_ID, obj_correct_ID );
}

void EdgeObjectDetection::lsd_to_lines( std::vector<Line> & vLines, ntuple_list lsd_format )
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
    if (isnan(line.rSlope))
    {
       std::cerr << "LINE SLOPE FOUND TO BE NOT A NUMBER???? " << std::endl;
       exit(-1);
    }
    normalize( line.v2Direction );
    vLines.push_back(line);
  }
}

void EdgeObjectDetection::lines_to_edges (Edges & edges, std::vector<Line> vLines)
{
   edges.virEdges.clear();
   edges.virEdges.reserve (vLines.size()*EDGELET_LENGTH_TEST);
  for( unsigned int no_line = 0; no_line < vLines.size(); ++ no_line )
  {
    Line & line = vLines[no_line];
    for( int idx = 0; idx < line.no_pixel; ++idx )
    {
      TooN::Vector<2> v2 = line.v2Start + static_cast<REAL_TYPE>(idx)*line.v2Direction;
      if( v2[0] >= 0 && v2[1] >0  && v2[0] < lsd_image->xsize && v2[1] < lsd_image->ysize )
        edges.virEdges.push_back( CVD::ir(v2) );
    }
  }

}

int EdgeObjectDetection::addCodeBook(const std::vector<Line> & vLines, const Edges vEdges, CVD::ImageRef irTopLeft, CVD::ImageRef irButtomRight, const TooN::Vector<4> & v4ChainAngles, int object_ID, int obj_correct_ID)
{
  std::cout << "trying to add codebook " << std::endl;
  if (vEdges.virEdges.size() == 0)
     return 0;
  ++miTotalObjectViews;
  unsigned int no_class = 0;
  bool found_class = false;
  for( ; no_class < mvObjDescriptorClasses.size(); ++no_class )
  {
    if( mvObjDescriptorClasses[no_class].mv4ChainAngles == v4ChainAngles )
    {
      found_class = true;
      break;
    }
  }
  
  if( !found_class ){ // no_class will be equal to the target class.
    mvObjDescriptorClasses.resize( mvObjDescriptorClasses.size() + 1 );
    std::cout << "new codebook angles added" << std::endl;
  }
  std::cout << "and here too" << std::endl;

  ObjectsTemplate & mObjectsTemplate = mvObjDescriptorClasses[no_class].mObjectsTemplate;

  std::vector<CodeBookElement> & mvCodeBook = mvObjDescriptorClasses[no_class].mvCodeBook;

  TooN::Vector<4> & mv4ChainAngles = mvObjDescriptorClasses[no_class].mv4ChainAngles;

  mv4ChainAngles = v4ChainAngles;
  if(object_ID >= mObjectsTemplate.iNoObjects)
  { object_ID = mObjectsTemplate.iNoObjects;
    ++mObjectsTemplate.iNoObjects; 
    mObjectsTemplate.vObjectElements.resize( mObjectsTemplate.iNoObjects );
    mObjectsTemplate.vObjectElements[object_ID].iObject_ID = object_ID;
  }
  else if( object_ID < -1 )
    return -1;

  std::cout << "trying to add a new view to object " << object_ID << std::endl;
  ObjectElement & obj = mObjectsTemplate.vObjectElements[object_ID];
  obj.obj_correct_ID = obj_correct_ID;
  int view_ID = obj.iNoViews;
  std::cout << "this will be view number " << view_ID << " " << vEdges.virEdges.size() << std::endl;
  ++obj.iNoViews;
  obj.vViews.resize(obj.iNoViews);
  View & view = obj.vViews[view_ID];
  view.iView_ID = view_ID;
  view.allEdges = vEdges;
  view.iObject_ID = object_ID;
  std::cout << "resized" << std::endl;

  gen_Edgelets( view.vEdgelets, view.viLineIndexes, irTopLeft, irButtomRight, vLines, EDGELET_LENGTH_TRAIN );
  std::cout << "generated lines" << std::endl;
  gen_TestDirections( view.vEdgelets, view.viLineIndexes );
  std::cout << "generated edgelets " << view.vEdgelets.size() << " with mvTestDirection of size " << mvTestDirections.size() << std::endl;
  #ifdef ECM
  gen_EdgeletCells( view.vEdgelets );
  #endif
  gen_CertainDirectionsChain( mvCodeBook, view.vEdgelets, v4ChainAngles, view_ID, object_ID, obj_correct_ID );
  std::cout << "and codebook size " << mvCodeBook.size() << std::endl;
  mbEmpty = false;
  return object_ID;
}

void EdgeObjectDetection::load_cellDirectStructure (char* filename)
{
    FILE * if_codebook = fopen( filename, "r+b" );
    //std::ifstream if_codebook( filename, std::ios::binary );

    if (if_codebook == NULL){
           std::cout << "********* cell direct structure file cannot be opened " << std::endl;
           exit(-1);
    }
    unsigned int img_read, cellNo, NoOfCells;
    unsigned int count = 0;
    for (int i = 0; i < 768; i++)
    {
       for (int j = 0; j < 63; j++)
       {  
          //cds[i].cellVals[j].cells.reserve (20000);
          cds[i].cellVals[j].cellNo = i;
          cds[i].cellVals[j].dir = j;
        //  if_codebook >> NoOfCells;
 	  img_read = fread( &NoOfCells, 1, sizeof(NoOfCells), if_codebook );
          cds[i].cellVals[j].cells.reserve(NoOfCells);
       /*   if (NoOfCells == 0){
              std::cout << "for " << i << " and " << j << " cellNo is 0 " << std::endl;
              count++;
           }*/
          for (int k = 0; k < NoOfCells; k++)
          {
             img_read = fread ( &cellNo, 1, sizeof(cellNo), if_codebook);
             cds[i].cellVals[j].cells.push_back(cellNo);
          }
       }
    }
    //if_codebook.close();
    fclose( if_codebook );
}

void EdgeObjectDetection::clear_codebooks()
{
  std::cerr << "CLEARNING CODEBOOK" << std::endl;
  for( unsigned int no_class = 0; no_class < mvObjDescriptorClasses.size(); ++no_class )
  {
     std::vector<CodeBookElement> & mvCodeBook = mvObjDescriptorClasses[no_class].mvCodeBook;
     mvCodeBook.clear();
     ObjectsTemplate & mObjectsTemplate = mvObjDescriptorClasses[no_class].mObjectsTemplate;
     mObjectsTemplate.vObjectElements.clear();
     std::cout << "resizing " << mvCodeBook.size () << " " << mObjectsTemplate.vObjectElements.size() << std::endl;
     std::vector< std::vector<int> > & mvirSecondIndexes = mvObjDescriptorClasses[no_class].mvirSecondIndexes;
     for (int i = 0; i < mvirSecondIndexes.size(); i++)
     {
       std::vector<int> & thisList = mvirSecondIndexes[i];
       thisList.clear();
     }
  }  


}

void EdgeObjectDetection::initialise_codebook()
{

  mvObjDescriptorClasses.clear();
  unsigned int class_size = 1;
  mvObjDescriptorClasses.resize( class_size );
  for( unsigned int no_class = 0; no_class < class_size; ++no_class )
  {
    std::cerr << "In chain number " << no_class << std::endl;
    ObjectsTemplate & mObjectsTemplate = mvObjDescriptorClasses[no_class].mObjectsTemplate;
    std::vector<CodeBookElement> & mvCodeBook = mvObjDescriptorClasses[no_class].mvCodeBook;
    TooN::Vector<4> & mv4ChainAngles = mvObjDescriptorClasses[no_class].mv4ChainAngles;
   /*if (no_class == 0)
    {*/
      /* mv4ChainAngles[0] = -0.80693;
       mv4ChainAngles[1] = -2.17300;
       mv4ChainAngles[2] = 2.86770;
       mv4ChainAngles[3] = 2.73730;*/
       mv4ChainAngles[0] = 0.2;
       mv4ChainAngles[1] = 1;
       mv4ChainAngles[2] = 2;
       mv4ChainAngles[3] = 3;
   /* }
    else
    {
       mv4ChainAngles[0] = 0.53106;
       mv4ChainAngles[1] = 0.00133;
       mv4ChainAngles[2] = -1.7719;
       mv4ChainAngles[3] = 0.44997;
    }*/

    std::vector<CVD::ImageRef> & mvirFirstAngleIndexes = mvObjDescriptorClasses[no_class].mvirFirstAngleIndexes;
    mvirFirstAngleIndexes.resize(64);
    std::vector< std::vector<int> > & mvirSecondIndexes = mvObjDescriptorClasses[no_class].mvirSecondIndexes;
    mvirSecondIndexes.resize(163840);
    unsigned int no_objects = 0;
    std::cerr <<"Size of Code Book : "<< mvCodeBook.size() << std::endl;
  }
  returnedEdgelets.reserve (5); // for chains of size 5
  returnedCB.reserve (1000);
}

void EdgeObjectDetection::load_codebook(std::string filename)
{
#if 0
  std::ifstream if_codebook( filename.c_str(), std::ios::binary );
  mObjectsTemplate.vObjects.clear();
  mvCodeBook.clear();
  if_codebook >> mv4ChainAngles;
  int codebook_size;
  if_codebook >> codebook_size;
  mvCodeBook.resize( codebook_size );
  for( unsigned int i = 0; i < codebook_size; ++i )
  {
    if_codebook >> mvCodeBook[i].rAngle_0;
    if_codebook >> mvCodeBook[i].rAngle_1;
    if_codebook >> mvCodeBook[i].rRelativeDist_1;
    if_codebook >> mvCodeBook[i].rAngle_2;
    if_codebook >> mvCodeBook[i].rRelativeDist_2;
    if_codebook >> mvCodeBook[i].rAngle_3;
    if_codebook >> mvCodeBook[i].rRelativeDist_3;
    if_codebook >> mvCodeBook[i].iIdx_0;
    if_codebook >> mvCodeBook[i].iIdx_1;
    if_codebook >> mvCodeBook[i].iIdx_2;
    if_codebook >> mvCodeBook[i].iIdx_3;
    if_codebook >> mvCodeBook[i].iIdx_4;
    if_codebook >> mvCodeBook[i].iView_ID;
    if_codebook >> mvCodeBook[i].iObject_ID;
  }
  if_codebook >> mObjectsTemplate.iNoObjects;
  mObjectsTemplate.vObjects.resize( mObjectsTemplate.iNoObjects );
  for( int obj_ID = 0; obj_ID < mObjectsTemplate.iNoObjects; ++obj_ID )
  {
    Object & obj = mObjectsTemplate.vObjects[obj_ID];
    obj.iObject_ID = obj_ID;
    if_codebook >> obj.iNoViews;
    obj.vViews.resize( obj.iNoViews );
    for( int view_ID = 0; view_ID < obj.iNoViews; ++view_ID )
    {
      View & view = obj.vViews[view_ID];
      view.iView_ID = view_ID;
      view.iObject_ID = obj_ID;
      std::vector<Edgelet> & vEdgelets = view.vEdgelets;
      int noEdgelets;
      if_codebook >> noEdgelets;
      vEdgelets.resize( noEdgelets );
      for( int edg_ID = 0; edg_ID < noEdgelets; ++edg_ID )
      {
        if_codebook >> vEdgelets[edg_ID].v2Pose;
        if_codebook >> vEdgelets[edg_ID].rAngle;
        if_codebook >> vEdgelets[edg_ID].rSlope;
      }
    }
  }
  if_codebook.close();
#else
  mbEmpty = false;
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
  //  TooN::Vector<4> & mv4ChainMinDirections = mvObjDescriptorClasses[no_class].mv4ChainMinDirections;
  //  TooN::Vector<4> & mv4ChainMaxDirections = mvObjDescriptorClasses[no_class].mv4ChainMaxDirections;
    for( int i = 0; i < 4; ++i ) 
    {
        img_read = fread( &mv4ChainAngles[i], 1, sizeof(mv4ChainAngles[i]), if_codebook );
     //   mv4ChainMinDirections[i] = static_cast<float>(mv4ChainAngles[i] - EPS_DIR);
     //   mv4ChainMaxDirections[i] = static_cast<float>(mv4ChainAngles[i] + EPS_DIR);
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
    std::cerr << "LOADING OBJECTS OF SIZE " << no_objects << std::endl;
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
     //   bool reading = false;
        for( unsigned int edg_ID = 0; edg_ID < no_edgelets; ++ edg_ID )
        {
          img_read = fread( &vEdgelets[edg_ID].v2Pose[0], 1, sizeof( vEdgelets[edg_ID].v2Pose[0] ), if_codebook );
          img_read = fread( &vEdgelets[edg_ID].v2Pose[1], 1, sizeof( vEdgelets[edg_ID].v2Pose[1] ), if_codebook );
      //    if (obj_ID == 0 && view_ID == 0)
      //    {
      //       std::cout << vEdgelets[edg_ID].v2Pose[0] << " " << vEdgelets[edg_ID].v2Pose[1] << std::endl;
       //      reading = true;
      //    }
          img_read = fread( &vEdgelets[edg_ID].rSlope, 1, sizeof(REAL_TYPE), if_codebook );
        }
      //  if (reading)
      //      exit(-1);
        img_read = fread( &no_edgelets2, 1, sizeof(no_edgelets2), if_codebook);
        std::vector<CVD::ImageRef> & allEdges = view.allEdges.virEdges;
        allEdges.resize (no_edgelets2);
        for (unsigned int edg_ID = 0; edg_ID < allEdges.size(); edg_ID++)
        {
           img_read = fread(&allEdges[edg_ID].y, 1, sizeof (allEdges[edg_ID].y), if_codebook);
           img_read = fread(&allEdges[edg_ID].x, 1, sizeof (allEdges[edg_ID].x), if_codebook);
        }
      }
    }
    // if a new object is added (no prior detection), it adds it at the end
    latestObjID = no_objects;
    std::cerr <<"Size of Code Book : "<< mvCodeBook.size() << std::endl;
  }
  returnedEdgelets.reserve (5); // for chains of size 5
  returnedCB.reserve (1000);
  //thisCB.reserve(10000);
  fclose( if_codebook );
#endif
#if 0 
  for( unsigned int no_class = 0; no_class < mvObjDescriptorClasses.size(); ++no_class )
  {
    ObjectsTemplate & mObjectsTemplate = mvObjDescriptorClasses[no_class].mObjectsTemplate;
    std::vector<CodeBookElement> & mvCodeBook = mvObjDescriptorClasses[no_class].mvCodeBook;
#if 1
    std::cout <<"Code Book " << std::endl;
    for( unsigned int i = 0; i < mvCodeBook.size(); ++i )
    {
      printf("%f %f %f %f %f %f %f %d %d %d %d %d %d %d\n", mvCodeBook[i].rAngle_0, mvCodeBook[i].rAngle_1, mvCodeBook[i].rRelativeDist_1, 
        mvCodeBook[i].rAngle_2, mvCodeBook[i].rRelativeDist_2, mvCodeBook[i].rAngle_3, mvCodeBook[i].rRelativeDist_3, 
        mvCodeBook[i].iIdx_0, mvCodeBook[i].iIdx_1, mvCodeBook[i].iIdx_2, mvCodeBook[i].iIdx_3, mvCodeBook[i].iIdx_4, 
        mvCodeBook[i].iView_ID, mvCodeBook[i].iObject_ID);
    }
#endif
#if 1 
    std::cout <<"Edgelets " << std::endl;
    for( unsigned int obj_ID = 0; obj_ID < mObjectsTemplate.vObjects.size(); ++ obj_ID )
    {
      Object & obj = mObjectsTemplate.vObjects[obj_ID]; 
      for( unsigned int view_ID = 0; view_ID < obj.vViews.size(); ++view_ID )
      {
        View & view = obj.vViews[view_ID];
//        printf("%d %d\n", view.iView_ID, view.iObject_ID );
        std::vector<Edgelet> & edgelets = view.vEdgelets;
        for( unsigned int edl_ID = 0; edl_ID < edgelets.size(); ++ edl_ID )
        {
          printf("%f %f %f %f\n", edgelets[edl_ID].v2Pose[0], edgelets[edl_ID].v2Pose[1], edgelets[edl_ID].rSlope, edgelets[edl_ID].rAngle);
        }
#if 1
    for( unsigned int i = 0; i < mvCodeBook.size(); ++i )
    {
      std::cout << edgelets[mvCodeBook[i].iIdx_0].v2Pose << " "
        << edgelets[mvCodeBook[i].iIdx_1].v2Pose << " "
        << edgelets[mvCodeBook[i].iIdx_2].v2Pose << " "
        << edgelets[mvCodeBook[i].iIdx_3].v2Pose << " "
        << edgelets[mvCodeBook[i].iIdx_4].v2Pose << std::endl;
    }
#endif
      }
    }
#endif    
  }
#endif
#if 0
  for( unsigned int no_class = 0; no_class < mvObjDescriptorClasses.size(); ++ no_class )
  {
    std::vector<CVD::ImageRef> & mvirFirstAngleIndexes = mvObjDescriptorClasses[no_class].mvirFirstAngleIndexes;
    REAL_TYPE rTheta = -0.5*M_PI;
    REAL_TYPE rStep = M_PI/64.f;
    for( int i = 0; i < 64; ++i )
      std::cout << rTheta + (REAL_TYPE)i*rStep << " " << mvirFirstAngleIndexes[i] << std::endl;
  }
#endif
}

void EdgeObjectDetection::load_codebook(std::string filename, int no_class )
{
  FILE * if_codebook = fopen( filename.c_str(), "r+b" );
  unsigned int class_size;
  int img_read;
  img_read = fread( &class_size, 1, sizeof(class_size), if_codebook );
  if( no_class < mvObjDescriptorClasses.size() )
  {
    ObjectsTemplate & mObjectsTemplate = mvObjDescriptorClasses[no_class].mObjectsTemplate;
    std::vector<CodeBookElement> & mvCodeBook = mvObjDescriptorClasses[no_class].mvCodeBook;
    TooN::Vector<4> & mv4ChainAngles = mvObjDescriptorClasses[no_class].mv4ChainAngles;
    for( int i = 0; i < 4; ++i ) img_read = fread( &mv4ChainAngles[i], 1, sizeof(mv4ChainAngles[i]), if_codebook );
    unsigned int codebook_size;
    img_read = fread( &codebook_size, 1, sizeof(codebook_size), if_codebook );
    unsigned int no_existing_codebooks = mvCodeBook.size();
    mvCodeBook.resize( no_existing_codebooks + codebook_size );
    int no_existing_objects = mObjectsTemplate.iNoObjects;
    for( unsigned int i = no_existing_codebooks; i < codebook_size + no_existing_codebooks; ++i )
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
      mvCodeBook[i].iObject_ID += no_existing_objects;
    }
    std::vector<CVD::ImageRef> & mvirFirstAngleIndexes = mvObjDescriptorClasses[no_class].mvirFirstAngleIndexes;
    mvirFirstAngleIndexes.resize(64);
    for( int i = 0; i < 64; ++i )
    {
      img_read = fread( &mvirFirstAngleIndexes[i].x, 1, sizeof(mvirFirstAngleIndexes[i].x), if_codebook );
      img_read = fread( &mvirFirstAngleIndexes[i].y, 1, sizeof(mvirFirstAngleIndexes[i].y), if_codebook );
    }
    unsigned int no_objects;
    img_read = fread( &no_objects, 1, sizeof(no_objects), if_codebook );
    mObjectsTemplate.iNoObjects += no_objects;
    mObjectsTemplate.vObjectElements.resize(mObjectsTemplate.iNoObjects);
    for( unsigned int obj_ID = no_existing_objects; obj_ID < no_objects + no_existing_objects; ++ obj_ID )
    {
      std::cout << "reading object " << obj_ID << std::endl;
      ObjectElement & obj = mObjectsTemplate.vObjectElements[obj_ID];
      obj.iObject_ID = obj_ID;
      img_read = fread( &obj.obj_correct_ID, 1, sizeof (obj.obj_correct_ID), if_codebook);
      img_read = fread( &obj.iNoViews, 1, sizeof(obj.iNoViews), if_codebook );
      obj.vViews.resize( obj.iNoViews );
      for( int view_ID = 0; view_ID < obj.iNoViews; ++view_ID )
      {
        std::cout << "reading view " << view_ID << std::endl;
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
        /*std::vector<int> & viLineIndexes = view.viLineIndexes;
        unsigned int no_lineIndexes;
        img_read = fread( &no_lineIndexes, 1, sizeof(no_lineIndexes), if_codebook );
        viLineIndexes.resize( no_lineIndexes );
        for( unsigned int idx = 0; idx < no_lineIndexes; ++idx )
          img_read = fread( &viLineIndexes[idx], 1, sizeof(int), if_codebook );
        std::vector<TooN::Vector<3> > & vv3VisiblePoints = view.vv3VisiblePoints;
        unsigned int no_visiblePoints;
        img_read = fread( &no_visiblePoints, 1, sizeof(no_visiblePoints), if_codebook );
        vv3VisiblePoints.resize( no_visiblePoints );
        for( unsigned int idx = 0; idx < no_visiblePoints; ++idx )
        {
          img_read = fread( &vv3VisiblePoints[idx][0], 1, sizeof(vv3VisiblePoints[idx][0]), if_codebook );
          img_read = fread( &vv3VisiblePoints[idx][1], 1, sizeof(vv3VisiblePoints[idx][1]), if_codebook );
          img_read = fread( &vv3VisiblePoints[idx][2], 1, sizeof(vv3VisiblePoints[idx][2]), if_codebook );
        }
        TooN::Matrix<3> m3R;
        TooN::Vector<3> v3;
        for( int r = 0; r < 3; ++r )
          for( int c = 0; c < 3; ++c )
            img_read = fread( &m3R[r][c], 1, sizeof(m3R[r][c]), if_codebook );
        for( int idx = 0; idx < 3; ++idx )
          img_read = fread( &v3[idx], 1, sizeof(v3[idx]), if_codebook );
        view.se3W2C.get_rotation() = m3R;
        view.se3W2C.get_translation() = v3;*/
        img_read = fread( &no_edgelets2, 1, sizeof(no_edgelets2), if_codebook);
        std::vector<CVD::ImageRef> & allEdges = view.allEdges.virEdges;
        allEdges.resize (no_edgelets2);
        for (unsigned int edg_ID = 0; edg_ID < allEdges.size(); edg_ID++){
           img_read = fread(&allEdges[edg_ID].y, 1, sizeof (allEdges[edg_ID].y), if_codebook);
           img_read = fread(&allEdges[edg_ID].x, 1, sizeof (allEdges[edg_ID].x), if_codebook);
        }

      }
    }
/*
    std::vector<CodeBookElement>::iterator bCodeBookIterator = mvCodeBook.begin();
    std::vector<CodeBookElement>::iterator eCodeBookIterator = mvCodeBook.end();
    while( bCodeBookIterator != eCodeBookIterator )
    {
      unsigned short &iObject_ID = (*bCodeBookIterator).iObject_ID;
      unsigned short &iView_ID = (*bCodeBookIterator).iView_ID;
      mObjectsTemplate.vObjectElements[iObject_ID].vViews[iView_ID].vCodeBookIterators.push_back( bCodeBookIterator );
      ++bCodeBookIterator;
    }
*/    
//    std::cout <<"Size of Code Book : "<< mvCodeBook.size() << std::endl;
    fclose( if_codebook );
    sortCodeBooks();
   // return true;
  }
  else
  {
    fclose( if_codebook );
   // return false;
  }
}

void EdgeObjectDetection::save_codebook(std::string filename)
{
#if 0
  std::ofstream of_codebook( filename.c_str(), std::ios::trunc|std::ios::binary );
  // write fixed angles used for generating code books.
  of_codebook << mv4ChainAngles << std::endl;
  // write code books
  of_codebook << mvCodeBook.size() << std::endl;
  for( unsigned int i = 0; i < mvCodeBook.size(); ++ i )
    of_codebook << mvCodeBook[i].rAngle_0 << " " 
      << mvCodeBook[i].rAngle_1 << " " << mvCodeBook[i].rRelativeDist_1 << " " 
      << mvCodeBook[i].rAngle_2 << " " << mvCodeBook[i].rRelativeDist_2 << " " 
      << mvCodeBook[i].rAngle_3 << " " << mvCodeBook[i].rRelativeDist_3 << " " 
      << mvCodeBook[i].iIdx_0 << " " << mvCodeBook[i].iIdx_1 << " " 
      << mvCodeBook[i].iIdx_2 << " " << mvCodeBook[i].iIdx_3 << " " << mvCodeBook[i].iIdx_4 << " " 
      << mvCodeBook[i].iView_ID << " " << mvCodeBook[i].iObject_ID << std::endl;
  // write objects template
  int noObjects = mObjectsTemplate.iNoObjects;
  of_codebook << noObjects << std::endl;
  for( int obj_ID = 0; obj_ID < noObjects; ++obj_ID )
  {
    Object & obj = mObjectsTemplate.vObjects[obj_ID];
    int noViews = obj.iNoViews;
    of_codebook << noViews << std::endl;
    for( int view_ID = 0; view_ID < noViews; ++view_ID )
    {
      View & view = obj.vViews[view_ID];
      std::vector<Edgelet> & vEdgelets = view.vEdgelets;
      int noEdgelets = vEdgelets.size();
      of_codebook << noEdgelets << std::endl;
      for( int edg_ID = 0; edg_ID < noEdgelets; ++edg_ID )
        of_codebook << vEdgelets[edg_ID].v2Pose << " " << vEdgelets[edg_ID].rAngle << " " << vEdgelets[edg_ID].rSlope << std::endl; 
    }
  }
  of_codebook.close();
#else
  FILE *of_codebook;
  of_codebook = fopen(filename.c_str(), "w+b");
  unsigned int class_size = mvObjDescriptorClasses.size();
  std::cerr << " class size: *********************" << class_size << std::endl;
  fwrite( &class_size, 1, sizeof(class_size), of_codebook );
  for( unsigned int no_class = 0; no_class < class_size; ++no_class )
  {
    ObjectsTemplate & mObjectsTemplate = mvObjDescriptorClasses[no_class].mObjectsTemplate;
    std::vector<CodeBookElement> & mvCodeBook = mvObjDescriptorClasses[no_class].mvCodeBook;
    TooN::Vector<4> & mv4ChainAngles = mvObjDescriptorClasses[no_class].mv4ChainAngles;

    for( int i = 0; i < 4; ++i ) fwrite( &mv4ChainAngles[i], 1, sizeof(mv4ChainAngles[i]), of_codebook );
    unsigned int no_CodeBookElements = mvCodeBook.size();
    std::cerr << "Number of codebook elements " << no_CodeBookElements << std::endl;
    fwrite( &no_CodeBookElements, 1, sizeof(no_CodeBookElements), of_codebook );
    for( unsigned int i = 0; i < no_CodeBookElements; ++i )
    {
      fwrite( &mvCodeBook[i].rAngle_0, 1, sizeof(REAL_TYPE), of_codebook );
      fwrite( &mvCodeBook[i].rAngle_1, 1, sizeof(REAL_TYPE), of_codebook );
      fwrite( &mvCodeBook[i].rRelativeDist_1, 1, sizeof(REAL_TYPE), of_codebook );
      fwrite( &mvCodeBook[i].rAngle_2, 1, sizeof(REAL_TYPE), of_codebook );
      fwrite( &mvCodeBook[i].rRelativeDist_2, 1, sizeof(REAL_TYPE), of_codebook );
      fwrite( &mvCodeBook[i].rAngle_3, 1, sizeof(REAL_TYPE), of_codebook );
      fwrite( &mvCodeBook[i].rRelativeDist_3, 1, sizeof(REAL_TYPE), of_codebook );
      fwrite( &mvCodeBook[i].iIdx_0, 1, sizeof(unsigned short), of_codebook );
      fwrite( &mvCodeBook[i].iIdx_1, 1, sizeof(unsigned short), of_codebook );
      fwrite( &mvCodeBook[i].iIdx_2, 1, sizeof(unsigned short), of_codebook );
      fwrite( &mvCodeBook[i].iIdx_3, 1, sizeof(unsigned short), of_codebook );
      fwrite( &mvCodeBook[i].iIdx_4, 1, sizeof(unsigned short), of_codebook );
      fwrite( &mvCodeBook[i].iView_ID, 1, sizeof(unsigned short), of_codebook );
      fwrite( &mvCodeBook[i].iObject_ID, 1, sizeof(unsigned short), of_codebook );
      fwrite( &mvCodeBook[i].obj_correct_ID, 1, sizeof(unsigned short), of_codebook );
    }
    std::cout << "Wait until the codebook is written fully.................. " << std::endl;
    std::vector<CVD::ImageRef> & mvirFirstAngleIndexes = mvObjDescriptorClasses[no_class].mvirFirstAngleIndexes;
    for( int i = 0; i < 64; ++i )
    {
      fwrite( &mvirFirstAngleIndexes[i].x, 1, sizeof(mvirFirstAngleIndexes[i].x), of_codebook );
      fwrite( &mvirFirstAngleIndexes[i].y, 1, sizeof(mvirFirstAngleIndexes[i].y), of_codebook );
    }
    std::vector< std::vector<int> > & mvirSecondIndexes = mvObjDescriptorClasses[no_class].mvirSecondIndexes;
    for (int i = 0; i < 163840; i++) 
    {
       int sx = mvirSecondIndexes[i].size();
       fwrite (&sx, 1, sizeof (sx), of_codebook);
       std::vector<int> &vectors = mvirSecondIndexes[i];
       for (int j = 0; j < mvirSecondIndexes[i].size(); j++)
       {
          fwrite (&vectors[j], 1, sizeof (vectors[j]), of_codebook);
       }
    }
    unsigned int no_objects = mObjectsTemplate.iNoObjects;
    fwrite( &no_objects, 1, sizeof(no_objects), of_codebook );
    for( unsigned int obj_ID = 0; obj_ID < no_objects; ++obj_ID )
    {
      ObjectElement & obj = mObjectsTemplate.vObjectElements[obj_ID];
      unsigned int correct_id = obj.obj_correct_ID;
      fwrite( &correct_id, 1, sizeof(correct_id), of_codebook );
      unsigned int no_views = obj.iNoViews;
      fwrite( &no_views, 1, sizeof(no_views), of_codebook );
      for( unsigned int view_ID = 0; view_ID < no_views; ++view_ID )
      {
        View & view = obj.vViews[view_ID];
        std::vector<Edgelet> & vEdgelets = view.vEdgelets;
        unsigned int no_edgelets = vEdgelets.size();
        fwrite( &no_edgelets, 1, sizeof(no_edgelets), of_codebook );
        for( unsigned int edg_ID = 0; edg_ID < no_edgelets; ++edg_ID )
        {
          fwrite( &vEdgelets[edg_ID].v2Pose[0], 1, sizeof(vEdgelets[edg_ID].v2Pose[0]), of_codebook );
          fwrite( &vEdgelets[edg_ID].v2Pose[1], 1, sizeof(vEdgelets[edg_ID].v2Pose[1]), of_codebook );
//          fwrite( &vEdgelets[edg_ID].rAngle, 1, sizeof(REAL_TYPE), of_codebook );
          fwrite( &vEdgelets[edg_ID].rSlope, 1, sizeof(REAL_TYPE), of_codebook );
        }
        unsigned int no_edgelets2 = view.allEdges.virEdges.size();
        fwrite (&no_edgelets2, 1, sizeof (no_edgelets2), of_codebook);
        std::vector<CVD::ImageRef> & allEdges = view.allEdges.virEdges;
        for (unsigned int edg_ID = 0; edg_ID < allEdges.size(); edg_ID++){
           fwrite( &allEdges[edg_ID].y, 1, sizeof(allEdges[edg_ID].y), of_codebook);
           fwrite( &allEdges[edg_ID].x, 1, sizeof(allEdges[edg_ID].x), of_codebook);
        }
      }
    }
#if 0
    {
    std::vector<Edgelet> & edgelets = mObjectsTemplate.vObjects[0].vViews[0].vEdgelets;
    for( unsigned int i = 0; i < mvCodeBook.size(); ++i )
    {
/*      std::cout << edgelets[mvCodeBook[i].iIdx_0].v2Pose << " "
        << edgelets[mvCodeBook[i].iIdx_1].v2Pose << " "
        << edgelets[mvCodeBook[i].iIdx_2].v2Pose << " "
        << edgelets[mvCodeBook[i].iIdx_3].v2Pose << " "
        << edgelets[mvCodeBook[i].iIdx_4].v2Pose << std::endl;*/
      std::cout << mvCodeBook[i].iIdx_0 << " "
        << mvCodeBook[i].iIdx_1 << " "
        << mvCodeBook[i].iIdx_2 << " "
        << mvCodeBook[i].iIdx_3 << " "
        << mvCodeBook[i].iIdx_4 << std::endl;
    }
    }
#endif
#if 0 
    for( unsigned int no_class = 0; no_class < mvObjDescriptorClasses.size(); ++no_class )
    {
      ObjectsTemplate & mObjectsTemplate = mvObjDescriptorClasses[no_class].mObjectsTemplate;
      std::vector<CodeBookElement> & mvCodeBook = mvObjDescriptorClasses[no_class].mvCodeBook;
  #if 1
      std::cout <<"Code Book " << std::endl;
      for( unsigned int i = 0; i < mvCodeBook.size(); ++i )
      {
        printf("%f %f %f %f %f %f %f %d %d %d %d %d %d %d\n", mvCodeBook[i].rAngle_0, mvCodeBook[i].rAngle_1, mvCodeBook[i].rRelativeDist_1, 
          mvCodeBook[i].rAngle_2, mvCodeBook[i].rRelativeDist_2, mvCodeBook[i].rAngle_3, mvCodeBook[i].rRelativeDist_3, 
          mvCodeBook[i].iIdx_0, mvCodeBook[i].iIdx_1, mvCodeBook[i].iIdx_2, mvCodeBook[i].iIdx_3, mvCodeBook[i].iIdx_4, 
          mvCodeBook[i].iView_ID, mvCodeBook[i].iObject_ID);
      }
  #endif
    }
#endif
  }
  fclose(of_codebook);
#endif
}

void EdgeObjectDetection::project_edgelets( std::vector<Edgelet> & vProjectedEdgelets, const std::vector<Edgelet> & vEdgelets, const TooN::Matrix<3> & h, bool bWrite )
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
#if 0  
  if( bWrite )
  {
    mEdgeImage.fill(0);

    char filename[100];
    sprintf(filename,"%s/projected_edges_%d.%d.%d.xvy", folderName, master_file_id, projected_file_id, obj_correct_ID);
    std::ofstream out_xvy(filename);
    for( unsigned int i = 0; i < vProjectedEdgelets.size(); ++i )
    {
      mEdgeImage[ CVD::ir( vProjectedEdgelets[i].v2Pose )] = 255;
      out_xvy << vProjectedEdgelets[i].v2Pose << std::endl;
    }
    out_xvy.close();
    sprintf(filename,"%s/projected_edges_%d.%d.png", folderName, master_file_id, projected_file_id);
    ++projected_file_id;
    CVD::img_save(mEdgeImage, filename);
  }
#endif
}

void EdgeObjectDetection::generateOfflineCellDirectStruct (int h, int w, int cellH, int cellW, REAL_TYPE cellDirInterval)
{
    const int THIS_FIRST_DIST = 500;
    int noCellBins = 2*3.14/cellDirInterval + 1;
    int noHBins = h / cellH;
    int noWBins = w / cellW;
    bool checkedCells [768];
    unsigned int cellNo, sum=0;
    char filename [100];
    sprintf(filename, "cellOfflineStruct3_%d_%d.dat", h, w);
    //std::ofstream of_cds (filename, std::ios::binary);
    FILE *of_cds;
    of_cds = fopen(filename, "w+b");
    std::cout << "trying to build for " << noHBins << " " << noWBins << " " << noCellBins << std::endl;
    int cellDir;
    for (int i = 0; i < cellCount; i++)
    {
       for (int k = 0; k < noCellBins; k++)
       {
          sum = 0;
          for (int j = 0; j < cellCount; j++)
          {    
             cellDir = calculateTwoCellsDirection (i,j,0.1);
             if (cellDir == k){
                checkedCells[j] = 1;
                sum++;
             }
             else
                checkedCells[j] = 0;
          }
          fwrite (&sum, 1, sizeof(sum), of_cds);
          std::cout << sum << " ";
          for (int l = 0; l < cellCount; l++)
             if (checkedCells[l])
                fwrite (&l, 1, sizeof(l), of_cds);
       }
        std::cout << std::endl;
       fflush(of_cds);
    }
/*
    for (int j = 0; j < noHBins; j++)
          for (int i = 0; i < noWBins; i++)
          {
              std::cout << j << " " << i << std::endl;
              for (int k = 0; k < noCellBins; k++)
              {
                  int centreY = j*cellH+cellH/2;
                  int centreX = i*cellW+cellW/2;
                  sum = 0;
                  REAL_TYPE dir = k * cellDirInterval;
                  std::vector<Edgelet> pixelsWithin = findAllPixelsWithin (h, w, centreY, centreX, dir, THIS_FIRST_DIST);
                  for (int l = 0; l < pixelsWithin.size(); l++)
                  {
                     cellNo = pixelsWithin[l].v2Pose[1]/cellH * cellRow + pixelsWithin[l].v2Pose[0]/cellW;
                     checkedCells[cellNo] = 1;
                  }
                  
                  for (unsigned int l = 0; l < 768; l++)
                      if (checkedCells[l])
                          sum++;
                  fwrite (&sum, 1, sizeof(sum), of_cds);
                  std::cout << sum << " ";
                  for (unsigned int l = 0; l < 768; l++)
                  {
                      if (checkedCells[l])
                      {      
                           //of_cds << l << " ";
                           fwrite (&l, 1, sizeof(l), of_cds);
                      }
                  }
              }
              std::cout << std::endl;
              fflush(of_cds);
          }*/
    fclose(of_cds);
//    of_cds.close();
}

std::vector<Edgelet> EdgeObjectDetection::findAllPixelsWithin (int h, int w, int y, int x, REAL_TYPE dir, const int first_dist)
{
    std::vector<Edgelet> pixelsWithin;
    TooN::Vector<2> centre;
    REAL_TYPE dist, angle;
    centre[0] = x;
    centre[1] = y;
    pixelsWithin.reserve (h*w/2);
    for (int i = 0; i < w; i++)
       for (int j = 0; j < h; j++)
       {
         TooN::Vector<2> thisPixel;
         thisPixel[0] = i;
         thisPixel[1] = j;
         TooN::Vector<2>  v2Del = centre - thisPixel;
         dist = norm (v2Del);
         angle = fabs(atan2(v2Del[1], v2Del[0]) - dir);
         if (dist < first_dist && angle < EPS_DIR_TEST)
         {
            Edgelet ed;
            ed.v2Pose[0] = i;
            ed.v2Pose[1] = j;
            pixelsWithin.push_back (ed);
         }
       }
    return pixelsWithin;
}

void EdgeObjectDetection::project_edges( Edges & vProjectedEdges, std::vector<TooN::Vector<2> >& vv2Pnts, const Edges & vEdges, const TooN::Matrix<3> & h, bool bWrite )
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
    v2Pose[0] = v2Pose[0]*imgSizeRatio;
    v2Pose[1] = v2Pose[1]*imgSizeRatio;
    vv2Pnts.push_back (v2Pose);
  }
#if 0  
  if( bWrite )
  {
    mEdgeImage.fill(0);

    char filename[100];
    sprintf(filename,"%s/projected_all_edges_%d.%d.%d.xvy", folderName, master_file_id, projected_file_id, obj_correct_ID);
    std::ofstream out_xvy(filename);
    TooN::Vector<2> v2Pose;
    for( unsigned int i = 0; i < vProjectedEdgelets.size(); ++i )
    {
      v2Pose[0] = vProjectedEdgelets[i].x;
      v2Pose[1] = vProjectedEdgelets[i].y;
      mEdgeImage[CVD::ir(v2Pose)] = 255;
      out_xvy << v2Pose << std::endl;
    }
    out_xvy.close();
    sprintf(filename,"%s/projected_all_edges_%d.%d.png", folderName, master_file_id, projected_file_id);
    ++projected_file_id;
    CVD::img_save(mEdgeImage, filename);
  }
#endif
}

void EdgeObjectDetection::getChain_level1 (std::vector<int> &returnedEdgelets, std::vector<std::vector<CodeBookElement>::iterator> &returnedCB, const std::vector<Edgelet> testEdgelets, std::vector<TestDirection> Q, int edgeNo, TooN::Vector<4> mvTestDirections, EdgeCellMatrix* edgeletImageMatrix, REAL_TYPE prevDir, REAL_TYPE prevDist, std::vector<std::vector<CodeBookElement>::iterator> mvCBIterator, REAL_TYPE minDir, REAL_TYPE maxDir)
{
    REAL_TYPE angle = mvTestDirections[1]+prevDir + 2*pi;
    if (angle > 2*pi)
       angle = angle - 2*pi;
    REAL_TYPE direction = angle;
    std::cout << "entering to find cells in direction " << direction << " of " << angle << " from " << mvTestDirections[1] << " and prev direction " << prevDir << std::endl;
    std::vector<unsigned int> otherEdgelets = getEdgeletsInDirection(testEdgelets, edgeNo, edgeletImageMatrix, 10, 10, direction);
    if (otherEdgelets.size() == 0)
    {
       returnedCB = mvCBIterator;
       return;
    }

    int beginInd = mvirDirectionIndexes[testEdgelets[edgeNo].v2Pose[0]].x;
    if (beginInd == -1 || otherEdgelets.size() == 0)
    {
       returnedCB = mvCBIterator;
       std::cout << "invalid beginIndex " << std::endl;
       return;
    }
    int endInd = mvirDirectionIndexes[testEdgelets[edgeNo].v2Pose[0]].y;
    int addVal = edgeNo*(testEdgelets.size()-1);
    std::vector<int> vRandomSamples2(otherEdgelets.size());
    for (unsigned int i = 0; i < otherEdgelets.size(); i++)
        vRandomSamples2[i] = i;
    std::random_shuffle (vRandomSamples2.begin(), vRandomSamples2.end());
    for (int otherE = 0; otherE < vRandomSamples2.size(); otherE++)
    {
       REAL_TYPE thisAngle, thisDir, thisDist, dir2;
       int indexVal;
       int edgeNo2 = otherEdgelets[vRandomSamples2[otherE]];
       if (edgeNo2 < edgeNo)
          indexVal = addVal + edgeNo2;
       else
          indexVal = addVal + edgeNo2 -1;

       const TestDirection &dPair = Q[indexVal];
       thisAngle = dPair.rAngle;
       thisDir = dPair.rDir_2;
       thisDist = dPair.rDist;
       dir2 = thisDir + prevDir;
       if (dir2 > pi)
           dir2 = dir2 - 2*pi; // CHECK
     

       if (dir2 < minDir || dir2 > maxDir)
       { 
          std::cout << "the unmatched directions are " << edgeNo << " " << edgeNo2 << " " << dir2 << " " << minDir << " " << maxDir << std::endl; 
          std::cout << "this can be explained by (" << testEdgelets[edgeNo].v2Pose[0] << ", " << testEdgelets[edgeNo].v2Pose[1] << ") and (" << testEdgelets[edgeNo2].v2Pose[0] << ", " << testEdgelets[edgeNo2].v2Pose[1] << ") does it work?" << std::endl;
          continue; 
       }
       
      // std::cout << "matched correct direction " << edgeNo2 << std::endl;
       std::vector<std::vector<CodeBookElement>::iterator> passCodeBook;
       passCodeBook.reserve( mvCBIterator.size() );
       for( unsigned int pass_idx = 0; pass_idx < mvCBIterator.size(); ++pass_idx )
         if( fabs( mvCBIterator[pass_idx]->rAngle_1 - thisAngle ) < EPS_ANGLE && fabs( mvCBIterator[pass_idx]->rRelativeDist_1 - thisDist ) < EPS_DIST)
         {
          passCodeBook.push_back(mvCBIterator[pass_idx]);  
         }
      
       if( passCodeBook.size() == 0 ) 
       {
       //  std::cout << "No codebook matched " << std::endl;
         continue;
       }
         
       std::vector<int> returnedEdgelets2;
       returnedEdgelets2.reserve (3);
       std::cout << "entering level 2 with edgeNo " << edgeNo2 << " and code book size " << passCodeBook.size() << std::endl;
       getChain_level2 (returnedEdgelets2, returnedCB, testEdgelets, Q, edgeNo2, mvTestDirections, edgeletImageMatrix, dir2, thisDist, passCodeBook, minDir, maxDir);
       for (int j = 0; j < returnedEdgelets2.size(); j++)
          std::cout << returnedEdgelets2[j] << " ";
       std::cout << std::endl;

       if (returnedEdgelets2.size() > 0)
       {
          returnedEdgelets.push_back (edgeNo2);
          for (int j = 0; j < returnedEdgelets2.size(); j++)
             returnedEdgelets.push_back (returnedEdgelets2[j]);
        //  returnedCB = passCodeBook;
          break;
       }
   }
}

void EdgeObjectDetection::getChain_level2 (std::vector<int> &returnedEdgelets, std::vector<std::vector<CodeBookElement>::iterator> &returnedCB, const std::vector<Edgelet> testEdgelets, std::vector<TestDirection> Q, int edgeNo, TooN::Vector<4> mvTestDirections, EdgeCellMatrix* edgeletImageMatrix, REAL_TYPE prevDir, REAL_TYPE prevDist, std::vector<std::vector<CodeBookElement>::iterator> mvCBIterator, REAL_TYPE minDir, REAL_TYPE maxDir)
{
    REAL_TYPE angle = mvTestDirections[2]+prevDir + 2*pi;
    if (angle > 2*pi)
       angle = angle - 2*pi;
    REAL_TYPE direction = angle;
    std::vector<unsigned int> otherEdgelets = getEdgeletsInDirection(testEdgelets, edgeNo, edgeletImageMatrix, 10, 10, direction);
    if (otherEdgelets.size() == 0)
    {
       returnedCB = mvCBIterator;
     //  std::cout << "no edgelets in this direction " << direction << std::endl;
       return;
    }
    std::cout << "No of edgelets in this direction " << otherEdgelets.size() << std::endl;

    int beginInd = mvirDirectionIndexes[testEdgelets[edgeNo].v2Pose[0]].x;
    if (beginInd == -1 || otherEdgelets.size() == 0)
    {
       returnedCB = mvCBIterator;
       std::cout << "beginInd invalid or " << otherEdgelets.size() << std::endl;
       return;
    }
    int endInd = mvirDirectionIndexes[testEdgelets[edgeNo].v2Pose[0]].y;
    int addVal = edgeNo*(testEdgelets.size()-1);
    std::vector<int> vRandomSamples2(otherEdgelets.size());
    for (unsigned int i = 0; i < otherEdgelets.size(); i++)
        vRandomSamples2[i] = i;
    std::random_shuffle (vRandomSamples2.begin(), vRandomSamples2.end());
    for (int otherE = 0; otherE < vRandomSamples2.size(); otherE++)
    {
       REAL_TYPE thisAngle, thisDir, thisDist, dir2;
       int indexVal;
       int edgeNo2 = otherEdgelets[vRandomSamples2[otherE]];
       if (edgeNo2 < edgeNo)
          indexVal = addVal + edgeNo2;
       else
          indexVal = addVal + edgeNo2 -1;

       const TestDirection &dPair = Q[indexVal];
       thisAngle = dPair.rAngle;
       thisDir = dPair.rDir_1;
       thisDist = dPair.rDist;
       dir2 = thisDir + prevDir;
       if (dir2 > pi)
           dir2 = dir2 - 2*pi; // CHECK
     

       if (dir2 < minDir || dir2 > maxDir)
       {  
        //  std::cout << "invalid otherE " << otherEdgelets[vRandomSamples2[otherE]] << std::endl;
          continue; 
       }

       std::cout << "valid otherE " << otherEdgelets[vRandomSamples2[otherE]] << std::endl;       

       std::vector<std::vector<CodeBookElement>::iterator> passCodeBook;
       passCodeBook.reserve( mvCBIterator.size() );
       for( unsigned int pass_idx = 0; pass_idx < mvCBIterator.size(); ++pass_idx )
         if( fabs( mvCBIterator[pass_idx]->rAngle_2 - thisAngle ) < EPS_ANGLE && fabs( mvCBIterator[pass_idx]->rRelativeDist_2 - thisDist ) < EPS_DIST)
         {
          passCodeBook.push_back(mvCBIterator[pass_idx]);  
         }
      
       if( passCodeBook.size() == 0 ) 
       {
          std::cout << "No matching codebook stop 2 " << std::endl;
          continue;
       }
         
       std::vector<int> returnedEdgelets2;
       returnedEdgelets2.reserve (2);
       std::cout << "entering level 3 with edgeNo " << edgeNo2 << std::endl;
       getChain_level3 (returnedEdgelets2, returnedCB, testEdgelets, Q, edgeNo2, mvTestDirections, edgeletImageMatrix, dir2, thisDist, passCodeBook, minDir, maxDir);
       for (int j = 0; j < returnedEdgelets2.size(); j++)
          std::cout << returnedEdgelets2[j] << " ";
       std::cout << std::endl;

       if (returnedEdgelets2.size() > 0)
       {
          returnedEdgelets.push_back (edgeNo2);
          for (int j = 0; j < returnedEdgelets2.size(); j++)
             returnedEdgelets.push_back (returnedEdgelets2[j]);
        //  returnedCB = passCodeBook;
          break;
       }
   }
   std::cout << "finished without matching any codebooks" << std::endl;
}

void EdgeObjectDetection::getChain_level3 (std::vector<int> &returnedEdgelets, std::vector<std::vector<CodeBookElement>::iterator> &returnedCB, const std::vector<Edgelet> testEdgelets, std::vector<TestDirection> Q, int edgeNo, TooN::Vector<4> mvTestDirections, EdgeCellMatrix* edgeletImageMatrix, REAL_TYPE prevDir, REAL_TYPE prevDist, std::vector<std::vector<CodeBookElement>::iterator> mvCBIterator, REAL_TYPE minDir, REAL_TYPE maxDir)
{
    REAL_TYPE angle = mvTestDirections[3]+prevDir + 2*pi;
    if (angle > 2*pi)
       angle = angle - 2*pi;
    REAL_TYPE direction = angle;
    std::vector<unsigned int> otherEdgelets = getEdgeletsInDirection(testEdgelets, edgeNo, edgeletImageMatrix, 10, 10, direction);
    if (otherEdgelets.size() == 0)
    {
       returnedCB = mvCBIterator;
       return;
    }

    int beginInd = mvirDirectionIndexes[testEdgelets[edgeNo].v2Pose[0]].x;
    if (beginInd == -1 || otherEdgelets.size() == 0)
    {
       returnedCB = mvCBIterator;
       return;
    }
    int endInd = mvirDirectionIndexes[testEdgelets[edgeNo].v2Pose[0]].y;
    int addVal = edgeNo*(testEdgelets.size()-1);
    std::vector<int> vRandomSamples2(otherEdgelets.size());
    for (unsigned int i = 0; i < otherEdgelets.size(); i++)
        vRandomSamples2[i] = i;
    std::random_shuffle (vRandomSamples2.begin(), vRandomSamples2.end());
    for (int otherE = 0; otherE < vRandomSamples2.size(); otherE++)
    {
       REAL_TYPE thisAngle, thisDir, thisDist, dir2;
       int indexVal;
       int edgeNo2 = otherEdgelets[vRandomSamples2[otherE]];
       if (edgeNo2 < edgeNo)
          indexVal = addVal + edgeNo2;
       else
          indexVal = addVal + edgeNo2 -1;

       const TestDirection &dPair = Q[indexVal];
       thisAngle = dPair.rAngle;
       thisDir = dPair.rDir_1;
       thisDist = dPair.rDist;
       dir2 = thisDir + prevDir;
       if (dir2 > pi)
           dir2 = dir2 - 2*pi; // CHECK
     

       if (dir2 < minDir || dir2 > maxDir)
       {  
          continue; 
       }
       
       std::vector<std::vector<CodeBookElement>::iterator> passCodeBook;
       passCodeBook.reserve( mvCBIterator.size() );
       for( unsigned int pass_idx = 0; pass_idx < mvCBIterator.size(); ++pass_idx )
         if( fabs( mvCBIterator[pass_idx]->rAngle_3 - thisAngle ) < EPS_ANGLE && fabs( mvCBIterator[pass_idx]->rRelativeDist_3 - thisDist ) < EPS_DIST)
         {
          passCodeBook.push_back(mvCBIterator[pass_idx]);  
         }
      
       if( passCodeBook.size() == 0 ) 
       {
         continue;
       }
         
       std::vector<int> returnedEdgelets2;
       returnedEdgelets2.reserve (1);
      // getChain_level2 (returnedEdgelets2, returnedCB, testEdgelets, Q, edgeNo2, mvTestDirections, edgeletImageMatrix, dir2, thisDist, passCodeBook, minDir, maxDir);
      // for (int j = 0; j < returnedEdgelets2.size(); j++)
      //    std::cout << returnedEdgelets2[j] << " ";
      // std::cout << std::endl;

      // if (returnedEdgelets2.size() > 0)
      // {
          std::cout << "found the ending one edgeNo " << edgeNo2 << std::endl;
          returnedEdgelets.push_back (edgeNo2);
          returnedCB = passCodeBook;
          break;
      // }
   }
}

void EdgeObjectDetection::getChainInDirection_withCodebook (std::vector<int> &returnedEdgelets, std::vector<std::vector<CodeBookElement>::iterator> &returnedCB, const std::vector<Edgelet> testEdgelets, std::vector<TestDirection> Q, int edgeNo, TooN::Vector<4> mvTestDirections, unsigned int chainNo, EdgeCellMatrix* edgeletImageMatrix, REAL_TYPE prevDir, REAL_TYPE prevDist, std::vector<std::vector<CodeBookElement>::iterator> mvCBIterator, REAL_TYPE minDir, REAL_TYPE maxDir)
{  
   //std::cout << "entered recursion ************** chainNo " << chainNo << std::endl;
     if (mvTestDirections.size() == chainNo)
   {
       returnedCB = mvCBIterator;
       return;
   }
   REAL_TYPE angle = mvTestDirections[chainNo] + prevDir + 2*pi;
   if (angle > 2*pi)
      angle = angle - 2*pi;
   int direction = angle;
   std::vector<unsigned int> otherEdgelets = getEdgeletsInDirection(testEdgelets, edgeNo, edgeletImageMatrix, 10, 10, direction);
   int beginSndIndex = mvirDirectionIndexes[testEdgelets[edgeNo].v2Pose[0]].x;
   if( beginSndIndex == -1 || otherEdgelets.size() == 0) 
   {
      returnedCB = mvCBIterator;
      if (chainNo == 1)
         //std::cout << "No other edgelets from this angle " << angle << std::endl;
      return;
   }
//std::cout << "entering recursion from edge no " << edgeNo << " with angle " << mvTestDirections[chainNo] << ", prevDir " << prevDir << ", resultAngle " << angle << " and code elements size " << mvCBIterator.size() << std::endl;
   //std::cout << "Edgelets in this direction returned are " << otherEdgelets.size() << std::endl;
   int endSndIndex = mvirDirectionIndexes[testEdgelets[edgeNo].v2Pose[0]].y;
   int addVal = edgeNo * (testEdgelets.size()-1);
      std::vector<int> vRandomSamples2(otherEdgelets.size());
   for( unsigned int i = 0; i < otherEdgelets.size(); ++ i ) vRandomSamples2[i]=i;
   std::random_shuffle( vRandomSamples2.begin(), vRandomSamples2.end() );

   for(int otherE = 0; otherE < vRandomSamples2.size(); ++otherE)
   {
       REAL_TYPE thisAngle, thisDir, thisDist, dir2;
     
      int indexVal;
      if (otherEdgelets[vRandomSamples2[otherE]] < edgeNo)
      {
         indexVal = addVal + otherEdgelets[vRandomSamples2[otherE]];
      }
      else
      {
         indexVal = addVal + otherEdgelets[vRandomSamples2[otherE]] - 1;
      }
      const TestDirection & dPair = Q[indexVal];
      if (dPair.iIdx_0 != edgeNo || dPair.iIdx_1 != otherEdgelets[vRandomSamples2[otherE]])
      {
          std::cout << "I made a mistake in the mapping to vTestDirections " << std::endl;
          std::cout << "what it accesses is (" << dPair.iIdx_0 << ", " << dPair.iIdx_1 << ") - while it was looking for (" << edgeNo << ", " << otherEdgelets[vRandomSamples2[otherE]] << ")" << std::endl;
          exit(-1);
      }
      thisAngle =  dPair.rAngle;
      thisDir = dPair.rDir_1;
      thisDist = dPair.rDist;
      dir2 = thisDir + prevDir;
      if (dir2 > pi)
        dir2 = dir2 - 2*pi;
      std::vector<std::vector<CodeBookElement>::iterator> passCodeBook;
      passCodeBook.reserve( mvCBIterator.size() );

       if (dir2 > minDir && dir2 < maxDir)
       {
          // find the codebooks related to this
          if (chainNo == 0)
          { 
             for( unsigned int pass_idx = 0; pass_idx < mvCBIterator.size(); ++pass_idx )
               if( fabs( mvCBIterator[pass_idx]->rAngle_1 - thisAngle ) < EPS_ANGLE && fabs( mvCBIterator[pass_idx]->rRelativeDist_1 - thisDist ) < EPS_DIST){
                 passCodeBook.push_back(mvCBIterator[pass_idx]);  
               }
          }
          else if (chainNo == 1)
          {  
             for( unsigned int pass_idx = 0; pass_idx < mvCBIterator.size(); ++pass_idx )
               if( fabs( mvCBIterator[pass_idx]->rAngle_2 - thisAngle ) < EPS_ANGLE && fabs( mvCBIterator[pass_idx]->rRelativeDist_2 - thisDist ) < EPS_DIST){
                 passCodeBook.push_back(mvCBIterator[pass_idx]);
               }
          }
          else
          {  
            for( unsigned int pass_idx = 0; pass_idx < mvCBIterator.size(); ++pass_idx )
               if( fabs( mvCBIterator[pass_idx]->rAngle_3 - thisAngle ) < EPS_ANGLE && fabs( mvCBIterator[pass_idx]->rRelativeDist_3 - thisDist ) < EPS_DIST){
                 passCodeBook.push_back(mvCBIterator[pass_idx]);
               }
          }
       }
          if( passCodeBook.size() == 0 ) {
            if (chainNo == 1)
               std::cout << "pass code size is 0" << std::endl;
            continue;
          }
          std::vector<int> returnedEdgelets2;
          returnedEdgelets2.reserve (mvTestDirections.size());
 
               std::cout << "Before entering the recursion of chainNo " << (chainNo+1) << ", codebook is " << edgeNo << " " << otherEdgelets[vRandomSamples2[otherE]]   << " prevDir " << prevDir << " thisDir " << thisDir << " direction2 " << dir2 << " size of mvTestDirections " << chainNo << " and passed codebook size = " << passCodeBook.size() << std::endl;
               getChainInDirection_withCodebook (returnedEdgelets2, returnedCB, testEdgelets, Q, otherEdgelets[vRandomSamples2[otherE]], mvTestDirections, chainNo+1, edgeletImageMatrix, dir2, thisDist, passCodeBook, minDir, maxDir);
               std::cout << "after entering recursion of chainNo << " << (chainNo+1) << " codebook is " << returnedCB.size() << " returnedEdgelets2.size() = " << returnedEdgelets2.size();
               for (int j = 0; j < returnedEdgelets2.size(); j++)
                     std::cout << returnedEdgelets2[j] << " ";
               std::cout << std::endl;

          if (chainNo == (mvTestDirections.size()-1))
          {  returnedEdgelets.push_back (otherEdgelets[vRandomSamples2[otherE]]);
             for (int j = 0; j < returnedEdgelets2.size(); j++){
                returnedEdgelets.push_back(returnedEdgelets2[j]);
             }
             std::cout << "RESULTING CHAIN ////////////// " ;
             for (int j = 0; j < returnedEdgelets.size(); j++){
                std::cout << returnedEdgelets[j] << " ";
             }
             std::cout << std::endl;
          }
   }

}

int EdgeObjectDetection::getCellNo (Edgelet e)
{
   unsigned int cellH = 10;
   unsigned int cellW = 10;
   //ASSUMES FIXED CELL H AND CELL W
   unsigned int thisCellH = (e.v2Pose[0] + cellH/2)/cellH;
  if (thisCellH > 31)
      thisCellH = 31;
   unsigned int thisCellW = (e.v2Pose[1] + cellW/2)/cellW;
   if (thisCellW > 23)
     thisCellW = 23;
   unsigned int cellNo = thisCellW * cellRow + thisCellH;
   return cellNo;
}

EdgeCellMatrix* EdgeObjectDetection::getEdgeletImageMatrix (const std::vector<Edgelet> &testEdgelets, int h, int w, int cellH, int cellW)
{   
    for (int i = 0; i < cellCount; i++)
    {
         ecm[i].count = 0;
         ecm[i].cellEdgelets.clear();
    }

    for (int i = 0; i < testEdgelets.size(); i++)
    {
       unsigned int cellNo = getCellNo (testEdgelets[i]);
       if (cellNo > 768)
       {
           exit(-1);
       }
       ecm[cellNo].count++;
       ecm[cellNo].cellEdgelets.push_back (i);
    }
    return ecm;
}

std::vector<unsigned int> EdgeObjectDetection::getEdgeletsInDirection (const std::vector<Edgelet> testEdgelets, int edgeNo, EdgeCellMatrix* edgeletImageMatrix, int cellH, int cellW, REAL_TYPE dir)
{
   int thisCellH = (testEdgelets[edgeNo].v2Pose[1]+cellH/2)/cellH;
   int thisCellW = (testEdgelets[edgeNo].v2Pose[0]+cellW/2)/cellW;
   unsigned int cellNo = thisCellH * cellRow + thisCellW;
   unsigned int dirCell = dir/0.1;
   std::vector<unsigned int> otherCellsInDirection = cds[cellNo].cellVals[dirCell].cells;
   //std::cout << "No of cells to search within " << cds[cellNo].cellVals[dirCell].cells.size() << " for " << cellNo << " for edge " << edgeNo << " in position (" << testEdgelets[edgeNo].v2Pose[0] << ", " << testEdgelets[edgeNo].v2Pose[1] << ") mapped to cell (" << thisCellH << ", " << thisCellW << ") direction " << dirCell << std::endl;
   std::vector<unsigned int> otherEdgelets;
   int newEdgeNo;
   otherEdgelets.reserve (testEdgelets.size());
   for (int i = 0; i < cds[cellNo].cellVals[dirCell].cells.size(); i++)
   {
       if (edgeletImageMatrix[cds[cellNo].cellVals[dirCell].cells[i]].count > 0)
       {
          for (int j = 0; j < edgeletImageMatrix[cds[cellNo].cellVals[dirCell].cells[i]].count; j++)
          {   newEdgeNo = edgeletImageMatrix[cds[cellNo].cellVals[dirCell].cells[i]].cellEdgelets[j];
              if (newEdgeNo != edgeNo) 
              {
                 otherEdgelets.push_back (newEdgeNo);
              }
          }
       }
   }
   //std::cout << "No of edgelets in those cells " << otherEdgelets.size () << std::endl;
   return otherEdgelets;
}

REAL_TYPE EdgeObjectDetection::findClosestEdges( std::vector<MatchedPair> & vMatchedPairs, std::vector<int> & fullIndexes, std::vector<int> & viIndexes, const std::vector<Edgelet> & vMoveEdgelets, const std::vector<Edgelet> & vTargetEdgelets, REAL_TYPE rErr1, REAL_TYPE rErr2, REAL_TYPE rAngleThreshold )
{
  REAL_TYPE rLamda = 2.f;
  REAL_TYPE rError = 0.f;
  vMatchedPairs.clear();
  vMatchedPairs.reserve( vMoveEdgelets.size() );
  viIndexes.clear();
  viIndexes.reserve( vMoveEdgelets.size() );
  fullIndexes.clear();
  fullIndexes.reserve( vMoveEdgelets.size() );
  for( unsigned int iMove = 0; iMove < vMoveEdgelets.size(); ++iMove )
  {
    REAL_TYPE rMinDistance = rErr1;
    int iMinIndex = -1;
    REAL_TYPE rAngle;
    MatchedPair mp;
    for( unsigned int iTarget = 0; iTarget < vTargetEdgelets.size(); ++iTarget )
    {
      const Edgelet & eTarget = vTargetEdgelets[iTarget];
      const Edgelet & eMove = vMoveEdgelets[iMove];
      rAngle = fabs( atan( (eTarget.rSlope - eMove.rSlope)/(1.f + eTarget.rSlope*eMove.rSlope) ) );
      if( rAngle < rAngleThreshold )
      {
        REAL_TYPE rDistance = norm( eTarget.v2Pose - eMove.v2Pose ) + rLamda * rAngle; 
        if( rDistance < rMinDistance )
        {
          rMinDistance = rDistance;
          iMinIndex = iTarget;
        }
      }
    }
    if( iMinIndex != -1 )
    {
      rError += rMinDistance;
      mp.v2FstView = vMoveEdgelets[iMove].v2Pose;
      mp.v2SndView = vTargetEdgelets[iMinIndex].v2Pose;
      vMatchedPairs.push_back(mp);
      viIndexes.push_back(iMinIndex);
    }
    else
    {
      rError += rErr2;
    }
    fullIndexes.push_back(iMinIndex);
  }
  rError /= vMoveEdgelets.size();
  return rError;
}

REAL_TYPE EdgeObjectDetection::findClosestEdges_NP( std::vector<MatchedPair> & vMatchedPairs, std::vector<int> & fullIndexes, std::vector<int> & viIndexes, const std::vector<Edgelet> & vMoveEdgelets, const std::vector<Edgelet> & vTargetEdgelets, REAL_TYPE rErr1, REAL_TYPE rErr2, REAL_TYPE rAngleThreshold )
{
  REAL_TYPE rLamda = 2.f;
  REAL_TYPE rError = 0.f;
  int concealedEdgesCount = 0;
  vMatchedPairs.clear();
  vMatchedPairs.reserve( vMoveEdgelets.size() );
  viIndexes.clear();
  viIndexes.reserve( vMoveEdgelets.size() );
  fullIndexes.clear();
  fullIndexes.reserve( vMoveEdgelets.size() );
  for( unsigned int iMove = 0; iMove < vMoveEdgelets.size(); ++iMove )
  {
    REAL_TYPE rMinDistance = rErr1;
    int iMinIndex = -1;
    REAL_TYPE rAngle;
    MatchedPair mp;
    for( unsigned int iTarget = 0; iTarget < vTargetEdgelets.size(); ++iTarget )
    {
      const Edgelet & eTarget = vTargetEdgelets[iTarget];
      const Edgelet & eMove = vMoveEdgelets[iMove];
      rAngle = fabs( atan( (eTarget.rSlope - eMove.rSlope)/(1.f + eTarget.rSlope*eMove.rSlope) ) );
      if( rAngle < rAngleThreshold )
      {
        REAL_TYPE rDistance = norm( eTarget.v2Pose - eMove.v2Pose ) + rLamda * rAngle; 
        if( rDistance < rMinDistance )
        {
          rMinDistance = rDistance;
          iMinIndex = iTarget;
        }
      }
    }
    if( iMinIndex != -1 )
    {
      rError += rMinDistance;
      mp.v2FstView = vMoveEdgelets[iMove].v2Pose;
      mp.v2SndView = vTargetEdgelets[iMinIndex].v2Pose;
      vMatchedPairs.push_back(mp);
      viIndexes.push_back(iMinIndex);
    }
    else
    {
      // if that is not the case, try to assess whether an edge exists in this neighbourhood but was incorrectly penalised by canny threshold
    /*  if (concealedEdge(vMoveEdgelets[iMove]))
      {
          concealedEdgesCount++;
        std::cerr << "Concealed Edge Detected " << std::endl;
      //  rError += rErr2/10;
      }
      else*/
        rError += rErr2;

      
    }
    fullIndexes.push_back(iMinIndex);
  }
  rError /= (vMoveEdgelets.size()-concealedEdgesCount);
  return rError;
}

bool EdgeObjectDetection::concealedEdge (Edgelet e)
{
      CVD::ImageRef edgeLocation;
      edgeLocation.x = e.v2Pose[0];
      edgeLocation.y = e.v2Pose[1];
      REAL_TYPE edgeValue = mEdgeImage[edgeLocation];
      return (edgeValue > 0);
}

REAL_TYPE EdgeObjectDetection::findClosestEdges( std::vector<MatchedPair> & vMatchedPairs, std::vector<int> & viIndexes, const std::vector<Edgelet> & vMoveEdgelets, const std::vector<Edgelet> & vTargetEdgelets, REAL_TYPE rErr1, REAL_TYPE rErr2, REAL_TYPE rAngleThreshold )
{
  REAL_TYPE rLamda = 2.f;
  REAL_TYPE rError = 0.f;
  vMatchedPairs.clear();
  vMatchedPairs.reserve( vMoveEdgelets.size() );
  viIndexes.clear();
  viIndexes.reserve( vMoveEdgelets.size() );
  for( unsigned int iMove = 0; iMove < vMoveEdgelets.size(); ++iMove )
  {
    REAL_TYPE rMinDistance = rErr1;
    int iMinIndex = -1;
    REAL_TYPE rAngle;
    MatchedPair mp;
    for( unsigned int iTarget = 0; iTarget < vTargetEdgelets.size(); ++iTarget )
    {
      const Edgelet & eTarget = vTargetEdgelets[iTarget];
      const Edgelet & eMove = vMoveEdgelets[iMove];
      rAngle = fabs( atan( (eTarget.rSlope - eMove.rSlope)/(1.f + eTarget.rSlope*eMove.rSlope) ) );
      if( rAngle < rAngleThreshold )
      {
        REAL_TYPE rDistance = norm( eTarget.v2Pose - eMove.v2Pose ) + rLamda * rAngle; 
        if( rDistance < rMinDistance )
        {
          rMinDistance = rDistance;
          iMinIndex = iTarget;
        }
      }
    }
    if( iMinIndex != -1 )
    {
      rError += rMinDistance;
      mp.v2FstView = vMoveEdgelets[iMove].v2Pose;
      mp.v2SndView = vTargetEdgelets[iMinIndex].v2Pose;
      vMatchedPairs.push_back(mp);
      viIndexes.push_back(iMinIndex);
    }
    else
    {
      rError += rErr2;
    }
  }
  rError /= vMoveEdgelets.size();
  return rError;
}

REAL_TYPE EdgeObjectDetection::iterativeClosestEdges2_NP( TooN::Matrix<3> & h, std::vector<int>& fullIndexes, std::vector<int> & viFoundEdgeIndexes, std::vector<Edgelet> & vMoveEdgelets, std::vector<Edgelet> & vTargetEdgelets, std::vector<MatchedPair> & vMatchedPairs, REAL_TYPE rErr1, REAL_TYPE rErr2, REAL_TYPE rAngleThreshold )
{
  std::vector<int> vFstIndexes, vSndIndexes;
  std::vector<int> * vptrIndexes = &vFstIndexes;
 // std::vector<int> * vpfullIndexes = &vFstIndexes;
  //std::vector<MatchedPair> vMatchedPairs;
  bool bSwitch=true;
  REAL_TYPE error = 77777.77;
  REAL_TYPE clError = findClosestEdges_NP( vMatchedPairs, fullIndexes, *vptrIndexes, vMoveEdgelets, vTargetEdgelets, rErr1, rErr2, rAngleThreshold );
  
  if( vMatchedPairs.size() < 10 )
  {
    //std::cerr << "Less than 10 matched pairs" << std::endl;
    return error;
  }
  if(  !mHomography.estimateHomography( h, vMatchedPairs, vMatchedPairs.size() ) )
  {
    //std::cerr <<"Error from Homography !" << std::endl;
    return error;
  }

  std::vector<Edgelet> vNewPoseEdgelets;
  project_edgelets( vNewPoseEdgelets, vMoveEdgelets, h, false );
  const REAL_TYPE minDiffError=1.f;
  const int max_trial = 10;
  REAL_TYPE oldError = clError + 2.f;
  TooN::Matrix<3> H;
  int no_trial = 0;
  while( oldError - clError > minDiffError && no_trial < max_trial )
  {
    oldError = clError;
    H = h;
    if( bSwitch ) vptrIndexes = &vSndIndexes; else vptrIndexes = &vFstIndexes;
    bSwitch = !bSwitch; 
    clError = findClosestEdges( vMatchedPairs, fullIndexes, *vptrIndexes, vNewPoseEdgelets, vTargetEdgelets, rErr1, rErr2, rAngleThreshold );
    TooN::Matrix<3> dH;
    if( clError >= oldError || vMatchedPairs.size() < 4 || !mHomography.estimateHomography( dH, vMatchedPairs, vMatchedPairs.size() ) )
    {
      if( bSwitch ) vptrIndexes = &vSndIndexes; else vptrIndexes = &vFstIndexes;
      clError = oldError;
      h = H;
      break;
    }
    h = dH*h;
    project_edgelets( vNewPoseEdgelets, vMoveEdgelets, h, false );
    ++no_trial;
  }
  viFoundEdgeIndexes = *vptrIndexes;
//  fullIndexes = *vpfullIndexes;
//  std::cout <<"no iterations " << no_trial << std::endl;
  return clError/vMoveEdgelets.size();
}


REAL_TYPE EdgeObjectDetection::iterativeClosestEdges2( TooN::Matrix<3> & h, std::vector<int>& fullIndexes, std::vector<int> & viFoundEdgeIndexes, std::vector<Edgelet> & vMoveEdgelets, std::vector<Edgelet> & vTargetEdgelets, std::vector<MatchedPair> & vMatchedPairs, REAL_TYPE rErr1, REAL_TYPE rErr2, REAL_TYPE rAngleThreshold )
{
  std::vector<int> vFstIndexes, vSndIndexes;
  std::vector<int> * vptrIndexes = &vFstIndexes;
 // std::vector<int> * vpfullIndexes = &vFstIndexes;
  //std::vector<MatchedPair> vMatchedPairs;
  bool bSwitch=true;
  REAL_TYPE error = 77777.77;
  REAL_TYPE clError = findClosestEdges( vMatchedPairs, fullIndexes, *vptrIndexes, vMoveEdgelets, vTargetEdgelets, rErr1, rErr2, rAngleThreshold );
  
  if( vMatchedPairs.size() < 10 )
  {
    //std::cerr << "Less than 10 matched pairs" << std::endl;
    return error;
  }
  if(  !mHomography.estimateHomography( h, vMatchedPairs, vMatchedPairs.size() ) )
  {
    //std::cerr <<"Error from Homography !" << std::endl;
    return error;
  }

  std::vector<Edgelet> vNewPoseEdgelets;
  project_edgelets( vNewPoseEdgelets, vMoveEdgelets, h, false );
  const REAL_TYPE minDiffError=1.f;
  const int max_trial = 10;
  REAL_TYPE oldError = clError + 2.f;
  TooN::Matrix<3> H;
  int no_trial = 0;
  while( oldError - clError > minDiffError && no_trial < max_trial )
  {
    oldError = clError;
    H = h;
    if( bSwitch ) vptrIndexes = &vSndIndexes; else vptrIndexes = &vFstIndexes;
    bSwitch = !bSwitch; 
    clError = findClosestEdges( vMatchedPairs, fullIndexes, *vptrIndexes, vNewPoseEdgelets, vTargetEdgelets, rErr1, rErr2, rAngleThreshold );
    TooN::Matrix<3> dH;
    if( clError >= oldError || vMatchedPairs.size() < 4 || !mHomography.estimateHomography( dH, vMatchedPairs, vMatchedPairs.size() ) )
    {
      if( bSwitch ) vptrIndexes = &vSndIndexes; else vptrIndexes = &vFstIndexes;
      clError = oldError;
      h = H;
      break;
    }
    h = dH*h;
    project_edgelets( vNewPoseEdgelets, vMoveEdgelets, h, false );
    ++no_trial;
  }
  viFoundEdgeIndexes = *vptrIndexes;
//  fullIndexes = *vpfullIndexes;
//  std::cout <<"no iterations " << no_trial << std::endl;
  return clError/vMoveEdgelets.size();
}

REAL_TYPE EdgeObjectDetection::iterativeClosestEdges( TooN::Matrix<3> & h, std::vector<int> & viFoundEdgeIndexes, std::vector<Edgelet> & vMoveEdgelets, std::vector<Edgelet> & vTargetEdgelets, REAL_TYPE rErr1, REAL_TYPE rErr2, REAL_TYPE rAngleThreshold )
{
  std::vector<int> vFstIndexes, vSndIndexes;
  std::vector<int> * vptrIndexes = &vFstIndexes;
  std::vector<MatchedPair> vMatchedPairs;
  bool bSwitch=true;
//  double angle = 0;
  REAL_TYPE error = 77777.77;
  REAL_TYPE clError = findClosestEdges( vMatchedPairs, *vptrIndexes, vMoveEdgelets, vTargetEdgelets, rErr1, rErr2, rAngleThreshold );
  if( vMatchedPairs.size() < 10 )
  {
    //std::cerr << "Less than 10 matched pairs" << std::endl;
    return error;
  }
  if(  !mHomography.estimateHomography( h, vMatchedPairs, vMatchedPairs.size() ) )
  {
    //std::cerr <<"Error from Homography !" << std::endl;
    return error;
  }

  std::vector<Edgelet> vNewPoseEdgelets;
  project_edgelets( vNewPoseEdgelets, vMoveEdgelets, h, false );
  const REAL_TYPE minDiffError=1.f;
  const int max_trial = 10;
  REAL_TYPE oldError = clError + 2.f;
  TooN::Matrix<3> H;
  int no_trial = 0;
  while( oldError - clError > minDiffError && no_trial < max_trial )
  {
    oldError = clError;
    H = h;
    if( bSwitch ) vptrIndexes = &vSndIndexes; else vptrIndexes = &vFstIndexes;
    bSwitch = !bSwitch; 
    clError = findClosestEdges( vMatchedPairs, *vptrIndexes, vNewPoseEdgelets, vTargetEdgelets, rErr1, rErr2, rAngleThreshold );
    TooN::Matrix<3> dH;
    if( clError >= oldError || vMatchedPairs.size() < 4 || !mHomography.estimateHomography( dH, vMatchedPairs, vMatchedPairs.size() ) )
    {
      if( bSwitch ) vptrIndexes = &vSndIndexes; else vptrIndexes = &vFstIndexes;
      clError = oldError;
      h = H;
      break;
    }
    h = dH*h;
    project_edgelets( vNewPoseEdgelets, vMoveEdgelets, h, false );
    ++no_trial;
  }
 // angle = calculateRotationAngle (vMatchedPairs);
  viFoundEdgeIndexes = *vptrIndexes;
//  std::cout <<"no iterations " << no_trial << std::endl;
  return clError/vMoveEdgelets.size();
}

double EdgeObjectDetection::get_memory()
{
   char file_str[4096], dum_str[4096];
   int file_ptr = -1, file_len;
   file_ptr = open ("/proc/self/stat", O_RDONLY);
   if (file_ptr < 0)
     std::cout << "COULD NOT OPEN /proc/self/stat" << std::endl;
   else
     std::cout << "file pointer = " << file_ptr << std::endl;
   file_len = read (file_ptr, file_str, sizeof(file_str)-1);
  // close (file_ptr);
   file_str[file_len] = '\0';
   int dum_int;
   unsigned int dum_uint, vm_size, rss;
   static int page_size = getpagesize();
   int num_fields = sscanf(file_str,
                           "%d " // pid
                           "%s " // comm
                           "%c " // state
                           "%d %d %d %d %d " //ppid, pgrp, session, tty, tpgid
                           "%u %u %u %u %u " // flags, minflt, cminflt, majflt, cmajflt
                           "%d %d %d %d %d %d " // utime, stime, ctime, cstime, counter, priority
                           "%u %u " // timeout, itrealvalue
                           "%d " // starttime
                           "%u %u", //vsize, rss
                           &dum_int, dum_str, dum_str, &dum_int, &dum_int, &dum_int, &dum_int, &dum_int,
                           &dum_uint, &dum_uint, &dum_uint, &dum_uint, &dum_uint,
                           &dum_int, &dum_int, &dum_int, &dum_int, &dum_int, &dum_int,
                           &dum_uint, &dum_uint,
                           &dum_int,
                           &vm_size, &rss);
    double mem_size = (vm_size/page_size);
    return mem_size;
}

