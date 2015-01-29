/*
 * =====================================================================================
 *
 *       Filename:  EdgeObjectDetection.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  27/07/10 01:35:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Pished Bunnun (pbunnun), pbunnun@cs.bris.ac.uk
 *        Company:  Computer Science Department, University of Bristol
 *
 * =====================================================================================
 */

#ifndef EDGEOBJECTDETECTION_H
#define EDGEOBJECTDETECTION_H

#include <cvd/image.h>
#include <cvd/image_ref.h>
#include <cvd/image_io.h>
#include <TooN/TooN.h>
#include "DefineTypes.h"
#include "DataStructure4Detection.h"
#include "Homography.h"
#include "LowLevelImageData.h"
#include "CannyEdgeDetector.h"
#include "lsd.h"

class EdgeObjectDetection
{
  public:
    EdgeObjectDetection( CVD::ImageRef irSize, char* outputFolderName, char* outputFileName, bool isStar, bool lsd, double maxTime, int imgSizeRatio );
    ~EdgeObjectDetection();
    /**
      Pre calculate EdgeLnks from edge points specified by edges.
      @param[out] binOfEdgeLinks
      @param[in] edges contains edge points generated from Canny edge detector.
     */
    void preCalculateEdgeLnks( LowLevelImageData<CVD::byte> & llImage, int master_file_id);
    void preCalculateEdgeLnksWithRegions( LowLevelImageData<CVD::byte> & llImage, std::vector<SegmentedRegion> imageRegions, std::vector<SegmentedMask> imageMasks, int master_file_id);
    void preCalculateEdgeLnksWithRegions_MultiScale( LowLevelImageData<CVD::byte> & llImage, std::vector<SegmentedRegion> imageRegions, std::vector<SegmentedMask> imageMasks, int master_file_id);

    void preCalculateEdgeLnks_depth( LowLevelImageData<CVD::byte> & llImage, CVD::Image<CVD::byte> & depth_image, int master_file_id);
    void preCalculateEdgeLnksWithRegions_depth( LowLevelImageData<CVD::byte> & llImage, CVD::Image<CVD::byte> & depth_image, std::vector<SegmentedRegion> imageRegions, int master_file_id);
    void getDistanceTransformImage (CVD::Image<CVD::Rgb<CVD::byte> >& transformImage);
    void getEdgeImage (CVD::Image<CVD::Rgb<CVD::byte> >& transformImage);

    /**
      Detect an object.
      @param[out] H is a homography.
      @param[in] iObjectID is a specific id of a object to be detected.
      @param[in] binOfEdgeLinks is generated from preCalculateEdgeLnks().
      @return true if it find the object. 
     */
    int detect(int a);
    int detect_star(int a);
    int detect_constellation(int a);
    bool detect_from_previous (int a, int prev_obj_no);
    void setDepth (bool d) {isDepth = d;}
    void setMask (bool m) {isMask = m;}
//    inline const TooN::SE3<>& get_se3W2C( unsigned int iObjectID, int no_se3 = 0 ) { assert(iObjectID < mvObjectDescriptors.size()); return mvObjectDescriptors[iObjectID].vse3W2C[no_se3]; }
//    inline const TooN::SE3<>& get_se3W2COriginal( unsigned int iObjectID ) { assert(iObjectID < mvObjectDescriptors.size()); return mvObjectDescriptors[iObjectID].se3W2COriginal; }
//    inline const TooN::Matrix<3>& get_homography( unsigned int iObjectID ) { assert(iObjectID < mvObjectDescriptors.size()); return mvObjectDescriptors[iObjectID].m3H; }
//    inline const std::vector<TooN::Vector<3> >& get_visible_points( unsigned int iObjectID ) { assert(iObjectID < mvObjectDescriptors.size()); return mvObjectDescriptors[iObjectID].vv3VisiblePoints; }
//    inline const std::vector<TooN::Vector<3> >& get_template_points( unsigned int iObjectID ) { assert(iObjectID < mvObjectDescriptors.size()); return mvObjectDescriptors[iObjectID].vv3Templates; }
    void set_camera_parameters( TooN::Matrix<3> & K ) { mHomography.set_camera_parameters( K ); }

    std::vector<TooN::Vector<2> > mv2DetectedEdges;
    static const REAL_TYPE pi = 3.14;
    static const REAL_TYPE multiplyCell = 0.1;
    EdgeCellMatrix *ecm;
    CellDirectStruct* cds;

    /**
      Generate line segments from edge points.
      @param[out] vLines contains line segments.
      @param[in] edges contains edge points generated from Canny edge detector.
     */
    void get_lines(std::vector<Line> & vLines, const Edges & edges);
    void get_lines_DT(std::vector<Line> & vLines, const Edges & edges, const CVD::Image<REAL_TYPE> & mvEdgeDistanceTransformImage);
    inline int get_no_object( ) { return miTotalObjectViews; };
    void gen_Edgelets( std::vector<Edgelet> & vEdgelets, std::vector<int> & viLineIndexes, CVD::ImageRef irTopLeft, CVD::ImageRef irButtomRight, const std::vector<Line> & vLines, int iMaxLength = 5 );
    void gen_Edgelets( std::vector<Edgelet> & vEdgelets, std::vector<int> & viLineIndexes, const std::vector<Line> & vLines, int iMaxLength = 5 );
    void gen_Edgelets_with_Directional( std::vector<Edgelet> & vEdgelets, std::vector<int> & viLineIndexes, const std::vector<Line> & vLines, long oriSize = 76000, int iMaxLength = 5);
    void gen_Edgelets_with_Directional_Regions( std::vector<Edgelet> & vEdgelets, std::vector<int> & viLineIndexes, const std::vector<Line> & vLines, const std::vector<SegmentedRegion> regions, int maxValue, int iMaxLength = 5, long oriSize = 76000);
    void gen_Edgelets_with_Directional_Masks( std::vector<Edgelet> & vEdgelets, std::vector<int> & viLineIndexes, const std::vector<Line> & vLines, const std::vector<SegmentedMask> masks, long oriSize = 76000, int iMaxLength = 5);
    void gen_TestDirections( const std::vector<Edgelet> & vEdgelets, const std::vector<int>& viLineIndexes );
    void gen_TestDirections( const std::vector<Edgelet> & vEdgelets ); 
    void gen_EdgeletCells ( const std::vector<Edgelet> & vEdgelets );
    void gen_CertainDirectionsChain( std::vector<CodeBookElement> & mvCodeBook, const std::vector<Edgelet> & vEdgelets, const TooN::Vector<4> & v4ChainAngles, int view_ID, int object_ID, int obj_correct_ID);
    void gen_CertainDirectionsChain_star( std::vector<CodeBookElement> & mvCodeBook, const std::vector<Edgelet> & vEdgelets, const TooN::Vector<4> & v4ChainAngles, int view_ID, int object_ID, int obj_correct_ID);
    void getChainInDirection_withCodebook (std::vector<int> &returnedEdgelets, std::vector<std::vector<CodeBookElement>::iterator> &returnedCB, const std::vector<Edgelet> testEdgelets, std::vector<TestDirection> Q, int edgeNo, TooN::Vector<4> mvTestDirections, unsigned int chainNo, EdgeCellMatrix* edgeletImageMatrix, REAL_TYPE prevDir, REAL_TYPE prevDist, std::vector<std::vector<CodeBookElement>::iterator> mvCBIterator, REAL_TYPE minDir, REAL_TYPE maxDir);
    void getChain_level1 (std::vector<int> &returnedEdgelets, std::vector<std::vector<CodeBookElement>::iterator> &returnedCB, const std::vector<Edgelet> testEdgelets, std::vector<TestDirection> Q, int edgeNo, TooN::Vector<4> mvTestDirections, EdgeCellMatrix* edgeletImageMatrix, REAL_TYPE prevDir, REAL_TYPE prevDist, std::vector<std::vector<CodeBookElement>::iterator> mvCBIterator, REAL_TYPE minDir, REAL_TYPE maxDir);
    void getChain_level2 (std::vector<int> &returnedEdgelets, std::vector<std::vector<CodeBookElement>::iterator> &returnedCB, const std::vector<Edgelet> testEdgelets, std::vector<TestDirection> Q, int edgeNo, TooN::Vector<4> mvTestDirections, EdgeCellMatrix* edgeletImageMatrix, REAL_TYPE prevDir, REAL_TYPE prevDist, std::vector<std::vector<CodeBookElement>::iterator> mvCBIterator, REAL_TYPE minDir, REAL_TYPE maxDir);
   void getChain_level3 (std::vector<int> &returnedEdgelets, std::vector<std::vector<CodeBookElement>::iterator> &returnedCB, const std::vector<Edgelet> testEdgelets, std::vector<TestDirection> Q, int edgeNo, TooN::Vector<4> mvTestDirections, EdgeCellMatrix* edgeletImageMatrix, REAL_TYPE prevDir, REAL_TYPE prevDist, std::vector<std::vector<CodeBookElement>::iterator> mvCBIterator, REAL_TYPE minDir, REAL_TYPE maxDir);


    EdgeCellMatrix* getEdgeletImageMatrix (const std::vector<Edgelet> & vEdgelets, int h, int w, int cellH, int cellW);


    int addCodeBook(const std::vector<Line> & lines, const Edges edges, CVD::ImageRef irTopLeft, CVD::ImageRef irButtomRight, const TooN::Vector<4> & v4ChainAngles, int object_ID, int obj_correct_ID); // return the number of Objects

    void load_cellDirectStructure (char* filename);
    std::vector<DetectedObject> & get_detected_objects( ) { return mvDetectedObjects; };
    std::vector<DetectedObject> & get_previous_detected_objects() { return mvDetectedObjects_previousFrame; };
    std::vector<DetectedObject> & get_candidate_objects( ) { return mvCandidateObjects; };
    std::vector<CheckedChains> & get_checked_chains() { return mvCheckedChains; };
    int getLatestObjectID () { return latestObjID; };
    void setLatestObjectID (int i);
    int getLatestCorrectObjectID () { return latest_correctObjID; };
    void getFoundChains (std::vector<CheckedChains> &);
    void getMatchedConstellations (std::vector< TooN::Vector<2> > & thisChecked);
    inline bool isCodeBookEmpty() { return mbEmpty; };
    double get_memory ();

    void initialise_codebook();
    void load_codebook(std::string filename);
    void load_codebook(std::string filename, int class_no);
    void clear_codebooks();

    void save_codebook(std::string filename);
    // Helping function
    bool addView();
    bool addView(int objNo, SegmentedRegion region, bool AndSave);
    bool addView(int objNo, SegmentedRegion region, SegmentedMask mask, bool AndSave);
    bool isViewAdded () { return viewAdded; };
    int getNumberOfLearntObjects () { return learntObjNo; }
    int getNumberOfLearntViews () { return learntViewNo; }
    double getNumberOfLearntConstellations() { return mvObjDescriptorClasses[0].mvCodeBook.size(); }
    int max (int a, int b);
    int min (int a, int b);
    char pointcloudFilename [200];
    char viewFilename [200];
    char maskFilename [200];
    char objectFolderName [200];
    char bbFileName [200];
    char rangeEdgeFilename[200];
    int addCodeBook( LowLevelImageData<CVD::byte> & llImage, CVD::ImageRef irTopLeft, CVD::ImageRef irButtomRight, const TooN::Vector<4> & v4ChainAngles, int object_ID, int obj_correct_ID);
    void sortCodeBooks();

    void generateOfflineCellDirectStruct (int h, int w, int cellH, int cellW, REAL_TYPE cellDirInterval);
    std::vector<Edgelet> findAllPixelsWithin(int h, int w, int y, int x, REAL_TYPE dir, const int first_dist);
    void getLSDEdges (Edges &);
    void getEdges (Edges &);
    void getLSDEdgesInRegions (Edges & e, std::vector<SegmentedRegion> imR);
    int getCellNo (Edgelet e);

    CVD::ImageRef getViewTopLeft () {return thisIrTopLeft; };
    CVD::ImageRef getViewBottomRight () {return thisIrBottomRight; };
    CVD::Image<CVD::byte> getEdgeImage() {return mEdgeImage; };
    std::vector<Edgelet> getEdgelets () {return mvEdgelets;};

    void set_maxTime (REAL_TYPE new_max) {MAX_TIME = new_max;};


    int resetLatestObj ();

  private:
    bool calculateTwoEdgelets( TestDirection & testDirection, const std::vector<Edgelet> & vEdgelets, int iIdx_0, int iIdx_1 );
    int calculateTwoCellsDirection (int cell1, int cell2, double cellInterval);

    bool concealedEdge (Edgelet e);

    bool calculateEdgeletAndCell( CellDirection & cd, const std::vector<Edgelet> & vEdgelets, int edgeNo, int iIdx_1 );
    inline REAL_TYPE calculateAngle( REAL_TYPE theta2, REAL_TYPE theta1 )
    {
#if 0    
      REAL_TYPE dX1 = cos(theta1);
      REAL_TYPE dY1 = sin(theta1);
      REAL_TYPE dX2 = cos(theta2);
      REAL_TYPE dY2 = sin(theta2);
      REAL_TYPE result = acos( dX1*dX2 + dY1*dY2 );
      if( asin(dX1*dY2 - dX2*dY1) < 0 )
        result *= -1.f;
      return result;
#else      
      REAL_TYPE theta = theta2 - theta1;
      return ( sin(theta) < 0 ) ? -acos(cos(theta)) : acos(cos(theta)); // fource theta to be [-pi,pi];
#endif      
    }
    ;
    void lsd_to_lines( std::vector<Line> & vLines, ntuple_list lsd_format );
    void lines_to_edges(Edges & points, std::vector<Line> vLines); 
    void project_edgelets( std::vector<Edgelet> & vProjectedEdgelets, const std::vector<Edgelet> & vEdgelets, const TooN::Matrix<3> & h, bool bWrite=false );
    void project_edges( Edges & vProjectedEdgelets, std::vector<TooN::Vector<2> >& vv2Pnts, const Edges & vEdgelets, const TooN::Matrix<3> & h, bool bWrite=false );
    REAL_TYPE findClosestEdges( std::vector<MatchedPair> & vMatchedPairs, std::vector<int> & viIndexes, const std::vector<Edgelet> & vMoveEdgelets, const std::vector<Edgelet> & vTargetEdgelets, REAL_TYPE rErr1 = 15.f, REAL_TYPE rErr2 = 30.f, REAL_TYPE rAngleThreshold = 0.3 );
    REAL_TYPE findClosestEdges( std::vector<MatchedPair> & vMatchedPairs, std::vector<int> & fullIndexes, std::vector<int> & viIndexes, const std::vector<Edgelet> & vMoveEdgelets, const std::vector<Edgelet> & vTargetEdgelets, REAL_TYPE rErr1 = 15.f, REAL_TYPE rErr2 = 30.f, REAL_TYPE rAngleThreshold = 0.3 );
    REAL_TYPE findClosestEdges_NP( std::vector<MatchedPair> & vMatchedPairs, std::vector<int> & fullIndexes, std::vector<int> & viIndexes, const std::vector<Edgelet> & vMoveEdgelets, const std::vector<Edgelet> & vTargetEdgelets, REAL_TYPE rErr1 = 15.f, REAL_TYPE rErr2 = 30.f, REAL_TYPE rAngleThreshold = 0.3 );
    REAL_TYPE iterativeClosestEdges( TooN::Matrix<3> & h, std::vector<int> & viFoundEdgeIndexes, std::vector<Edgelet> & vMoveEdgelets, std::vector<Edgelet> & vTargetEdgelets, REAL_TYPE rErr1 = 15.f, REAL_TYPE rErr2 = 30.f, REAL_TYPE rAngleThreshold = 0.3 );
    REAL_TYPE iterativeClosestEdges2( TooN::Matrix<3> & h, std::vector<int> & fullIndexes, std::vector<int> & viFoundEdgeIndexes, std::vector<Edgelet> & vMoveEdgelets, std::vector<Edgelet> & vTargetEdgelets, std::vector<MatchedPair> & vMatchedPairs, REAL_TYPE rErr1 = 15.f, REAL_TYPE rErr2 = 30.f, REAL_TYPE rAngleThreshold = 0.3 );
    REAL_TYPE iterativeClosestEdges2_NP( TooN::Matrix<3> & h, std::vector<int> & fullIndexes, std::vector<int> & viFoundEdgeIndexes, std::vector<Edgelet> & vMoveEdgelets, std::vector<Edgelet> & vTargetEdgelets, std::vector<MatchedPair> & vMatchedPairs, REAL_TYPE rErr1 = 15.f, REAL_TYPE rErr2 = 30.f, REAL_TYPE rAngleThreshold = 0.3 );


    // added to speed up the search within cells
   std::vector<unsigned int> getEdgeletsInDirection (const std::vector<Edgelet> testEdgelets, int edgeNo, EdgeCellMatrix* edgeletImageMatrix, int cellH, int cellW, REAL_TYPE dir);

    /**
     Return a side of point (v2Trd) relatived to a line ( v2Fst->v2Snd ) as a boolean value.
     @param[in] v2Check = v2Trd - v2Fst.
     @param[in] v2Line = v2Snd - v2Fst.
     */
    inline bool checkSide( TooN::Vector<2> & v2Check, TooN::Vector<2> & v2Line )
    {
      return ( v2Line[0]*v2Check[1]-v2Line[1]*v2Check[0] > 0 );
    };

    int miTotalObjectViews;

    bool mbEmpty;

    Homography mHomography;

    std::vector<int> returnedEdgelets;
    std::vector<std::vector<CodeBookElement>::iterator> returnedCB;
    std::vector<std::vector<CodeBookElement>::iterator> thisCB;
    
    CannyEdgeDetector mCannyEdgeDetector;
 public:
    std::vector<ObjDescriptorClass> mvObjDescriptorClasses;
    std::vector<ObjDescriptorClass> mvIndObjDescriptorClasses;

 private:
    std::vector<CVD::Image<CVD::byte> > mvEdgeOrientationImages;
    std::vector<CVD::Image<REAL_TYPE> > mvEdgeOrientationDistTransformImages;
    std::vector<Orientation> mvEdgeOrientations;
    std::vector<Orientation> mvEdgeOrientationDistTransforms;
    CVD::Image<REAL_TYPE> mvEdgeDistanceTransformImage;
    CVD::Image<REAL_TYPE> mvEdgeDistanceTransformImageIntegral;

    // share parameters
    CVD::Image<CVD::byte> mEdgeImage;
    CVD::ImageRef mirImageSize;
    std::vector<Edgelet> mvEdgelets;
    Edges edges;
    std::vector<int> mviLineIndexes;
    std::list<TestDirection> mlTestDirections;
    std::vector<TestDirection> mvTestDirections;
    std::vector<std::vector<TestDirection> > mvTestDirectionsPreviousObjects;
    std::vector<CellDirection> mvCellDirections;
    std::vector<TestDirection> vCellTestDirections;
    std::vector<CVD::ImageRef> mvirDirectionIndexes;
    std::vector<DetectedObject> mvDetectedObjects;
    std::vector<DetectedObject> mvDetectedObjects_previousFrame;
    std::vector<DetectedObject> mvDetectedObjects_undetected;
    std::vector<DetectedObject> mvCandidateObjects;
    std::vector<CheckedChains> mvCheckedChains;

    // for keeping the latest detected object ID (and correct object ID), in case a view is to be added
    int latestObjID;
    int latest_correctObjID;
    LowLevelImageData<CVD::byte> thisLlImage;
    CVD::ImageRef thisIrTopLeft;
    CVD::ImageRef thisIrBottomRight;
    bool viewAdded;
    int imgSizeRatio;
    bool isDepth;
    bool isMask;
    bool isStar;
    int learntObjNo, learntViewNo, learntConstNo;
    int lineMaxLen; 
    REAL_TYPE MAX_TIME;
};

#endif
