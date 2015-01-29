/*
 * =====================================================================================
 *
 *       Filename:  EdgeExtraction.h
 *
 * =====================================================================================
 */

#ifndef EDGEEXTRACTION_H
#define EDGEEXTRACTION_H

#include <cvd/image.h>
#include <cvd/image_ref.h>
#include <TooN/TooN.h>
#include "DefineTypes.h"
#include "DataStructure4Detection.h"
#include "Homography.h"
#include "LowLevelImageData.h"
#include "CannyEdgeDetector.h"
#include "lsd.h"
#include <sensor_msgs/Image.h>
#include "LLImagePyramid.h"
#include <sensor_msgs/image_encodings.h>
#include <cv_bridge/cv_bridge.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

class EdgeExtraction
{
  public:
    EdgeExtraction(CVD::ImageRef irSize, bool lsd=true);
    ~EdgeExtraction();
    /**
      Pre calculate EdgeLnks from edge points specified by edges.
      @param[out] binOfEdgeLinks
      @param[in] edges contains edge points generated from Canny edge detector.
     */
    void preCalculateEdgeLnks( LowLevelImageData<CVD::byte> & llImage);
    /**
      Generate line segments from edge points.
      @param[out] vLines contains line segments.
      @param[in] edges contains edge points generated from Canny edge detector.
     */
    void get_lines(std::vector<Line> & vLines, const Edges & edges);
    void gen_Edgelets_with_Directional( std::vector<Edgelet> & vEdgelets, std::vector<int> & viLineIndexes, const std::vector<Line> & vLines, long oriSize = 76000, int iMaxLength = 5);
    void genDistanceTransform (const std::vector<Line> & vLines, long oriSize=76800);
    void getLSDEdges (Edges &);
    void load_codebook(std::string filename);
    View getView (int chainNo, int classNo, int viewNo) const;

    double getConfidenceScore (View v, TooN::Matrix<3> homography);

    /**
     * Special version of getConfidenceScore() for use by external
     * libraries, e.g. cognito_tracking
     */
    double getConfidenceScoreModified (View v, TooN::Matrix<3> homography, std::vector<Edgelet>& projectedEdgelets);

    LLImagePyramid<CVD::byte> getImageFromKinect (const sensor_msgs::ImageConstPtr& msg);

    void saveDistanceTransformImages(const std::string& directory);

  private:
    void lsd_to_lines( std::vector<Line> & vLines, ntuple_list lsd_format );
    void lines_to_edges(Edges & points, std::vector<Line> vLines); 
    void project_edgelets( std::vector<Edgelet> & vProjectedEdgelets, const std::vector<Edgelet> & vEdgelets, const TooN::Matrix<3> & h );
    void project_edges( Edges & vProjectedEdgelets, std::vector<TooN::Vector<2> >& vv2Pnts, const Edges & vEdgelets, const TooN::Matrix<3> & h );

    CannyEdgeDetector mCannyEdgeDetector;

    std::vector<CVD::Image<CVD::byte> > mvEdgeOrientationImages;
    std::vector<CVD::Image<REAL_TYPE> > mvEdgeOrientationDistTransformImages;
    std::vector<Orientation> mvEdgeOrientations;
    std::vector<Orientation> mvEdgeOrientationDistTransforms;

    // share parameters
    CVD::Image<CVD::byte> mEdgeImage;
    CVD::ImageRef mirImageSize;
    std::vector<Edgelet> mvEdgelets;
    Edges edges;
    std::vector<int> mviLineIndexes;

    std::vector<ObjDescriptorClass> mvObjDescriptorClasses;
    LLImagePyramid<CVD::byte> mLLImagePyramid;
    cv_bridge::CvImageConstPtr rgb_image;
    CVD::Image<CVD::byte> mGrayImageF;
    CVD::Image<CVD::byte> mGrayImage;

};

#endif
