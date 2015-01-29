#include "MODetector.h"

#include <gvars3/instances.h>
#include <cvd/image_io.h>
#include "GetTransformDist.h"
#include "BitMacros.h"
#include "CannyEdgeDetector.h"
#include "EdgeObjectDetection.h"
#include <cvd/vector_image_ref.h>
#include <QTimer>
#include <ros/ros.h>
#include "CameraSource.h"
#include <QtGui/QApplication>

#include "MainWindow.h"
#include <cvd/timer.h>
#include <stdio.h>

using namespace GVars3;
using namespace std;


namespace mod_base
{

MODetector::MODetector ()
{
   run_yet = false;
}

void MODetector::Start (ros::Publisher detectionImage_pub_, ros::Publisher detectedO_, std::string colorImgSrc, std::string cf, std::string lf, std::string of, int imgSize, int learnObjNo, bool wI, bool iS, bool aP, bool aF, int scNo)
{ 
   if (!run_yet)
   {
       codebookFilename = cf;
       loadCodebookFilename = lf;
       withInterface = wI;
       andPublish = aP;
       isStar = iS;
       andWriteToFile = aF;
       outputFilename = of;
       scenarioNo = scNo;
       detected_objs_ = detectedO_;
       run(detectionImage_pub_, colorImgSrc, imgSize, learnObjNo);
       run_yet = true;
   }
}

void MODetector::run(ros::Publisher detectionImage_pub_, std::string colorImgSrc, int imgSize, int learnObjNo)
{
   runDetector(detectionImage_pub_, colorImgSrc, imgSize, learnObjNo, 0, NULL);
}

void MODetector::runDetector (ros::Publisher detectionImage_pub_, std::string colorImgSrc, int imgSize, int learnObjNo, int argc, char** argv)
{
  int img_h = 320, img_w = 240;
  int window_x = 220, window_y = 300;
  char codebookName [100];
  char gtFileName [100];
  sprintf (gtFileName, "/home/cognito2/Dataset_Bristol2010/TestImages_undistorted/yellow_gt.txt");
  sprintf(codebookName, "seqCB.eod.bin");
  char folderName[100];
  char outputFolderName[100];
  sprintf(folderName, "../img-ccg4/");
  char dbFName[100];
  sprintf (dbFName, "Dataset_Bristol2010");
  char bbFilename[100];
  char outputFileName[100];
  sprintf(outputFileName, "data-log2.txt");
  sprintf(outputFolderName, "/home/staff/damen/ubuntu/ros/cognito/cognito_vision/cognito_mod/f3");
  char imgtype [10];
  int imgSrc;
  sprintf(imgtype, "png");
  int lastFrame = 100;
  int firstFrame = 1;
  double maxTime = 0.1;
  bool lsd = false;
  const int maxObjNo = 20;
  int includeObjects [maxObjNo];
  for (int i = 0; i < 20; i++)
     includeObjects[i] = 0;

  bool createCodebook = false;

  LLImagePyramid<CVD::byte> mLLImagePyramid( CVD::ImageRef(img_h,img_w), 1 );
  LLImagePyramid<CVD::byte> mLLImagePyramidK( CVD::ImageRef(img_h*2, img_w*2), 2);
  if (!createCodebook){
    ros::init(argc, argv, "live_video_rec");
    CameraSource * pRealCamera ( new CameraSource );
    std::cout << imgSize;
    if (imgSize == 1)
      pRealCamera->setImgSize (CVD::ImageRef (1280, 960));
    else if (imgSize == 2)
      pRealCamera->setImgSize (CVD::ImageRef (640, 480));
    pRealCamera->setSources (colorImgSrc);
    pRealCamera->connect ( "firewire" );

    QApplication a(argc, argv);
    MainWindow w (pRealCamera, outputFolderName, outputFileName, maxTime, gtFileName, window_x, window_y, 0, 0, 1, isStar, lsd);
   // w.setKinectImgSize (imgSize);
    w.setPublishers (detectionImage_pub_);
    w.setCodebookFilename (codebookFilename);
    if (loadCodebookFilename != "")
    {
       w.loadCodebookFromFile (loadCodebookFilename);
    }
    else
    {
        // here we need to create an empty codebook and enable multiple paths to be added to it!
        w.initialiseCodebook();
        std::cerr << "Initialised Empty Codebook " << std::endl;
    }
    w.setUndistort ( false );
    w.setLearn(0, learnObjNo);
    w.setMask (false);
    w.setPublishing (andPublish);
    if (andWriteToFile)
    {
		w.setPrintToFile (outputFilename);
	}
    w.setObjsPublisher (detected_objs_);
    w.setScenarioNo (scenarioNo);
    if (withInterface)
    {
       w.show ();
    }
    bool result = a.exec();
    delete pRealCamera;

  }
}
}
