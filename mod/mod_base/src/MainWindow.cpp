/********************************************************************************/
/**  Computer Vision Group, Department of Computer Science. University of      **/
/**  Bristol. This code can not be used or copied without express permission.  **/ 
/********************************************************************************/
#include "MainWindow.h"
#include "ui_mainwindow.h"
#include <QTimer>
#include <QKeyEvent>
#include <QFileDialog>
//#include "LinearMath/btTransform.h"
#include <TooN/TooN.h>
#include <gvars3/instances.h>
#include <cvd/image_convert.h>
#include <tag/threepointpose.h>
#include <cv_bridge/cv_bridge.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <geometry_msgs/PoseStamped.h>
#include <image_geometry/pinhole_camera_model.h>
#include <mod_base/DetectedObject.h>
#include <mod_base/DetectedObjectArray.h>

#include <EdgeExtraction.h>
#include <fcntl.h>

// Share image data
CVD::Image<REAL_TYPE> imRealShareImage;
std::vector<gtItem> allGT;
LLImagePyramid<CVD::byte> mLLImagePyramidForDetection(CVD::ImageRef(640,480),2);
//LLImagePyramid<CVD::byte> mLLImagePyramidForDetection(CVD::ImageRef(1280,960),3);

double thisMaxTime;
int counter = 0;
bool fromPrevious = false;
bool andAlign = false;
int pubCounter = 0;
double detectedTimes [1000];
int foundValues = 0, false_positive = 0, false_negative = 0;
int frameNo = 0;
CVD::ImageRef beginW (300,530);

MainWindow::MainWindow( CameraSource * pCameraSource, char* outputFolderName, char* outputFileName, double maxTime, char* gtFileName, int w_x, int w_y, bool k, bool g, bool fw, bool iS, bool lsd, QWidget *parent ) 
:
  QMainWindow( parent ), ui( new Ui::MainWindow ), mpCameraSource( pCameraSource ), mRgbImage( mpCameraSource->getSize() ), mRgbImageF( mpCameraSource->getSize() ),
    mRgbImageforSaving( mpCameraSource->getSize() ), mGrayImage (mpCameraSource->getSize() ), mGrayImageE (mpCameraSource->getSize() ), mGrayImageF (mpCameraSource->getSize() ), mGrayImageDummy( mpCameraSource->getSize() ), mEdgeTransformImage( mpCameraSource->getSize() ), mLLImagePyramid( pCameraSource->getSize(), 1), mLLImagePyramidForDetection_TemplateImage( pCameraSource->getSize(), 2 ), miNoFrame( 0 ), mbUndistortImage( true ), mbSaveImages( false ), mbShowMessages( true ), mbShowEdges( false ), mbShowLSDEdges( false ), mbShowChains( true), mbShowAllChains( false), mbShowLines( false), mbRunEdgeTracker( true ), frameNo2 (0),
        mbContinueImageFrame( true ), mbNextImageFrame( true ), mpSaveImageThread( NULL ), mbDetection(true), mbSnapToEdges( false), mbLearningObject( false), miDetectedObjectNo(-1), learnAndSave(false), isKinect ( false )
{
  ui->setupUi( this );
 // Set up timer
  mpTimer = new QTimer( this );
  mbPublishers = false;
  isDepth = false;
  asService = false;
  isStar = iS;
 // codebookFileName = "myExtraCodeBook.eod.bin";
  const CVD::ImageRef imgScreen (640, 480);
  imgSize = mpCameraSource->getSize();
  imgSizeRatio = imgSize.x/imgScreen.x * 2;
  regionWidth = 15*imgSizeRatio/2;
  
  //const CVD::ImageRef imgScreen (1280, 960);
  mGrayImage.resize(imgScreen);
  mRgbImage.resize (imgScreen);
  mEdgeTransformImage.resize(imgScreen);
  //beginW.x = w_x;
  //beginW.y = w_y;
  beginW.x = 0;
  beginW.y = 0;
  connect( mpTimer, SIGNAL( timeout() ), this, SLOT( Update() ) );
  mpTimer->start();
  mTime.start();
  thisMaxTime = maxTime;
  allGT.reserve (100);
  kinect = k;
  gaze = g;
  firewire = fw;
  std::cout << "kinect inside set to " << kinect << std::endl;
  mpEdgeObjectDetection = new EdgeObjectDetection(CVD::ImageRef(320,240), outputFolderName, outputFileName, isStar, lsd, maxTime, imgSizeRatio);
  mpEdgeExtraction = new EdgeExtraction(CVD::ImageRef(320,240));
}

MainWindow::~MainWindow()
{
  if (mpEdgeObjectDetection->isViewAdded())
  {
     mpEdgeObjectDetection->save_codebook(codebookFilename);
     std::cerr << "saved the codebook correctly" << std::endl;
     std::cerr << "you can now end the nodelet manager" << std::endl;
  }
  output.close();
  delete mpTimer;
  delete ui;
  if( mpSaveImageThread != NULL ) delete mpSaveImageThread;
  delete mpEdgeObjectDetection;
}


void MainWindow::setPrintToFile (std::string filename)
{
	output.open (filename.c_str());
	if (output.fail())
	{
		std::cerr << "ERROR:: FILE " << filename << " could not be opened";
		return;
	}
	andPrintToFile = true;
}

void MainWindow::addViewLearning(int objNo, SegmentedRegion region, SegmentedMask mask)
{
  if (mpEdgeObjectDetection->addView (objNo, region, mask, learnAndSave))
  {
    if (learnAndSave)
    {
      const CVD::ImageRef smallSize2 (640, 480);
      const CVD::ImageRef smallSize (320,240);
      CVD::Image<CVD::byte> mGrayImageS2;
      if (imgSizeRatio == 2)
      {
         mGrayImageS2.resize(smallSize);
         CVD::halfSample(mGrayImage, mGrayImageS2);
      }
      else
      {
         CVD::Image<CVD::byte> mGrayImageS;
         mGrayImageS.resize(smallSize2);
         mGrayImageS2.resize(smallSize);
         CVD::halfSample(mGrayImage, mGrayImageS);
         CVD::halfSample(mGrayImageS, mGrayImageS2);
      }
      CVD::img_save (mGrayImageS2, mpEdgeObjectDetection->viewFilename);
      // saving the mask as an image
      cv_bridge::CvImageConstPtr rgb_image;
      rgb_image = cv_bridge::toCvCopy(mask.image, sensor_msgs::image_encodings::MONO8);
      CVD::BasicImage<CVD::byte> cvd_rgb_image_b (rgb_image->image.data, CVD::ImageRef(rgb_image->image.cols, rgb_image->image.rows));
      CVD::Image<CVD::byte> maskSmallImage;
      maskSmallImage.resize(smallSize);
      CVD::halfSample(cvd_rgb_image_b, maskSmallImage);
      CVD::img_save (maskSmallImage, mpEdgeObjectDetection->maskFilename);
      //mpCameraSource->savePointCloud (mpEdgeObjectDetection->pointcloudFilename, xyz1.x, xyz2.x, xyz1.y, xyz2.y);
    }
  }
  
  mbLearningObject = false;
}

void MainWindow::addViewLearning(int objNo, SegmentedRegion region)
{
  std::cerr << "entering addView with object No "<< objNo << std::endl;
  if (mpEdgeObjectDetection->addView (objNo, region, learnAndSave))
  {
    if (learnAndSave)
    {
      const CVD::ImageRef smallSize2 (640, 480);
      const CVD::ImageRef smallSize (320,240);
      CVD::Image<CVD::byte> mGrayImageS2;
      if (imgSizeRatio == 2)
      {
         mGrayImageS2.resize(smallSize);
         CVD::halfSample(mGrayImage, mGrayImageS2);
      }
      else
      {
         CVD::Image<CVD::byte> mGrayImageS;
         mGrayImageS.resize(smallSize2);
         mGrayImageS2.resize(smallSize);
         CVD::halfSample(mGrayImage, mGrayImageS);
         CVD::halfSample(mGrayImageS, mGrayImageS2);
      }
      CVD::img_save (mGrayImageS2, mpEdgeObjectDetection->viewFilename);
      //mpCameraSource->savePointCloud (mpEdgeObjectDetection->pointcloudFilename, xyz1.x, xyz2.x, xyz1.y, xyz2.y);
    }
  }
  mbLearningObject = false;
}

void MainWindow::on_mpFpsSlider_valueChanged(int value)
{
   REAL_TYPE new_max = 1.0f/value;
   std::cerr << "SLIDER VALUE CHANGED " << value << " " << new_max << std::endl;
   mpEdgeObjectDetection->set_maxTime (new_max);
}

void MainWindow::on_mpNextObjButton_clicked()
{
   learnObjNo++;
   std::cerr << "Now learning object number " << learnObjNo << std::endl;
}

void MainWindow::on_mpPrevObjButton_clicked()
{
  learnObjNo--;
  std::cerr << "Previous learning object number " << learnObjNo << std::endl;
}

void MainWindow::on_mpIsLearnCheckBox_stateChanged( int state )
{
   if (state == Qt::Checked)
      isLearn = true;
   else
      isLearn = false;
}

void MainWindow::initialiseCodebook ()
{
    mpEdgeObjectDetection->initialise_codebook();
}


void MainWindow::loadCodebookFromFile (std::string filename)
{
    std::cerr << "entering here with filename " << filename << std::endl;
    QString filenameQ = QString::fromStdString(filename);
    QFile file(filenameQ);
    if( file.exists() )
    {
      std::cerr << "file exists" << std::endl;
      file.close();
      if( mpEdgeObjectDetection->isCodeBookEmpty() )
        mpEdgeObjectDetection->load_codebook( filename );
      else
        mpEdgeObjectDetection->load_codebook( filename, 0); 
      // now we need to load individual codebooks?
      mbDetection = true;
    }
}

int MainWindow::get_memory()
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
    int mem_size = (vm_size/page_size)/1024;
    return mem_size;
}

void MainWindow::setPublishers (ros::Publisher dI_pub)
{
   detectionImage_pub_ = dI_pub;
}

std::vector<DetectedObject> MainWindow::cleanDetectedObjects (std::vector<DetectedObject> mvB4DetectedObjects)
{
   bool whichObjects [mvB4DetectedObjects.size()];
   std::vector<DetectedObject> mvDetectedObjects;
   mvDetectedObjects.reserve (mvB4DetectedObjects.size());
   for (int i = 0; i < mvB4DetectedObjects.size(); i++)
   {
      whichObjects[i] = true;
      for (int j = i+1; j < mvB4DetectedObjects.size(); j++)
      {
         if (mvB4DetectedObjects[i].iObject_ID == mvB4DetectedObjects[j].iObject_ID && mvB4DetectedObjects[i].error > mvB4DetectedObjects[j].error)
         {
            if (getPascalOverlap (mvB4DetectedObjects[i].irTopLeft.x, mvB4DetectedObjects[i].irTopLeft.y, mvB4DetectedObjects[i].irBottomRight.x, mvB4DetectedObjects[i].irBottomRight.y, mvB4DetectedObjects[j].irTopLeft.x, mvB4DetectedObjects[j].irTopLeft.y, mvB4DetectedObjects[j].irBottomRight.x, mvB4DetectedObjects[j].irBottomRight.y) > 0.2)
            {
                whichObjects[i] = false;
                whichObjects[j] = true;
            } 
         }
      }
   }
   for (int i = 0; i < mvB4DetectedObjects.size(); i++)
     if (whichObjects[i]){
       mvDetectedObjects.push_back (mvB4DetectedObjects[i]);
     }
     else{
     }
   return mvDetectedObjects;
}


void MainWindow::on_mpShowMessagesCheckBox_stateChanged( int state )
{
  if( state == Qt::Checked ) mbShowMessages = true;
  else mbShowMessages = false;
}

void MainWindow::on_mpShowEdgesCheckBox_stateChanged( int state )
{
  if( state == Qt::Checked ) mbShowEdges = true;
  else mbShowEdges = false;
}

void MainWindow::on_mpShowChainCheckBox_stateChanged (int state)
{
  if( state == Qt::Checked ) mbShowChains = true;
  else mbShowChains = false;
}

void MainWindow::on_mpShowAllChainsCheckBox_stateChanged (int state)
{
  if( state == Qt::Checked ) mbShowAllChains = true;
  else mbShowAllChains = false;
}

void MainWindow::on_mpShowLinesCheckBox_stateChanged (int state)
{
  if( state == Qt::Checked ) mbShowLines = true;
  else mbShowLines = false;
}
// clearing codebook from all the data

void MainWindow::autoStartDetector ()
{
}

void MainWindow::loadAllGT()
{
   loadGT ("../Dataset_Bristol2010/TestImages_undistorted/plier_gt.txt", 3);
   loadGT ("../Dataset_Bristol2010/TestImages_undistorted/hammer_gt.txt", 2);
   loadGT ("../Dataset_Bristol2010/TestImages_undistorted/box_gt.txt", 5);
   loadGT ("../Dataset_Bristol2010/TestImages_undistorted/stapler_gt.txt", 7);
   loadGT ("../Dataset_Bristol2010/TestImages_undistorted/charger_gt.txt", 10);
   loadGT ("../Dataset_Bristol2010/TestImages_undistorted/driver_gt.txt", 4);
   loadGT ("../Dataset_Bristol2010/TestImages_undistorted/screwdriver_gt.txt", 1);
   loadGT ("../Dataset_Bristol2010/TestImages_undistorted/wood_gt.txt", 6);
   loadGT ("../Dataset_Bristol2010/TestImages_undistorted/wrench_gt.txt", 8);
   loadGT ("../Dataset_Bristol2010/TestImages_undistorted/yellow_gt.txt", 9);

}

void MainWindow::loadGT (char* gt, int no)
{
   char gtFileName[100];
   sprintf (gtFileName, gt);
   std::ifstream if_codebook;
   if_codebook.open(gtFileName);
   std::string line;
   std::getline (if_codebook, line);
   while (line.size() > 0)
   {
     std::stringstream ss (line);
     std::string field;
     std::getline (ss, field, ',');
     std::stringstream fs (field);
     int frameNo = 0, x1, x2, y1, y2;
     fs >> frameNo;
     std::getline (ss, field, ',');
     std::stringstream fs2 (field);
     fs2 >> x1;
     std::getline (ss, field, ',');
     std::stringstream fs3 (field);
     fs3 >> y1;
     std::getline (ss, field, ',');
     std::stringstream fs4 (field);
     fs4 >> x2;
     std::getline (ss, field, ',');
     std::stringstream fs5 (field);
     fs5 >> y2;
     std::getline (if_codebook, line);
     gtItem gti;
     gti.frameNo = frameNo;
     gti.x1 = x1;
     gti.y1 = y1;
     gti.x2 = x2;
     gti.y2 = y2;
     gti.objNo = no;
     allGT.push_back (gti);
   }
}

bool MainWindow::compareToGT (int fNo, int objNo, int x1, int y1, int x2, int y2)
{
   std::cout << "comparing to " << fNo << " " << x1 << " " << y1 << " " << x2 << " " << y2 << std::endl;
   for (int i = 0; i < allGT.size(); i++)
   {
      if (allGT[i].frameNo == fNo && objNo == allGT[i].objNo)
      {
         std::cout << allGT[i].frameNo << " " << allGT[i].x1 << " " << allGT[i].y1 << " " << allGT[i].x2 << " " << allGT[i].y2 << std::endl;
         if (getPascalOverlap (allGT[i], x1, y1, x2, y2) > 0.5)
             return true;
      }
   }
   return false;
}

bool MainWindow::compareToGT2 (int frameNo, int objNo)
{
    if (frameNo <= 13 && objNo == 3)
       return true;
    return false;
}

double MainWindow::getPascalOverlap (gtItem g, int x1, int y1, int x2, int y2)
{
   double outerBB = (max(x2, g.x2)-min (x1, g.x1))*(max(y2,g.y2)-min(y1,g.y1));
   double innerBB = (min(x2, g.x2)-max (x1, g.x1))*(min(y2,g.y2)-max(y1,g.y1));
   return innerBB/outerBB;
}

double MainWindow::getPascalOverlap (int o1_x1, int o1_y1, int o1_x2, int o1_y2, int o2_x1, int o2_y1, int o2_x2, int o2_y2)
{
   double outerBB = (max(o1_x2, o2_x2)-min (o1_x1, o2_x1))*(max(o1_y2,o2_y2)-min(o1_y1,o2_y1));
   double innerBB = (min(o1_x2, o2_x2)-max (o1_x1, o2_x1))*(min(o1_y2,o2_y2)-max(o1_y1,o2_y1));
   return innerBB/outerBB;
}

int MainWindow::max (int a, int b)
{
   if (a > b) return a;
   return b;
}

int MainWindow::min (int a, int b)
{
   if (a < b) return a;
   return b;
}

void MainWindow::keyPressEvent( QKeyEvent *e )
{
  switch( e->key() )
  {
    case Qt::Key_Space:
      ui->mpIsLearnCheckBox->toggle();
      break;
    case Qt::Key_N:
      on_mpNextObjButton_clicked();
      break;
    case Qt::Key_P:
      on_mpPrevObjButton_clicked();
      break;
    default:
      break;
  }
}

void MainWindow::publishResults()
{
        std::vector<DetectedObject> & mvDetectedObjects = mpEdgeObjectDetection->get_detected_objects();
        if (mvDetectedObjects.size() > 0)
        {
             mod_base::DetectedObjectArray obj_array;
             obj_array.header = header_;
             obj_array.objects.resize(mvDetectedObjects.size());
             for (int objNo = 0; objNo < mvDetectedObjects.size(); objNo++)
             {
				obj_array.objects[objNo].objClass = mvDetectedObjects[objNo].iObject_ID;
				 
                obj_array.objects[objNo].roi.x_offset = mvDetectedObjects[objNo].irTopLeft.x*imgSizeRatio;
                obj_array.objects[objNo].roi.y_offset = mvDetectedObjects[objNo].irTopLeft.y*imgSizeRatio;
                obj_array.objects[objNo].roi.height = mvDetectedObjects[objNo].irBottomRight.y*imgSizeRatio - mvDetectedObjects[objNo].irTopLeft.y*imgSizeRatio;
                obj_array.objects[objNo].roi.width = mvDetectedObjects[objNo].irBottomRight.x*imgSizeRatio - mvDetectedObjects[objNo].irTopLeft.x*imgSizeRatio;
                
                obj_array.objects[objNo].points.resize(mvDetectedObjects[objNo].vv2Pnts.size());
                
                for (int ipt = 0; ipt < mvDetectedObjects[objNo].vv2Pnts.size(); ipt++)
                {
				   TooN::Vector<2> point = mvDetectedObjects[objNo].vv2Pnts[ipt];
				   obj_array.objects[objNo].points[ipt].x = point[0];
				   obj_array.objects[objNo].points[ipt].y = point[1];
                }
             }
             objects_publisher_.publish(obj_array);
        }   
}

void MainWindow::printToFile()
{
        std::vector<DetectedObject> & mvDetectedObjects = mpEdgeObjectDetection->get_detected_objects();
        if (mvDetectedObjects.size() > 0)
        {
			for (int objNo = 0; objNo < mvDetectedObjects.size(); objNo++)
            {
				output << header_.seq << " " << header_.stamp << " ";
				output << mvDetectedObjects[objNo].iObject_ID << " ";
				output << (mvDetectedObjects[objNo].irTopLeft.x*imgSizeRatio) << " " << (mvDetectedObjects[objNo].irTopLeft.y*imgSizeRatio) << " " << (mvDetectedObjects[objNo].irBottomRight.x*imgSizeRatio - mvDetectedObjects[objNo].irTopLeft.x*imgSizeRatio) << " " << (mvDetectedObjects[objNo].irBottomRight.y*imgSizeRatio - mvDetectedObjects[objNo].irTopLeft.y*imgSizeRatio) << std::endl;
		    }
		}
}

void MainWindow::Update()
{
  QString foundObjs = "";
  if (isLearn)
      foundObjs += "LEARNING: ";
  // Get image from CameraSource.
  bool newFrame = false;
  if( mbContinueImageFrame || mbNextImageFrame )
  {
      bool haveNew = mpCameraSource->getMonoAndRgbImageAllFrame (mGrayImage, mRgbImage);
      if (!haveNew)
         return;
      newFrame = true;
      header_ = mpCameraSource->getHeader();
      imageRegions.clear();
      imageMasks.clear();
      imageObjectNos.clear();
      std::vector<cv::Rect> thisImageRegions = mpCameraSource->getRegions();
      if (isLearn)
      {
                thisImageRegions.clear();
                thisImageRegions.push_back(cv::Rect (70, 70, 480, 300));
      }
      else
      {
                thisImageRegions.clear();
                thisImageRegions.push_back(cv::Rect (10, 10, 600, 420));
      }
      std::vector<sensor_msgs::Image> thisImageMasks = mpCameraSource->getMasks();
      std::vector<int> thisImageObjects = mpCameraSource->getObjNos();
      std::vector<DetectedObject> & mvDetectedObjects_previous = mpEdgeObjectDetection->get_previous_detected_objects();
      for (int r = 0; r < thisImageRegions.size(); r++)
      {
                SegmentedRegion sr;
                sr.irBottomRight.x = min(thisImageRegions[r].x+thisImageRegions[r].width+regionWidth, beginW.x+imgSize.x) - beginW.x;
                sr.irTopLeft.x = max (thisImageRegions[r].x-regionWidth, beginW.x) - beginW.x;
                sr.irBottomRight.y = min(thisImageRegions[r].y+thisImageRegions[r].height+regionWidth, beginW.y+imgSize.y) - beginW.y;
                sr.irTopLeft.y = max(thisImageRegions[r].y-regionWidth, beginW.y) - beginW.y;
                if (sr.irBottomRight.x > sr.irTopLeft.x && sr.irBottomRight.y > sr.irTopLeft.y)
                {
                   imageRegions.push_back (sr);
                }
      }
      if (!isMask)
      {
                  for (int r = 0; r < imageRegions.size(); r++)
                  {
                     SegmentedRegion sr = imageRegions[r];
              	     // draw four lines representing the region
	             std::vector<TooN::Vector<2> > regionBoundary;
                     TooN::Vector<2> c1, c2, c3, c4;
                     c1[0] = sr.irTopLeft.x;
                     c1[1] = sr.irTopLeft.y;
                     c2[0] = sr.irTopLeft.x;
                     c2[1] = sr.irBottomRight.y;
                     c3[0] = sr.irBottomRight.x;
                     c3[1] = sr.irBottomRight.y;
                     c4[0] = sr.irBottomRight.x;
                     c4[1] = sr.irTopLeft.y; 
                     regionBoundary.push_back (c1);
                     regionBoundary.push_back (c2);
                     regionBoundary.push_back (c3);
                     regionBoundary.push_back (c4);
                     regionBoundary.push_back (c1);
                     ui->mpGLWidget->addFeature(regionBoundary,  QColor(250,250,200), LINE, 1.0f );
                  }
      }
      else
      {
               for (int r = 0; r < thisImageMasks.size(); r++)
               {
                  SegmentedMask sm;
                  imageObjectNos.push_back (thisImageObjects[r]);
                  sm.image = thisImageMasks[r];
                  std::vector<TooN::Vector<2> > maskBoundary;
                  //>> how to draw the mask points
                  for (int mW = 1; mW < sm.image.width-1; mW++)
                     for (int mH = 1; mH < sm.image.height-1; mH++)
                         if (sm.image.data[mH*sm.image.width+mW] > 0)
                         { 
                            if (sm.image.data[(mH-1)*sm.image.width+mW] == 0 || sm.image.data[(mH+1)*sm.image.width+mW] == 0)
                            {
                                TooN::Vector<2> bP;
                                bP[0] = mW;
                                bP[1] = mH;
                                maskBoundary.push_back(bP);
                                continue;
                            }
                            if (sm.image.data[mH*sm.image.width+mW+1] == 0 || sm.image.data[mH*sm.image.width+mW-1] == 0)
                            {
                                TooN::Vector<2> bP;
                                bP[0] = mW;
                                bP[1] = mH;
                                maskBoundary.push_back(bP);
                                continue;
                            }
                            if (sm.image.data[(mH-1)*sm.image.width+mW+1] == 0 || sm.image.data[(mH+1)*sm.image.width+mW+1] == 0)
                            {
                                TooN::Vector<2> bP;
                                bP[0] = mW;
                                bP[1] = mH;
                                maskBoundary.push_back(bP);
                                continue;
                            }
                            if (sm.image.data[(mH-1)*sm.image.width+mW-1] == 0 || sm.image.data[(mH+1)*sm.image.width+mW-1] == 0)
                            {
                                TooN::Vector<2> bP;
                                bP[0] = mW;
                                bP[1] = mH;
                                maskBoundary.push_back(bP);
                                continue;
                            }                            
                         }
                  ui->mpGLWidget->addFeature(maskBoundary,  QColor(215,215,215), POINT, 1.f );
                  imageMasks.push_back (sm);
               }
      }
      mLLImagePyramid.SetImageAtLevelZero( mGrayImage );
      mLLImagePyramid.Reset();
 
      if( mbRunEdgeTracker )
      {
        if( mbDetection && newFrame )
        {
        //  static bool isEdgeObjectDetectionRunning = false;
           // isEdgeObjectDetectionRunning = true;
            mLLImagePyramidForDetection.SetImageAtLevelZero( mGrayImage );
            mLLImagePyramidForDetection.Reset();
            if (imageRegions.size() > 0 || imageMasks.size() > 0)
            {
           //   mpEdgeObjectDetection->preCalculateEdgeLnksWithRegions (mLLImagePyramidForDetection.GetImage(2), imageRegions);
              mpEdgeObjectDetection->preCalculateEdgeLnksWithRegions_MultiScale (mLLImagePyramidForDetection.GetImage(1), imageRegions, imageMasks, mpCameraSource->getSeqNo());
              if (isLearn)
              {
                 if (imageRegions.size() > 0)
                 {
                   for (int r = 0; r < imageRegions.size(); r++)
                   {
                      addViewLearning (learnObjNo, imageRegions[r]);
                   }
                 }
              }

            // checking previous edgelets??
           if (fromPrevious)
           {
            std::vector<DetectedObject> & mvDetectedObjects_previous = mpEdgeObjectDetection->get_previous_detected_objects();
            for (int no_prev_detect = 0; no_prev_detect < mvDetectedObjects_previous.size(); no_prev_detect++)
            {
   
               bool detected = mpEdgeObjectDetection->detect_from_previous (0, no_prev_detect);
                   std::vector<TooN::Vector<2> > prevObjBoundary;
                   TooN::Vector<2> c1, c2, c3, c4;
                   c1[0] = mvDetectedObjects_previous[no_prev_detect].irTopLeft.x * 2;
                   c1[1] = mvDetectedObjects_previous[no_prev_detect].irTopLeft.y * 2;
                   c2[0] = mvDetectedObjects_previous[no_prev_detect].irTopLeft.x * 2;
                   c2[1] = mvDetectedObjects_previous[no_prev_detect].irBottomRight.y * 2;
                   c3[0] = mvDetectedObjects_previous[no_prev_detect].irBottomRight.x * 2;
                   c3[1] = mvDetectedObjects_previous[no_prev_detect].irBottomRight.y * 2;
                   c4[0] = mvDetectedObjects_previous[no_prev_detect].irBottomRight.x * 2;
                   c4[1] = mvDetectedObjects_previous[no_prev_detect].irTopLeft.y * 2; 
                   prevObjBoundary.push_back (c1);
                   prevObjBoundary.push_back (c2);
                   prevObjBoundary.push_back (c3);
                   prevObjBoundary.push_back (c4);
                   prevObjBoundary.push_back (c1);
                   if (detected)
                     ui->mpGLWidget->addFeature(prevObjBoundary, Qt::green, LINE, 2.f );
                   else
                     ui->mpGLWidget->addFeature(prevObjBoundary, Qt::red, LINE, 4.f );

            }
           }
//            if (!mbLearningObject && mpEdgeObjectDetection->detect_constellation(0))
            if (!mbLearningObject && ( (isStar && mpEdgeObjectDetection->detect_star(0)) || (!isStar && mpEdgeObjectDetection->detect(0)) ) )
            {  std::vector<DetectedObject> & mvDetectedObjects = mpEdgeObjectDetection->get_detected_objects();
              
              for( int no_detected_object = 0; no_detected_object < mvDetectedObjects.size(); ++ no_detected_object )
              {
                 int detected_Object_ID = mvDetectedObjects[no_detected_object].iObject_ID;
                 switch(detected_Object_ID)
                 {
                    case 0:
                       if (scenarioNo == 1)
                            foundObjs += "box ";
                       else if (scenarioNo == 2)
                            foundObjs += "bottle ";
                       else
                            foundObjs += "1 ";
                       ui->mpGLWidget->addFeature(mvDetectedObjects[no_detected_object].vv2Pnts, Qt::cyan, POINT, 4.f );
                       break;
                    case 1:
                       if (scenarioNo == 1)
                            foundObjs += "screwdriver ";
                       else if (scenarioNo == 2)
                            foundObjs += "marker pen ";
                       else
                            foundObjs += "2 ";
                      ui->mpGLWidget->addFeature(mvDetectedObjects[no_detected_object].vv2Pnts, Qt::red, POINT, 4.f );
                       break;
                    case 2:
                       if (scenarioNo == 1)
                            foundObjs += "hammer ";
                       else if (scenarioNo == 2)
                            foundObjs += "tape dispenser ";
                       else
                            foundObjs += "3 ";
                      ui->mpGLWidget->addFeature(mvDetectedObjects[no_detected_object].vv2Pnts, Qt::blue, POINT, 4.f );
                       break;
                    case 3:
                       if (scenarioNo == 1)
                            foundObjs += "baton ";
                       else if (scenarioNo == 2)
                            foundObjs += "box ";
                       else
                            foundObjs += "4 ";
                      ui->mpGLWidget->addFeature(mvDetectedObjects[no_detected_object].vv2Pnts, Qt::magenta, POINT, 4.f );
                       break;
                    case 4:
                       if (scenarioNo == 1)
                            foundObjs += "baton ";
                       else if (scenarioNo == 2)
                            foundObjs += "box ";
                       else
                            foundObjs += "5 ";
                      ui->mpGLWidget->addFeature(mvDetectedObjects[no_detected_object].vv2Pnts, Qt::cyan, POINT, 4.f );
                       break;
                    case 5:
                       if (scenarioNo == 1)
                            foundObjs += "baton ";
                       else if (scenarioNo == 2)
                            foundObjs += "box ";
                       else
                            foundObjs += "6 ";
                      ui->mpGLWidget->addFeature(mvDetectedObjects[no_detected_object].vv2Pnts, Qt::yellow, POINT, 4.f );
                       break;
                    case 6:
                       if (scenarioNo == 1)
                            foundObjs += "baton ";
                       else if (scenarioNo == 2)
                            foundObjs += "box ";
                       else
                            foundObjs += "7 ";
                      ui->mpGLWidget->addFeature(mvDetectedObjects[no_detected_object].vv2Pnts, Qt::yellow, POINT, 4.f );
                       break;
                    default:
                        foundObjs += "Object New";
                        ui->mpGLWidget->addFeature(mvDetectedObjects[no_detected_object].vv2Pnts, Qt::cyan, POINT, 4.f );
                 }
                 if (mbShowChains)
                 {
                   ui->mpGLWidget->addFeature(mvDetectedObjects[no_detected_object].foundChainPnts, Qt::green, LINE, 3.f );
                 }

            }
          }
          else
          {
            /*   if (isLearn)
               {
                  std::vector<Edgelet> currentEdgelets = mpEdgeObjectDetection->getEdgelets();
                  std::vector<TooN::Vector<2> > vv2Pnts;
                  for (int i = 0; i < currentEdgelets.size(); i++){
                     TooN::Vector<2> scaledPoint;
                     scaledPoint[0] = currentEdgelets[i].v2Pose[0]*2;
                     scaledPoint[1] = currentEdgelets[i].v2Pose[1]*2;
                     vv2Pnts.push_back (scaledPoint);
                  }
                  ui->mpGLWidget->addFeature(vv2Pnts, Qt::red, POINT, 4.f );
               }*/
          }

        }
        else // mbRun
        {
          if( mbShowEdges )
          {
            mEdges.clear();
            mpEdgeObjectDetection->getEdges (mEdges);
          }
          if( mbShowLSDEdges )
          {
            mEdges2.clear();
            if (imageRegions.size() > 0)
              mpEdgeObjectDetection->getLSDEdgesInRegions (mEdges2, imageRegions);
            else
              mpEdgeObjectDetection->getLSDEdges (mEdges2);
          }
          if( mbShowAllChains )
          {
            checkedChains.clear ();
            checkedConstellations.clear ();
            mpEdgeObjectDetection->getFoundChains (checkedChains);
            mpEdgeObjectDetection->getMatchedConstellations (checkedConstellations);
            std::cerr << "should have read chains here " << checkedChains.size() << std::endl;
          }
          mLLImagePyramid.GetImage(0).ComputeSmoothedImage( 1.6f );
          imRealShareImage = mLLImagePyramid.GetImage(0).mSmoothedImage;
          //          std::cout <<"Tracking Time : " << CVD::timer.get_time() - start_time << std::endl;
        }
        }
        ++miNoFrame;
      }
      if( mbNextImageFrame ) mbNextImageFrame = false;
    }
    else
      usleep( 30000 );
  // Display background image whether undistorted or distorted image.
    if (! mbShowEdges)
       ui->mpGLWidget->Display_Distorted_Image( mRgbImage );
    else
      {
         mpEdgeObjectDetection->getEdgeImage(mEdgeTransformImage);
      //   mpEdgeObjectDetection->getDistanceTransformImage(mEdgeTransformImage);
         ui->mpGLWidget->Display_Distorted_Image(mEdgeTransformImage);
      }
    ui->mpGLWidget->Start_Screen_Coordinate_System( mRgbImage.size(), false );
    if( mbShowEdges )
    {
            mEdges.clear();
            mpEdgeObjectDetection->getEdges (mEdges);
     //  copy (mGrayImage, mGrayImageE);
     //  mCannyEdgeDetector.compute( mGrayImageE, mEdges, mLLImagePyramid.GetImage(0), 8, 0.5, 0.3 );
	
   //    ui->mpGLWidget->Draw_Features( mEdges.virEdges, QColor(249, 249, 180) , POINT, 1.f );
                  ui->mpGLWidget->Draw_Features( mEdges.virEdges, Qt::yellow, POINT, 3.f );
                  std::vector<Edgelet> currentEdgelets = mpEdgeObjectDetection->getEdgelets();
                  std::vector<TooN::Vector<2> > vv2Pnts;
                  for (int i = 0; i < currentEdgelets.size(); i++){
                     TooN::Vector<2> scaledPoint;
                     scaledPoint[0] = currentEdgelets[i].v2Pose[0]*2;
                     scaledPoint[1] = currentEdgelets[i].v2Pose[1]*2;
                     vv2Pnts.push_back (scaledPoint);
                  }
                  ui->mpGLWidget->addFeature(vv2Pnts, Qt::red, POINT, 4.f );
      }
     if( mbShowLSDEdges )
     {
       ui->mpGLWidget->Draw_Features( mEdges2.virEdges, Qt::yellow, POINT, 3.f );
     }
     if( mbShowAllChains )
     {
        QColor vColor[11] = {Qt::black, Qt::red, Qt::darkRed, Qt::darkGreen, Qt::blue, Qt::darkBlue, Qt::darkCyan, Qt::darkMagenta, Qt::white, Qt::darkYellow, Qt::darkGray};
        for( int no_chain = 0; no_chain < checkedChains.size(); ++ no_chain )
        {
	  QColor myColor = vColor[no_chain%11];
          #if 1
          bool outOfCoords = false;
          for (int j = 0; j < 5; j++)
          {
             TooN::Vector<2> v;
             v = checkedChains[no_chain].chains[j];
             if (v[0] > 640 || v[1] > 480)
             {
                outOfCoords = true;
                break;
             }   
          }
          #endif
          ui->mpGLWidget->addFeature(checkedChains[no_chain].chains, myColor, LINE, 1.f );
        }
        checkedConstellations.clear ();
        mpEdgeObjectDetection->getMatchedConstellations (checkedConstellations);
        ui->mpGLWidget->addFeature(checkedConstellations, Qt::cyan, POINT, 3.f );
     }

     ui->mpGLWidget->Draw_Features(); // actual draw


    // ui->mpGLWidget->Draw_Features( mLLImagePyramid.GetImage(0).mvMaxCorners, Qt::green, POINT, 5.f );
    ui->mpGLWidget->Stop_Screen_Coordinate_System();
  ui->mpGLWidget->Start_Screen_Coordinate_System( mRgbImage.size(), false );
  {
//    ui->mpGLWidget->Draw_Features( vEdgeLetsSegments, Qt::red, 1.f ); 
    // Draw text messages.
//    ui->mpGLWidget->Draw_Features( mEdges.virEdges, Qt::red, POINT, 2.f );

    if( mbShowMessages )
    {
      QString framePerSecond;
      framePerSecond.setNum( floor( 0.5f + 1000.f / mTime.elapsed() ), 'f', 2 );
      mTime.start();
      QString view_number_str;
      QString obj_number_str;
      QString const_number_str;
      obj_number_str.setNum( mpEdgeObjectDetection->getNumberOfLearntObjects() );
      view_number_str.setNum( mpEdgeObjectDetection->getNumberOfLearntViews() );
      //const_number_str.setNum( mpEdgeObjectDetection->getNumberOfLearntConstellations() );
      framePerSecond += " fps ";
      if (isLearn)
        framePerSecond += "learnt: " + obj_number_str + " objects " + view_number_str + " views";
      ui->mpGLWidget->qglColor( Qt::red );
      ui->mpGLWidget->renderText( 10, 470, framePerSecond );
      ui->mpGLWidget->qglColor( Qt::yellow );
//      if (!mpEdgeObjectDetection->isViewAdded())
      if (isLearn)
      {
        QFont font;
        font.setPointSize(60);
        ui->mpGLWidget->renderText( 60, 100, foundObjs, font );
      }
      else
      {
        QFont font;
        font.setPointSize(20);
        ui->mpGLWidget->renderText( 15, 20, foundObjs, font );
      }
    }
  }
  
  ui->mpGLWidget->Stop_Screen_Coordinate_System();

  // Expose GL buffer.
  ui->mpGLWidget->updateGL();


  // Save result image.
  if( mbSaveImages )
  {
    ui->mpGLWidget->Get_Drawbuffer( mRgbImageforSaving );
    if(mpSaveImageThread == NULL )
    {
      char stringName [100];
      sprintf (stringName, "/space/space/damen/test1");
      QString test = stringName;
      mpSaveImageThread = new SaveImageThread( CVD::ImageRef(640, 480), test, this );
    }
  //  const CVD::ImageRef smallSize (640,480);
  //  CVD::Image<CVD::Rgb<CVD::byte> > mS;
  //  mS.resize(smallSize);
  //  CVD::halfSample(mRgbImageforSaving, mS);
    mpSaveImageThread->Add_Image_into_Saving_Queue( mRgbImageforSaving );
  }

  ui->mpGLWidget->Get_Drawbuffer( mRgbImageforSaving );

  cv_bridge::CvImagePtr d_image(new cv_bridge::CvImage);
  d_image->encoding = "rgb8";
  CVD::ImageRef imSize = mRgbImageforSaving.size();
  IplImage* imageIp = cvCreateImageHeader (cvSize(imSize.x, imSize.y), 8, 3);
  cvSetData (imageIp, mRgbImageforSaving.data(), imageIp->widthStep);
  cvFlip (imageIp, NULL, 0);
  d_image->image = imageIp;

  if (andPublish)
  {
    detectionImage_pub_.publish (d_image->toImageMsg());
    publishResults();
  }
  if (andPrintToFile)
  {
	  printToFile();
  }


  frameNo2++;
  usleep(10000);
}
