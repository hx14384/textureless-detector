/********************************************************************************/
/**  Computer Vision Group, Department of Computer Science. University of      **/
/**  Bristol. This code can not be used or copied without express permission.  **/ 
/********************************************************************************/
#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTime>
#include "CameraSource.h"
#include <cvd/image.h>
#include <cvd/image_io.h>
#include "LLImagePyramid.h"
#include "SaveImageThread.h"
#include "DefineTypes.h"
#include "EdgeObjectDetection.h"
#include "CannyEdgeDetector.h"
#include "EdgeExtraction.h"

//#include <sensor_msgs/CameraInfo.h>
//#include <sensor_msgs/PointCloud2.h>
//#include <ros/ros.h>

class QTimer;

struct CreateState
{
  enum Type
  {
    Vertex = 0, Polygon, Rectangle, Circle, Cylinder, Extrude, VirtualPlane
  };
};

struct InputType
{
  enum Type
  {
    None = 0, First, Second
  };
};

struct InputState
{
  enum Type
  {
    Create = 0, AdjustPoint, SecondClicked, MovePointAlongEpipolarLine, MovePoint, SelectPoint, SelectPlane
  };
};

struct MouseEventState
{
  enum Type
  {
    Press = 0, PrePress, Move, Release, PreRelease
  };
};

struct MouseButton
{
  enum Type
  {
    None = 0, Left, Middle, Right, Wheel
  };
};

struct MouseEvent
{
    MouseButton::Type button;
    double scale;
    MouseEventState::Type state;
    bool bMouseMove;
};

struct Select3DReturn
{
  enum Type
  {
    NoThing = 0, BeingBuildModelPoint, BuildModelPoint
  };
};

namespace Ui
{
  class MainWindow;
}
;

class MainWindow: public QMainWindow
{
  Q_OBJECT
  public:
    MainWindow( CameraSource * pCameraSource, char* outputFolderName, char* outputFileName, double maxTime, char* gtFileName, int w_x, int w_y, bool kinect, bool gaze, bool firewire, bool isStar = false, bool lsd = true, QWidget *parent = 0 );
    ~MainWindow();
    void publishResults();
    void printToFile ();
    void setPublishers (ros::Publisher detectionImage_pub_);
    void setCodebookFilename (std::string fn) {codebookFilename = fn; };
    void setDatasetFolder (std::string df);
    void setDepth (bool d) {isDepth = d; mpEdgeObjectDetection->setDepth (d);}
    void setLearn (bool l, int objNo) {isLearn = l; learnObjNo = objNo;}
    void setObjsPublisher (ros::Publisher op) {objects_publisher_ = op;}
    void setMask (bool m) {isMask = m; mpEdgeObjectDetection->setMask (m);}
    void setService (bool aS) {asService = aS;}
    void setPublishing (bool aP) {andPublish = aP;}
    void setStar (bool iS) {isStar = iS;}
    void setPrintToFile (std::string filename);
    void setBackground (bool wB) {withBG = wB;}
    void setCameraKinect (bool iK) { if (iK) isKinect = true; else isKinect = false; }
    void setUndistort (bool dis) { mbUndistortImage = dis; }
    void setScenarioNo (int sNo) {scenarioNo = sNo;}
    void addViewLearning (int objNo, SegmentedRegion region, SegmentedMask mask);
    void addViewLearning (int objNo, SegmentedRegion region);
    void loadCodebookFromFile (std::string filename);
    void initialiseCodebook ();
    void autoStartDetector();

    EdgeObjectDetection *mpEdgeObjectDetection;

    private Q_SLOTS:
    void Update();
    void on_mpShowMessagesCheckBox_stateChanged( int state );
    void on_mpShowEdgesCheckBox_stateChanged( int state );
    void on_mpIsLearnCheckBox_stateChanged (int state);
    void on_mpShowChainCheckBox_stateChanged (int state );
    void on_mpShowAllChainsCheckBox_stateChanged (int state );
    void on_mpShowLinesCheckBox_stateChanged (int state );
    void on_mpFpsSlider_valueChanged( int value ); 

    // Button
    void loadGT (char* gtFileName, int objNo);
    void loadAllGT ();
    bool compareToGT (int fNo, int objNo, int x1, int y1, int x2, int y2);
    bool compareToGT2 (int fNo, int objNo);
    int max (int a, int b);
    int min (int a, int b);

    double getPascalOverlap (gtItem g, int x1, int y1, int x2, int y2);
    double getPascalOverlap (int o1_x1, int o1_y1, int o1_x2, int o1_y2, int o2_x1, int o2_y1, int o2_x2, int o2_y2);
    std::vector<DetectedObject> cleanDetectedObjects (std::vector<DetectedObject> mvB4DetectedObjects);
    int get_memory ();
    
    void on_mpNextObjButton_clicked();
    void on_mpPrevObjButton_clicked();

  protected:
    void keyPressEvent( QKeyEvent *e );
    std_msgs::Header header_;

  private:
    Ui::MainWindow *ui;
    CameraSource * mpCameraSource;
    CVD::Image<CVD::Rgb<CVD::byte> > mRgbImage;
    CVD::Image<CVD::Rgb<CVD::byte> > mRgbImageF;
    CVD::Image<CVD::Rgb<CVD::byte> > mRgbImageforSaving;
    CVD::Image<CVD::byte> mGrayImage;
    CVD::Image<CVD::byte> mGrayImageE;
    CVD::Image<CVD::byte> mGrayImageF;
    CVD::Image<CVD::byte> mGrayImageDummy;
    CVD::Image<CVD::Rgb<CVD::byte> > mEdgeTransformImage;
    TooN::SE3<> mse3W2C;
    LLImagePyramid<CVD::byte> mLLImagePyramid;
    LLImagePyramid<CVD::byte> mLLImagePyramidForDetection_TemplateImage;
    QTimer *mpTimer; // Timer for calling Update()
    QTime mTime; // For measuring time interval
    int miNoFrame;
    // control flags
    bool mbUndistortImage;
    bool mbPublishers;
    bool mbSaveImages;
    bool mbShowMessages;
    bool mbShowEdges;
    bool mbShowLSDEdges;
    bool mbShowChains;
    bool mbShowAllChains;
    bool mbShowLines;
    bool mbRunEdgeTracker;
    bool mbContinueImageFrame;
    bool mbNextImageFrame;
    bool mbSnapToEdges;

    ros::Publisher objects_publisher_;
    bool mbLearningObject;
    bool isKinect;
    bool isDepth;
    bool isLearn;
    bool learnAndSave;
    bool isMask;
    bool withBG;
    int learnObjNo;
    int scenarioNo;
    bool asService;
    // threads
    SaveImageThread * mpSaveImageThread;
    // NEW WIREFRAME CLASS
    std::vector<TooN::Vector<2> > mvv2ImagePointFeatures;
    // Edge Tracker
    GVars3::gvar3<REAL_TYPE> mgvrEdgeThreshold;
    static const REAL_TYPE mgvrDefaultParams;
    // Edge Object Detection
    Edges mEdges;
    Edges mEdges2;
    std::vector<CheckedChains> checkedChains;
    std::vector< TooN::Vector<2> > checkedConstellations;
    std::vector<SegmentedRegion> imageRegions;
    std::vector<SegmentedMask> imageMasks;
    std::vector<int> imageObjectNos;
    EdgeExtraction *mpEdgeExtraction;
    bool mbDetection;
    bool kinect;
    bool andPublish;
    bool andPrintToFile;
    bool isStar;
    bool gaze;
    bool firewire;
    int miDetectedObjectNo;
    QString mObjectFilename;
    QString mObject4DetectionFilename;
    CannyEdgeDetector mCannyEdgeDetector;

    ros::Publisher detectionImage_pub_;

    std::string codebookFilename;
    int imgSizeRatio;
    int regionWidth;
    CVD::ImageRef imgSize;
    
    std::ofstream output;

    int frameNo2;

};
#endif // MAINWINDOW_H
