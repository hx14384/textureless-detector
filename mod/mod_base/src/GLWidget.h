/*
  Computer Vision Group, Department of Computer Science. University of      
  Bristol. This code can not be used or copied without express permission.   
*/

#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QtOpenGL/QGLWidget>
#include <cvd/image_ref.h>
#include <gvars3/gvars3.h>
#include <TooN/TooN.h>
#include <cvd/image.h>
#include <cvd/image_io.h>
#include "DefineTypes.h"

class MainWindow;

struct PointerEvent
{
  TooN::Vector<2> v2Pose;
  bool bPress;
};

struct PointerState
{
  TooN::Vector<2> v2Pose;
  bool bDraw;
};

enum Feature_Type
{
  POINT=0,
  LINE,
  LINE_STIPPLE
};

struct Feature
{
  std::vector<TooN::Vector<2> > vv2Pnts;
  QColor colour;
  Feature_Type featureType;
  REAL_TYPE pntSize;
};

struct PointFeature
{
  TooN::Vector<2> v2Pose;
  QColor colour;
  REAL_TYPE pntSize;
};

struct AlphaFeature
{
  std::vector<TooN::Vector<2> > vv2Poses;
  QColor colour;
  REAL_TYPE alpha;
};
/**
 @file GLWidget.h
 @brief GLWidget is a helper class for GL drawing.
*/
class GLWidget : public QGLWidget
{
    Q_OBJECT

public:
    GLWidget(QWidget *parent = 0, size_t width=640, size_t height=480);
//    GLWidget(QWidget *parent = 0, size_t width=1280, size_t height=960);
    ~GLWidget();

    void Display_Distorted_Image(CVD::Image<CVD::Rgb<CVD::byte> > &image);
    void Display_Float_Image(CVD::Image<REAL_TYPE> &image);

    void Begin_Undistort();
    void End_Undistort();
    void Set_Image(CVD::Image<CVD::Rgb<CVD::byte> > &image);
    void Get_Drawbuffer(CVD::Image<CVD::Rgb<CVD::byte> > &image);

    void Start_Screen_Coordinate_System( int width, int height, bool bFlip=true );
    inline void Start_Screen_Coordinate_System( CVD::ImageRef size, bool bFlip=true ) { Start_Screen_Coordinate_System(size.x, size.y, bFlip);};
    void Stop_Screen_Coordinate_System();

    /**
      All Draw_Features functions will not start drawing, they just add all features to its container. After calling Draw_Features,
      GLWidget will draw every features from the container.
    */
    void Draw_Features(std::vector<TooN::Vector<2> >& vv2Pnts, QColor colour, Feature_Type fType, REAL_TYPE pntSize = 2.f);
    void Draw_Features(std::vector<CVD::ImageRef> & virPnts, QColor colour, Feature_Type fType, REAL_TYPE pntSize = 2.f);
    void Draw_Features(std::vector<std::pair<TooN::Vector<2>, TooN::Vector<2> > > & vpv2Lines, QColor colour, REAL_TYPE lnSize = 2.f);
    void Draw_Features();
    QSize minimumSizeHint() const;
    QSize sizeHint() const;
    void resizeWindow(int width, int height);
    void addPointer(TooN::Vector<2>& v2PointerPose);
    void addFeature(std::vector<TooN::Vector<2> >& vv2Pnts, QColor colour, Feature_Type fType, REAL_TYPE pntSize = 2.f);
    void addFeature(std::vector<TooN::Vector<2> >& vv2Pnts, int size, QColor colour, Feature_Type fType, REAL_TYPE pntSize = 2.f);
    void addAlphaFeature(std::vector<TooN::Vector<2> >& vv2Pnts, QColor colour, REAL_TYPE alpha=0.3f);
    void addPointFeature(TooN::Vector<2> & v2Pose, QColor colour, REAL_TYPE pntSize = 2.f);
protected:
    void initializeGL();
    void drawPointer(TooN::Vector<2> & v2PointerPose, GLuint &textureID);
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void clearFeatures();
    void drawPoint(TooN::Vector<2>& v2Pose, QColor& colour, REAL_TYPE& pntSize);
    void drawPoints(std::vector<TooN::Vector<2> >& vv2Pose, QColor& colour, REAL_TYPE& pntSize);
    void drawPoints(std::vector<CVD::ImageRef> & virPose, QColor& colour, REAL_TYPE& pntSize);
    void drawLines(std::vector<TooN::Vector<2> >& vv2Pose, QColor& colour, REAL_TYPE& lnSize, bool bStipple=false);
    void drawLines(std::vector<CVD::ImageRef> & virPose, QColor& colour, REAL_TYPE& lnSize, bool bStipple=false);
    void drawLines(std::vector<std::pair<TooN::Vector<2>, TooN::Vector<2> > > &vpv2Lines, QColor&colour, REAL_TYPE& lnSize);
    void drawAlphas(std::vector<TooN::Vector<2> >& vv2Pose, QColor& colour, REAL_TYPE alpha=0.3f);
private:
    GLuint muiBackgroundTextureID;
    GLuint muiFrameBuffer;
    GLuint muiFrameBufferTextureID;
    GLuint muiFrameBufferTextureResultID;
    GLint miCurrentDrawBuffer;

    int miWidth, miHeight, miOrgWidth, miOrgHeight, miImageWidth, miImageHeight;
    GLuint muiShaderProgram_Undist;
    GLint miTex2DLoc_Undist;
// Camera Parameter
    static const TooN::Vector<6> mvDefaultParams;
    GVars3::gvar3<TooN::Vector<6> > mgvvCameraParams;

// how to undistort image ....
// If a graphic card suport OpenGL 2.0/GLSL, then
// it will use GPU to undistort image otherwise
// it will use a cpu to do that with the following class
    bool mbGPU;
    CVD::Image<CVD::Rgb<CVD::byte> > mRgbImage;
// Mouse/Touch pointer event
    PointerEvent mPointerEvent;
    PointerState mPointerState;
    GLuint mgluRedPointer;
    GLuint mgluGreenPointer;
// Scale
    TooN::Vector<2> mv2WindowScale;
    friend class MainWindow;
// Features
    std::vector<Feature> mvFeatures;
    std::vector<PointFeature> mvPointFeatures;
    std::vector<AlphaFeature> mvAlphaFeatures;
};

#endif
