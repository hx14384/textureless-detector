/*
  Computer Vision Group, Department of Computer Science. University of      
  Bristol. This code can not be used or copied without express permission.   
*/

#include <GL/glew.h>
#include <QtGui>
#include <QtOpenGL>
#include "GLWidget.h"
#include <gvars3/instances.h>
#include <cvd/image_io.h>

using namespace CVD;
using namespace TooN;
using namespace GVars3;

GLWidget::GLWidget(QWidget *parent,size_t width, size_t height) :
    QGLWidget(QGLFormat(QGL::DoubleBuffer|QGL::DirectRendering|QGL::AlphaChannel), parent),
    miWidth(width), miHeight(height), miOrgWidth(width), miOrgHeight(height), mbGPU(true),
    mRgbImage(ImageRef(width,height))
{
  mPointerState.bDraw = false;
  mPointerEvent.bPress = false;
//  setCursor(QCursor(Qt::BlankCursor));
}

GLWidget::~GLWidget()
{
    glDeleteTextures(1, &muiBackgroundTextureID);
}

void GLWidget::initializeGL()
{
    glewInit();

    /*GV2.Register(mgvvCameraParams, "Camera.Parameters", mvDefaultParams, HIDDEN|FATAL_IF_NOT_DEFINED);
       float K[4]={(*mgvvCameraParams)[0], (*mgvvCameraParams)[1], (*mgvvCameraParams)[2], (*mgvvCameraParams)[3]};
       float D[2]={(*mgvvCameraParams)[4], (*mgvvCameraParams)[5]};
    if( K[2] > 200 ) // the image size mostly will be (640,480)
    {*/
      //  miImageWidth = 1280; miImageHeight = 960;
        miImageWidth=640; miImageHeight=480;
    /*}
    else
    {
        miImageWidth = 320; miImageHeight = 240;
        mRgbImage.resize(ImageRef(320,240));
    }*/

    if( !GLEW_VERSION_2_0 )
    {
        std::cerr << "ERROR: OpenGL 2.0 not supported." << std::endl;
        mbGPU = false;
        glViewport(0, 0, miOrgWidth, miOrgHeight);
        glGenTextures(1, &muiBackgroundTextureID);
        glBindTexture(GL_TEXTURE_2D, muiBackgroundTextureID);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, miImageWidth, miImageHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);

        return;
    }

    glViewport(0, 0, miOrgWidth, miOrgHeight);
    glGenTextures(1, &muiBackgroundTextureID);
    glBindTexture(GL_TEXTURE_RECTANGLE_ARB, muiBackgroundTextureID);
    glTexParameterf(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameterf(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, GL_RGB, miImageWidth, miImageHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);


// define Shader Programs
    const char *vsrc =
        "#version 120\n"
        "void main()\n"
        "{\n"
        "    gl_TexCoord[0] = gl_MultiTexCoord0;\n"
        "    gl_Position = ftransform();\n"
        "}\n\0";
    GLuint vshader = glCreateShaderObjectARB(GL_VERTEX_SHADER);
    glShaderSourceARB(vshader, 1, (const GLchar**)&vsrc, NULL);
    glCompileShaderARB(vshader);

    const char *fsrc =
        "#version 120\n"
        "#extension GL_ARB_texture_rectangle : enable\n"
        "uniform vec4 K;\n"
        "uniform vec2 D;\n"
        "uniform sampler2DRect tex2D;\n"
        "void main()\n"
        "{\n"
        "    vec2 X = gl_TexCoord[0].st;\n"
        "    vec2 xn = (X-K.zw)/K.xy;\n"
        "    float r2 = dot(xn,xn);\n"
        "    float r4 = r2*r2;\n"
        "    float coef = 1.0 + D.x*r2 + D.y*r4;\n"
        "    vec2 Xd = vec2(coef*xn.x, coef*xn.y);\n"
        "    vec2 result = (Xd*K.xy  + K.zw);\n"
        "    gl_FragColor = texture2DRect(tex2D, result);\n"
        "}\n\0";
    GLuint fshader = glCreateShaderObjectARB(GL_FRAGMENT_SHADER);
    glShaderSourceARB(fshader, 1, (const GLchar**)&fsrc, NULL);
    glCompileShaderARB(fshader);

    muiShaderProgram_Undist = glCreateProgramObjectARB();
    glAttachObjectARB(muiShaderProgram_Undist, vshader);
    glAttachObjectARB(muiShaderProgram_Undist, fshader);
    glLinkProgramARB(muiShaderProgram_Undist);

    glUseProgramObjectARB(muiShaderProgram_Undist); // enable undistortion shader program.
    GLint KLoc, DLoc;
    //KLoc = glGetUniformLocationARB(muiShaderProgram_Undist, "K");
    //DLoc = glGetUniformLocationARB(muiShaderProgram_Undist, "D");
    miTex2DLoc_Undist = glGetUniformLocationARB(muiShaderProgram_Undist, "tex2D");

   // glUniform4fvARB(KLoc, 1, K);
   // glUniform2fvARB(DLoc, 1, D);

// create off-screen framebuffer
    glGenFramebuffersEXT(1, &muiFrameBuffer);
// create framebuffer texture
// first
    glGenTextures(1, &muiFrameBufferTextureID);
    glBindTexture(GL_TEXTURE_RECTANGLE_ARB, muiFrameBufferTextureID);
    glTexParameterf(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameterf(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, GL_RGB, miImageWidth, miImageHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
// second
    glGenTextures(1, &muiFrameBufferTextureResultID);
    glBindTexture(GL_TEXTURE_RECTANGLE_ARB, muiFrameBufferTextureResultID);
    glTexParameterf(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameterf(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, GL_RGB, miImageWidth, miImageHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

// attach them to the off-screen framebuffer
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, muiFrameBuffer);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_RECTANGLE_ARB, muiFrameBufferTextureID, 0);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT1_EXT, GL_TEXTURE_RECTANGLE_ARB, muiFrameBufferTextureResultID, 0);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

    glUseProgramObjectARB(0); // disable the undistortion shader program.

// Pointer
  	glGenTextures(1, &mgluRedPointer);
	  glBindTexture(GL_TEXTURE_2D, mgluRedPointer);
  	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
	  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
  	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_REPEAT);
	  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_REPEAT);
  	QImage pnt(tr("./iconRed.png"));
  	glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, pnt.width(), pnt.height(), 0, GL_BGRA, GL_UNSIGNED_BYTE, pnt.bits() );
	  gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA, pnt.width(), pnt.height(), GL_BGRA, GL_UNSIGNED_BYTE, pnt.bits());

  	glGenTextures(1, &mgluGreenPointer);
	  glBindTexture(GL_TEXTURE_2D, mgluGreenPointer);
  	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
	  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
	  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_REPEAT);
  	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_REPEAT);
  	pnt.load(tr("./iconGreen.png"));
  	glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, pnt.width(), pnt.height(), 0, GL_BGRA, GL_UNSIGNED_BYTE, pnt.bits() );
  	gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGBA, pnt.width(), pnt.height(), GL_BGRA, GL_UNSIGNED_BYTE, pnt.bits());
   // mv2WindowScale = TooN::makeVector( (double)width()/1280.f, (double)height()/960.f);
    mv2WindowScale = TooN::makeVector( (double)width()/640.f, (double)height()/480.f);
}

void GLWidget::Display_Distorted_Image(CVD::Image<CVD::Rgb<CVD::byte> > &image)
{
    if( !mbGPU )
    {
        glBindTexture(GL_TEXTURE_2D, muiBackgroundTextureID);
        glTexSubImage2D( GL_TEXTURE_2D, 0, 0, 0, image.size().x, image.size().y, GL_RGB, GL_UNSIGNED_BYTE, (void*)image.data());
        glPushMatrix();
        // draw background image
        glLoadIdentity();
        glDisable(GL_BLEND);
        glColor3f( 1.f, 1.f, 1.f );
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
        glEnable(GL_TEXTURE_2D);
        glBegin( GL_QUADS );
        {
            glTexCoord2i(0,0); 	glVertex2i(-1,1);
            glTexCoord2i(0,1); 	glVertex2i(-1,-1);
            glTexCoord2i(1,1); 	glVertex2i(1,-1);
            glTexCoord2i(1,0); 	glVertex2i(1,1);
        }
        glEnd();
        glDisable(GL_TEXTURE_2D);
        // end draw background image
        glPopMatrix();
        return;
    }
    qglColor(Qt::white);
    glBindTexture(GL_TEXTURE_RECTANGLE_ARB, muiBackgroundTextureID);
    glTexSubImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, 0, 0, image.size().x, image.size().y, GL_RGB, GL_UNSIGNED_BYTE, (void*)image.data());
    glViewport(0, 0, miWidth, miHeight);
    glEnable(GL_TEXTURE_RECTANGLE_ARB);
    glBegin( GL_QUADS );
    {
      glTexCoord2i(0,0);                      glVertex2i(-1,1);
      glTexCoord2i(0,miImageHeight);            glVertex2i(-1,-1);
      glTexCoord2i(miImageWidth,miImageHeight);   glVertex2i(1,-1);
      glTexCoord2i(miImageWidth,0);             glVertex2i(1,1);
    }
    glEnd();
    glDisable(GL_TEXTURE_RECTANGLE_ARB);
}

void GLWidget::Display_Float_Image(CVD::Image<REAL_TYPE> &image)
{
    if( !mbGPU )
    {
        glBindTexture(GL_TEXTURE_2D, muiBackgroundTextureID);
        glTexSubImage2D( GL_TEXTURE_2D, 0, 0, 0, image.size().x, image.size().y, GL_RGB, GL_UNSIGNED_BYTE, (void*)image.data());
        glPushMatrix();
        // draw background image
        glLoadIdentity();
        glDisable(GL_BLEND);
        glColor3f( 1.f, 1.f, 1.f );
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
        glEnable(GL_TEXTURE_2D);
        glBegin( GL_QUADS );
        {
            glTexCoord2i(0,0); 	glVertex2i(-1,1);
            glTexCoord2i(0,1); 	glVertex2i(-1,-1);
            glTexCoord2i(1,1); 	glVertex2i(1,-1);
            glTexCoord2i(1,0); 	glVertex2i(1,1);
        }
        glEnd();
        glDisable(GL_TEXTURE_2D);
        // end draw background image
        glPopMatrix();
        return;
    }
    qglColor(Qt::white);
    glBindTexture(GL_TEXTURE_RECTANGLE_ARB, muiBackgroundTextureID);
    glTexSubImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, 0, 0, image.size().x, image.size().y, GL_RGB, GL_UNSIGNED_BYTE, (void*)image.data());
    glViewport(0, 0, miWidth, miHeight);
    glEnable(GL_TEXTURE_RECTANGLE_ARB);
    glBegin( GL_QUADS );
    {
      glTexCoord2i(0,0);                      glVertex2i(-1,1);
      glTexCoord2i(0,miImageHeight);            glVertex2i(-1,-1);
      glTexCoord2i(miImageWidth,miImageHeight);   glVertex2i(1,-1);
      glTexCoord2i(miImageWidth,0);             glVertex2i(1,1);
    }
    glEnd();
    glDisable(GL_TEXTURE_RECTANGLE_ARB);
}

void GLWidget::Begin_Undistort()
{
    if( !mbGPU )
        return;
    qglColor(Qt::white);
    glGetIntegerv(GL_DRAW_BUFFER, &miCurrentDrawBuffer);
// enable framebuffer
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, muiFrameBuffer);
    glDrawBuffer(GL_COLOR_ATTACHMENT0_EXT);
}

void GLWidget::End_Undistort()
{
    if( !mbGPU )
        return;
// disable framebuffer
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
// draw screen buffer
// enable shader
    glUseProgramObjectARB(muiShaderProgram_Undist);
    glUniform1iARB(miTex2DLoc_Undist,0);

    glDrawBuffer(miCurrentDrawBuffer);
    glBindTexture(GL_TEXTURE_RECTANGLE_ARB, muiFrameBufferTextureResultID);
    glViewport(0, 0, miWidth, miHeight);
    glEnable(GL_TEXTURE_RECTANGLE_ARB);
    glBegin( GL_QUADS );
    {
        glTexCoord2i(0,0);                      glVertex2i(-1,-1);
        glTexCoord2i(0,miImageHeight);            glVertex2i(-1,1);
        glTexCoord2i(miImageWidth,miImageHeight); 	glVertex2i(1,1);
        glTexCoord2i(miImageWidth,0);             glVertex2i(1,-1);
    }
    glEnd();
    glDisable(GL_TEXTURE_RECTANGLE_ARB);
// disable shader
    glUseProgramObjectARB(0);
}

void GLWidget::Set_Image(CVD::Image<CVD::Rgb<CVD::byte> > &mRgbImage)
{
    if( !mbGPU )
    {
        glBindTexture(GL_TEXTURE_2D, muiBackgroundTextureID);
        glTexSubImage2D( GL_TEXTURE_2D, 0, 0, 0, mRgbImage.size().x, mRgbImage.size().y, GL_RGB, GL_UNSIGNED_BYTE, (void*)mRgbImage.data());
        glPushMatrix();
        // draw background image
        glLoadIdentity();
        glDisable(GL_BLEND);
        glColor3f( 1.f, 1.f, 1.f );
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
        glEnable(GL_TEXTURE_2D);
        glBegin( GL_QUADS );
        {
            glTexCoord2i(0,0); 	glVertex2i(-1,1);
            glTexCoord2i(0,1); 	glVertex2i(-1,-1);
            glTexCoord2i(1,1); 	glVertex2i(1,-1);
            glTexCoord2i(1,0); 	glVertex2i(1,1);
        }
        glEnd();
        glDisable(GL_TEXTURE_2D);
        // end draw background image
        glPopMatrix();
        return;
    }
    glBindTexture(GL_TEXTURE_RECTANGLE_ARB, muiFrameBufferTextureID);
    glTexSubImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, 0, 0, mRgbImage.size().x, mRgbImage.size().y, GL_RGB, GL_UNSIGNED_BYTE, (void*)mRgbImage.data());
    glDrawBuffer(GL_COLOR_ATTACHMENT1_EXT);
    glViewport(0,0,miImageWidth,miImageHeight);
    glEnable(GL_TEXTURE_RECTANGLE_ARB);
    glBegin( GL_QUADS );
    {
        glTexCoord2i(0,0);                      glVertex2i(-1,1);
        glTexCoord2i(0,miImageHeight);            glVertex2i(-1,-1);
        glTexCoord2i(miImageWidth,miImageHeight); 	glVertex2i(1,-1);
        glTexCoord2i(miImageWidth,0);             glVertex2i(1,1);
    }
    glEnd();
    glDisable(GL_TEXTURE_RECTANGLE_ARB);
}

void GLWidget::Get_Drawbuffer(CVD::Image<CVD::Rgb<CVD::byte> > &image)
{
    GLint idCurrentDrawbuffer;
    glGetIntegerv(GL_DRAW_BUFFER, &idCurrentDrawbuffer);
    glReadBuffer(idCurrentDrawbuffer);
    glReadPixels(0, 0, miOrgWidth, miOrgHeight, GL_RGB, GL_UNSIGNED_BYTE, (void*)image.data());
}

QSize GLWidget::minimumSizeHint() const
{
  //return QSize(1280, 960);
  return QSize(640,480);
}

QSize GLWidget::sizeHint() const
{
 // return QSize(1280, 960);
  return QSize(640,480);
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
  mPointerEvent.bPress = true;
  mPointerEvent.v2Pose = TooN::makeVector( (double)(event->pos().x())/mv2WindowScale[0], (double)(event->pos().y())/mv2WindowScale[1] );
  QWidget::mousePressEvent(event);
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
  mPointerEvent.v2Pose = TooN::makeVector((double)(event->pos().x())/mv2WindowScale[0], (double)(event->pos().y())/mv2WindowScale[1]);
  QWidget::mouseMoveEvent(event);
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event)
{
    mPointerEvent.bPress = false;
    QWidget::mouseReleaseEvent(event);
}

void GLWidget::Start_Screen_Coordinate_System( int width, int height, bool bFlip )
{
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    if( bFlip )
        glOrtho(0, width, 0, height, 0.0, -1.0);
    else
        glOrtho(0, width, height, 0, 0.0, -1.0);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
}

void GLWidget::Stop_Screen_Coordinate_System()
{
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}

void GLWidget::Draw_Features(std::vector<std::pair<TooN::Vector<2>, TooN::Vector<2> > > & vpv2Lines, QColor colour, REAL_TYPE lnSize )
{
  glEnable( GL_BLEND );
  glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
  drawLines( vpv2Lines, colour, lnSize );
  glDisable( GL_BLEND );
}

void GLWidget::Draw_Features(std::vector<TooN::Vector<2> >& vv2Pnts, QColor colour, Feature_Type fType, REAL_TYPE pntSize)
{
  glEnable( GL_BLEND );
  glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
  switch( fType )
  {
    case POINT:
      drawPoints(vv2Pnts, colour, pntSize);
      break;
    case LINE:
      drawLines(vv2Pnts, colour, pntSize, false);
      break;
    case LINE_STIPPLE:
      drawLines(vv2Pnts, colour, pntSize, true);
      break;
    default:
      break;
  }
  glDisable( GL_BLEND );
}

void GLWidget::Draw_Features(std::vector<CVD::ImageRef>& virPnts, QColor colour, Feature_Type fType, REAL_TYPE pntSize)
{
  glEnable( GL_BLEND );
  glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
  switch( fType )
  {
    case POINT:
      drawPoints(virPnts, colour, pntSize);
      break;
    case LINE:
      drawLines(virPnts, colour, pntSize, false);
      break;
    case LINE_STIPPLE:
      drawLines(virPnts, colour, pntSize, true);
      break;
    default:
      break;
  }
  glDisable( GL_BLEND );
}

void GLWidget::Draw_Features()
{
  glEnable( GL_BLEND );
  glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );    
// Draw Pointer
  if( mPointerState.bDraw )
    drawPointer(mPointerEvent.v2Pose, mgluGreenPointer);
// Draw Point/Line Features
  if( !mvFeatures.empty() )
  {
    std::vector<Feature>::iterator ftIter = mvFeatures.begin();
    std::vector<Feature>::iterator endIter = mvFeatures.end();
    for( ; ftIter != endIter; ++ftIter )
    {
      switch( (*ftIter).featureType )
      {
        case POINT:
          drawPoints((*ftIter).vv2Pnts, (*ftIter).colour, (*ftIter).pntSize);
          break;
        case LINE:
          drawLines((*ftIter).vv2Pnts, (*ftIter).colour, (*ftIter).pntSize, false);
          break;
        case LINE_STIPPLE:
          drawLines((*ftIter).vv2Pnts, (*ftIter).colour, (*ftIter).pntSize, true);
          break;
        default:
          break;
      }
    }
  }
  if( !mvPointFeatures.empty() ) 
  {
    std::vector<PointFeature>::iterator ptIter = mvPointFeatures.begin();
    std::vector<PointFeature>::iterator endIter = mvPointFeatures.end();
    for( ; ptIter != endIter; ++ptIter )
    {
      drawPoint( (*ptIter).v2Pose, (*ptIter).colour, (*ptIter).pntSize );
    }
  }
  if( !mvAlphaFeatures.empty() )
  {
    std::vector<AlphaFeature>::iterator ptIter = mvAlphaFeatures.begin();
    std::vector<AlphaFeature>::iterator endIter = mvAlphaFeatures.end();
    for( ; ptIter != endIter; ++ptIter )
    {
      drawAlphas( (*ptIter).vv2Poses, (*ptIter).colour, (*ptIter).alpha );
    }
  }
  glDisable( GL_BLEND );
  clearFeatures();
}

void GLWidget::addPointer( TooN::Vector<2>& v2PointerPose )
{
  mPointerState.v2Pose = v2PointerPose;
  mPointerState.bDraw = true;
}

void GLWidget::addAlphaFeature(std::vector<TooN::Vector<2> >& vv2Pnts, QColor colour, REAL_TYPE alpha)
{
  AlphaFeature alphaFeature;
  alphaFeature.vv2Poses = vv2Pnts;
  alphaFeature.colour = colour;
  alphaFeature.alpha = alpha;
  mvAlphaFeatures.push_back(alphaFeature);
}

void GLWidget::clearFeatures()
{
  mvFeatures.clear();
  mvPointFeatures.clear();
  mvAlphaFeatures.clear();
  mPointerState.bDraw = false;
}

void GLWidget::resizeWindow(int width, int height)
{
  resize(width, height);
//  mv2WindowScale = TooN::makeVector( (double)width/1280.f, (double)height/960.f);
  mv2WindowScale = TooN::makeVector( (double)width/640.f, (double)height/480.f);
  glViewport(0, 0, width, height);
}

void GLWidget::addFeature(std::vector<TooN::Vector<2> >& vv2Pnts, QColor colour, Feature_Type fType, REAL_TYPE pntSize)
{
  mvFeatures.resize( mvFeatures.size() + 1 );
  mvFeatures.back().vv2Pnts = vv2Pnts;
  mvFeatures.back().colour = colour;
  mvFeatures.back().featureType = fType;
  mvFeatures.back().pntSize = pntSize;
/*
  Feature feature;
  feature.vv2Pnts = vv2Pnts;
  feature.colour = colour;
  feature.featureType = fType;
  feature.pntSize = pntSize;
  mvFeatures.push_back(feature);
*/  
}

void GLWidget::addFeature(std::vector<TooN::Vector<2> >& vv2Pnts, int size, QColor colour, Feature_Type fType, REAL_TYPE pntSize)
{
  mvFeatures.resize( mvFeatures.size() + 1 );
  mvFeatures.back().vv2Pnts.resize( size );
  std::copy( vv2Pnts.begin(), vv2Pnts.begin() + size, mvFeatures.back().vv2Pnts.begin() );
  mvFeatures.back().colour = colour;
  mvFeatures.back().featureType = fType;
  mvFeatures.back().pntSize = pntSize;
/*
  Feature feature;
  feature.vv2Pnts = vv2Pnts;
  feature.colour = colour;
  feature.featureType = fType;
  feature.pntSize = pntSize;
  mvFeatures.push_back(feature);
*/  
}

void GLWidget::addPointFeature(TooN::Vector<2> & v2Pose, QColor colour, REAL_TYPE pntSize)
{
  PointFeature ptFeature;
  ptFeature.v2Pose = v2Pose;
  ptFeature.colour = colour;
  ptFeature.pntSize = pntSize;
  mvPointFeatures.push_back(ptFeature);
}

void GLWidget::drawPointer( TooN::Vector<2> & v2PointerPose, GLuint & textureID )
{
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, textureID);
//	glEnable(GL_BLEND);
//	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glColor3f(1,1,1);
	glBegin(GL_QUADS);
	{
		glTexCoord2f(0.f,0.f); glVertex2f(mv2WindowScale[0]*(v2PointerPose[0]-10.f),mv2WindowScale[1]*(v2PointerPose[1]-10.f));
		glTexCoord2f(0.f,1.f); glVertex2f(mv2WindowScale[0]*(v2PointerPose[0]-10.f),mv2WindowScale[1]*(v2PointerPose[1]+10.f));
		glTexCoord2f(1.f,1.f); glVertex2f(mv2WindowScale[0]*(v2PointerPose[0]+10.f),mv2WindowScale[1]*(v2PointerPose[1]+10.f));
		glTexCoord2f(1.f,0.f); glVertex2f(mv2WindowScale[0]*(v2PointerPose[0]+10.f),mv2WindowScale[1]*(v2PointerPose[1]-10.f));
	}
	glEnd();
//	glDisable(GL_BLEND);
	glDisable(GL_TEXTURE_2D);
}

void GLWidget::drawPoint(TooN::Vector<2>& v2Pose, QColor& colour, REAL_TYPE& pntSize)
{
  qglColor(colour);
  glEnable( GL_POINT_SMOOTH );
  glPointSize( pntSize );
  glBegin( GL_POINTS );
  {
    glVertex2f( v2Pose[0], v2Pose[1] );
  }
  glEnd();
  glDisable( GL_POINT_SMOOTH );
}

void GLWidget::drawPoints(std::vector<TooN::Vector<2> >& vv2Pnts, QColor& colour, REAL_TYPE& pntSize)
{
  qglColor(colour);
  glEnable( GL_POINT_SMOOTH );
  glPointSize( pntSize );
  glBegin( GL_POINTS );
  {
    for( std::vector<TooN::Vector<2> >::iterator pntIter = vv2Pnts.begin(); pntIter != vv2Pnts.end(); ++pntIter )
      glVertex2f((*pntIter)[0], (*pntIter)[1]);
  }
  glEnd();
  glDisable( GL_POINT_SMOOTH );
}

void GLWidget::drawPoints(std::vector<CVD::ImageRef>& virPnts, QColor& colour, REAL_TYPE& pntSize)
{
  qglColor(colour);
  glEnable( GL_POINT_SMOOTH );
  glPointSize( pntSize );
  glBegin( GL_POINTS );
  {
    for( std::vector<CVD::ImageRef>::iterator pntIter = virPnts.begin(); pntIter != virPnts.end(); ++pntIter )
      glVertex2i((*pntIter).x, (*pntIter).y);
  }
  glEnd();
  glDisable( GL_POINT_SMOOTH );
}

void GLWidget::drawLines(std::vector<std::pair<TooN::Vector<2>, TooN::Vector<2> > > & vpv2Lines, QColor &colour, REAL_TYPE& lnSize)
{
  qglColor( colour );
  glLineWidth( lnSize );
  glBegin(GL_LINES);
  {
    std::vector<std::pair<TooN::Vector<2>, TooN::Vector<2> > >::iterator pntIter = vpv2Lines.begin();
    std::vector<std::pair<TooN::Vector<2>, TooN::Vector<2> > >::iterator endIter = vpv2Lines.end();
    for( ; pntIter != endIter; ++pntIter )
    {
      TooN::Vector<2> &First = (*pntIter).first;
      TooN::Vector<2> &Second = (*pntIter).second;
      glVertex2f(First[0], First[1]);
      glVertex2f(Second[0], Second[1]);
    }
  }
  glEnd();
}

void GLWidget::drawLines(std::vector<TooN::Vector<2> >& vv2Pnts, QColor& colour, REAL_TYPE& lnSize, bool bStipple)
{
  qglColor( colour );
  glLineWidth(lnSize);
  if( bStipple )
  {
    glEnable( GL_LINE_STIPPLE );
    glLineStipple( 1, 0x0101 );
  }
  glBegin(GL_LINE_STRIP);
  {
    for( std::vector<TooN::Vector<2> >::iterator pntIter = vv2Pnts.begin(); pntIter != vv2Pnts.end(); ++pntIter )
      glVertex2f((*pntIter)[0], (*pntIter)[1]);
  }
  glEnd();
  if( bStipple )
    glDisable( GL_LINE_STIPPLE );
}

void GLWidget::drawLines(std::vector<CVD::ImageRef>& virPnts, QColor& colour, REAL_TYPE& lnSize, bool bStipple)
{
  qglColor( colour );
  glLineWidth(lnSize);
  if( bStipple )
  {
    glEnable( GL_LINE_STIPPLE );
    glLineStipple( 1, 0x0101 );
  }
  glBegin(GL_LINE_STRIP);
  {
    for( std::vector<CVD::ImageRef>::iterator pntIter = virPnts.begin(); pntIter != virPnts.end(); ++pntIter )
      glVertex2i((*pntIter).x, (*pntIter).y);
  }
  glEnd();
  if( bStipple )
    glDisable( GL_LINE_STIPPLE );
}

void GLWidget::drawAlphas(std::vector<TooN::Vector<2> >&vv2Poses, QColor& colour, REAL_TYPE alpha)
{
  colour.setAlpha(alpha*255);
  qglColor(colour);
  glBegin(GL_POLYGON);
  {
    std::vector<TooN::Vector<2> >::iterator pntIter = vv2Poses.begin();
    std::vector<TooN::Vector<2> >::iterator endIter = vv2Poses.end();
    for(; pntIter != endIter; ++pntIter )
    {
      glVertex2f((*pntIter)[0], (*pntIter)[1]);
    }
  }
  glEnd();
}

const TooN::Vector<6> GLWidget::mvDefaultParams = TooN::makeVector( -380.0, -380.0, 320.0, 240.0, 0.0, 0.0);
