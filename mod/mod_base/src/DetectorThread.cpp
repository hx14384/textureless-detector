/*
  Computer Vision Group, Department of Computer Science. University of      
  Bristol. This code can not be used or copied without express permission.   
*/

#include "DetectorThread.h"
extern LLImagePyramid<CVD::byte> mLLImagePyramidForDetection;

DetectorThread::DetectorThread( EdgeObjectDetection * pDetector, QObject *parent ) :
  QThread( parent ), mbAbort( false ), mpDetector( pDetector )
{
}

DetectorThread::~DetectorThread()
{
  mbAbort = true;
  msemWaitForNewFrame.release();
  wait();
}

void DetectorThread::Start()
{
  start();
}

void DetectorThread::run()
{
  forever 
  {
    msemWaitForNewFrame.acquire();
    if( !mbAbort )
    {
      mpDetector->preCalculateEdgeLnks( mLLImagePyramidForDetection.GetImage(1), 0);
      mpDetector->detect(0);
    }
    msemReadyForNextFrame.release();
    if( mbAbort )
      return;
  }
}
