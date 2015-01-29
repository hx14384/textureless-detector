/*
  Computer Vision Group, Department of Computer Science. University of      
  Bristol. This code can not be used or copied without express permission.   
*/

#ifndef DETECTORTHREAD_H
#define DETECTORTHREAD_H

#include <QThread>
#include <QSemaphore>
#include <TooN/TooN.h>
#include "EdgeObjectDetection.h"
/**
 @file DetectorThread.h
 @brief This will make a detector running in a different thread from the main thread.
 */
class DetectorThread: public QThread
{
  Q_OBJECT

  public:
    DetectorThread( EdgeObjectDetection * pDetector, QObject *parent = 0 );
    ~DetectorThread();
    void Start();
    bool IsReady()
    {
      return msemReadyForNextFrame.tryAcquire();
    }
    ;
    bool IsBlockReady()
    {
      msemReadyForNextFrame.acquire();
      return true;
    }
    ;
    void NewFrameIsReady()
    {
      msemWaitForNewFrame.release();
    }
    ;
  protected:
    void run();
    bool mbAbort;
    bool mbFound;
    EdgeObjectDetection * mpDetector;
    QSemaphore msemReadyForNextFrame;
    QSemaphore msemWaitForNewFrame;
};

#endif // DETECTORTHREAD_H
