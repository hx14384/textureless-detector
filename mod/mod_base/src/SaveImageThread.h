/*
  Computer Vision Group, Department of Computer Science. University of      
  Bristol. This code can not be used or copied without express permission.   
*/

#ifndef SAVEIMAGETHREAD_H
#define SAVEIMAGETHREAD_H

#include <QThread>
#include <QQueue>
#include <QString>
#include <cvd/image.h>
#include <cvd/image_ref.h>
#include <cvd/image_io.h>
/**
 @file SaveImageThread.h
 @brief It will run in a separate thread for saving image sequence.
 */
class SaveImageThread: public QThread
{
  Q_OBJECT

  public:
    SaveImageThread( CVD::ImageRef irImageSize, QString Directory, QObject* parent = 0);
    ~SaveImageThread();
    void Set_Saving_Directory( QString sSavingDirectory );
    void Add_Image_into_Saving_Queue( CVD::Image<CVD::Rgb<CVD::byte> > & image );
    int getFrameNo () { return miNoImageFrame; };

  protected:
    void run();
    int miNoImageFrame;
    bool mbAbort;
    QString msSavingDirectory;
    CVD::ImageRef mirImageSize;
    QQueue<CVD::Image<CVD::Rgb<CVD::byte> >*> mpImageQueue;
};

#endif // SAVEIMAGETHREAD_H
