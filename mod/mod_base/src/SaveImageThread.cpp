/*
  Computer Vision Group, Department of Computer Science. University of      
  Bristol. This code can not be used or copied without express permission.   
*/

#include "SaveImageThread.h"
#include <cvd/utility.h>

SaveImageThread::SaveImageThread( CVD::ImageRef irImageSize, QString savingDirectory, QObject *parent) :
  QThread( parent ), miNoImageFrame( 0 ), mbAbort( false ), mirImageSize(
      irImageSize )
{
  msSavingDirectory = savingDirectory;
  start();
}

SaveImageThread::~SaveImageThread()
{
  mbAbort = true;
  wait();
}

void SaveImageThread::Set_Saving_Directory( QString sSavingDirectory )
{
  msSavingDirectory = sSavingDirectory;
}

void SaveImageThread::Add_Image_into_Saving_Queue( CVD::Image<CVD::Rgb<CVD::byte> > &image )
{
  CVD::Image<CVD::Rgb<CVD::byte> > * image4queue = new CVD::Image<CVD::Rgb<CVD::byte> >( image );
  image4queue->make_unique();
/*  QString imageNo;
  std::cout << "4 " << miNoImageFrame << " " << std::endl;
  imageNo.sprintf("%.5d",miNoImageFrame++);
  QString filename = msSavingDirectory + "/image-" + imageNo + ".jpg";
  std::cout << "5" << std::endl;
  CVD::img_save(*image4queue, filename.toStdString());*/

  mpImageQueue.enqueue( image4queue );
}

void SaveImageThread::run()
{
  std::cout << "a1" << std::endl;
  while(true) 
  {
    if( !mpImageQueue.isEmpty() )
    {
      CVD::Image<CVD::Rgb<CVD::byte> > * image4save = mpImageQueue.dequeue();
      QString imageNo;
      imageNo.sprintf("%.5d",miNoImageFrame++);
      QString filename = msSavingDirectory + "/image-" + imageNo + ".jpg";
      CVD::img_save(*image4save, filename.toStdString());
      delete image4save;
    }
    else
    {
      if( mbAbort )
      return;
      sleep(1);
    }
  }
}

