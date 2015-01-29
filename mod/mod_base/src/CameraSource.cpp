/*
  Computer Vision Group, Department of Computer Science. University of      
  Bristol. This code can not be used or copied without express permission.   
*/

#include "CameraSource.h"

#include <ros/ros.h>
#include <boost/thread/mutex.hpp>
#include <sensor_msgs/Image.h>
#include <sensor_msgs/image_encodings.h>
#include <sensor_msgs/CameraInfo.h>
#include <cv_bridge/cv_bridge.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <iostream>
#include <cvd/colourspacebuffer.h>
#include <cvd/videosource.h>
#include <cvd/vision.h>
#include <cvd/Linux/v4lcontrol.h>
#include <cvd/Linux/v4lbuffer.h>
#include <cvd/Linux/dvbuffer3.h>
#include <cvd/image_convert.h>

#include <sys/mman.h>
#include <sys/stat.h>        /* For mode constants */
#include <fcntl.h>           /* For O_* constants */
#include <algorithm>
#include <cstring>
#include <stdio.h>

int reachedLearning = 0;

CameraSource::CameraSource() :
  mbIsConnected( false ), CAMERA_READ ( false )
{
}

CameraSource::~CameraSource()
{

}

void CameraSource::rgb_image_msg_cb2 (const sensor_msgs::ImageConstPtr& msg)
{
   m.lock();
   rgb_image_msg = msg;
   m.unlock();
}

void CameraSource::connect( const std::string& source, int argc, char *argv[])
{
  std::cerr << "in connection: " << source << std::endl;
  if( !mbIsConnected )
  {
    mbIsConnected = true;

    imgSub.subscribe(node_handle, colorImgSrc, 40);
    imgSub.registerCallback(&CameraSource::rgb_image_msg_cb2, this);
  }
}

void CameraSource::disconnect()
{
  mbIsConnected = false;
}

std::vector<cv::Rect> CameraSource::getRegions()
{
   return regions;
}

std::vector<sensor_msgs::Image> CameraSource::getMasks()
{ 
   return masks;
}

bool CameraSource::getMonoAndRgbImageAllFrame( CVD::Image<CVD::byte> &mono, CVD::Image<CVD::Rgb<CVD::byte> > &rgb )
{

  if( mbIsConnected )
  {
    if(node_handle.ok())
    {
    	   cv_bridge::CvImageConstPtr rgb_image;
           IplImage* imgCV;
           ros::spinOnce();

          if (rgb_image_msg && rgb_image_msg != old_rgb_image_msg)
          {
			 
            old_rgb_image_msg = rgb_image_msg;
            seqNo = rgb_image_msg->header.seq;

   	    try
	    {
  	          rgb_image = cv_bridge::toCvShare(rgb_image_msg, sensor_msgs::image_encodings::MONO8);
              header_ = rgb_image_msg->header;
              CVD::BasicImage<CVD::byte> cvd_rgb_image_b (rgb_image->image.data, CVD::ImageRef(rgb_image->image.cols, rgb_image->image.rows));
              CVD::ImageRef sX = cvd_rgb_image_b.size();
              convert_image (cvd_rgb_image_b, rgb);
              convert_image (cvd_rgb_image_b, mono);
              objNos.clear();
              masks.clear();
              regions.clear();
              objNos.push_back(0);
              regions.push_back(cv::Rect (10, 10, sX.x-10, sX.y-10));
              return true;
	    }
	    catch (cv_bridge::Exception& e)
	    {
	      ROS_ERROR("cv_bridge exception: %s", e.what());
              return false;
	    }

         }
         else
         {
          return false;
         }
    }
  }
    
  return true;
}

const bool& CameraSource::isConnected() const
{
  return mbIsConnected;
}

