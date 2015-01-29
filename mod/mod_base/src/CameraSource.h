/*
  Computer Vision Group, Department of Computer Science. University of      
  Bristol. This code can not be used or copied without express permission.   
*/

// Modefied from btl/VideoSource
#ifndef CAMERASOURCE_H_
#define CAMERASOURCE_H_

#include <string>
#include <cvd/videosource.h>
#include <cvd/Linux/v4lcontrol.h>
#include <ros/ros.h>
#include <boost/thread/mutex.hpp>
#include <sensor_msgs/Image.h>
#include <message_filters/subscriber.h>
#include <message_filters/synchronizer.h>
#include <message_filters/sync_policies/approximate_time.h>
#include <cv_bridge/cv_bridge.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>


/**
 @file CameraSource.h
 @brief CameraSource is a class for handling a hardware camera or image sequences.
*/
class CameraSource
{
  public:
    CameraSource();
    ~CameraSource();

    /// Connects a camera using the CVD videosource format
    void connect( const std::string& source, int argc = 0, char *argv[] = NULL);
    void connect()
    {
      std::string source( "dc1394:///0" );
      connect( source );
    }
    ;
    void disconnect();
    void update();

    CVD::ImageRef getSize()
    {
      return mirImageSize;
    }

    int getSeqNo() { return seqNo; }

    void setImgSize (CVD::ImageRef r)
    {
       mirImageSize = r;
    }
    /**
      @param[out] mono is an output image in Gray format.
      @param[out] rgb is an output image in RGB format.
    */
    bool getMonoAndRgbImageAllFrame( CVD::Image<CVD::byte> &mono, CVD::Image<CVD::Rgb<CVD::byte> > &rgb );

    std::vector<cv::Rect> getRegions();
    std::vector<sensor_msgs::Image> getMasks();
    std::vector<int> getObjNos () {return objNos;};

    CVD::Image<CVD::byte> getDepthEdges();

    void rgb_image_msg_cb2 (const sensor_msgs::ImageConstPtr& msg);

    std_msgs::Header getHeader () {return header_;}

    void setSources (std::string colorImg) {colorImgSrc = colorImg;}

    const bool& isConnected() const;

    inline bool isImageSequence()
    {
      return mbIsImageSequence;
    }
    ;
  private:

    bool mbIsConnected;

    CVD::ImageRef mirImageSize;
    sensor_msgs::ImageConstPtr rgb_image_msg, old_rgb_image_msg;
    ros::NodeHandle node_handle;
    ros::Subscriber rgb_image_subscriber;

    typedef message_filters::Subscriber<sensor_msgs::Image> ImageSubscriber;
    boost::shared_ptr<ImageSubscriber> image_sub_;

    message_filters::Subscriber<sensor_msgs::Image> imgSub;

    // internal mutex
    boost::mutex mutex_; 
 
    std_msgs::Header header_;

    std::vector<cv::Rect> regions;
    std::vector<sensor_msgs::Image> masks;
    std::vector<int> objNos;
    bool CAMERA_READ;

    boost::mutex m;
    std::string colorImgSrc;

    bool mbIsImageSequence;

    int seqNo;
};

#endif /*CAMERASOURCE_H_*/
