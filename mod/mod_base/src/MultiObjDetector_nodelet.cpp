#include "MultiObjDetector_nodelet.h"

#include <pluginlib/class_list_macros.h>

#include <mod_base/DetectedObjectArray.h>

namespace mod_base
{

void MultiObjDetectorNodelet::childInit(ros::NodeHandle& nh)
{
   
}

void MultiObjDetectorNodelet::rgb_image_msg_cb (const sensor_msgs::ImageConstPtr& msg)
{
   if (!started)
   {
       runNodelet();
       started = true;
   }
}

void MultiObjDetectorNodelet::subscribeTopics(ros::NodeHandle& nh)
{
    imgSub.subscribe(nh, colorImgSrc, 40);
    imgSub.registerCallback(&MultiObjDetectorNodelet::rgb_image_msg_cb, this);
}

void MultiObjDetectorNodelet::advertiseTopics(ros::NodeHandle& nh)
{
    detected_objs_ = nh.advertise<mod_base::DetectedObjectArray>("detectedObjects", 100);
    detectionImage_pub_ = nh.advertise<sensor_msgs::Image> ("detectionImage", 20);
}

void MultiObjDetectorNodelet::runNodelet()
{
   // here I need to find a way to run the detector and retrieve information from it
   counter++;
   modetector.Start(detectionImage_pub_, detected_objs_, colorImgSrc, codebookFilename, loadCodebookFilename, outputFilename, kinectImgSize, learnObjNo, withInterface, isStar, andPublish, andWriteToFile, scenarioNo);
}
}

/**
* Pluginlib declaration. This is needed for the nodelet to be dynamically loaded/unloaded
*/
typedef mod_base::MultiObjDetectorNodelet MultiObjDetectorNodelet;
PLUGINLIB_DECLARE_CLASS (mod_base, MultiObjDetectorNodelet, MultiObjDetectorNodelet, nodelet::Nodelet);
