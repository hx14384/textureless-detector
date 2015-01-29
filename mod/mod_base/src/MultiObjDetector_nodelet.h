#ifndef MULTIOBJDETECTOR_NODELET_H
#define MULTIOBJDETECTOR_NODELET_H

#include <ros/ros.h>
#include <nodelet/nodelet.h>
#include <sensor_msgs/Image.h>
#include <boost/thread/mutex.hpp>
#include <sensor_msgs/Image.h>
#include <message_filters/subscriber.h>


#include "MODetector.h"


namespace mod_base {

/**
 * This class implements the interface for a background segmenter nodelet.
 */
class MultiObjDetectorNodelet : public nodelet::Nodelet
{
public:
	MultiObjDetectorNodelet() { loadCodebookFilename = ""; isKinect = true; andPublish = true;  started=false;};

	virtual void onInit()
	{
		ros::NodeHandle& private_nh = getMTPrivateNodeHandle();
		childInit(private_nh);
		initConfigureService(private_nh);
                private_nh.getParam ("CodebookFilename", codebookFilename);
                private_nh.getParam ("LoadCodebookFilename", loadCodebookFilename);

                private_nh.getParam ("rgb_image", colorImgSrc);
          
                std::string scenarioStr;
                scenarioNo = 0;

                std::string imageSizeStr;
                private_nh.getParam ("ImageSize", imageSizeStr);
                if (imageSizeStr == "VGA") 
                   kinectImgSize = 2;
                else 
                   kinectImgSize = 1;
                std::string pubStr;
                private_nh.getParam ("PublishResults", pubStr);
                if (pubStr == "PUBLISH")
                    andPublish = true;
                else
                    andPublish = false;
                std::string starStr;
                private_nh.getParam ("PathShape", pubStr);
                if (pubStr == "STAR")
                    isStar = true;
                else
                    isStar = false;
                std::string interfaceStr;
                private_nh.getParam ("Interface", interfaceStr);
                if (interfaceStr == "SHOW")
                    withInterface = true;
                else
                    withInterface = false;
                private_nh.getParam ("WriteToFile", outputFilename);
                if (outputFilename == "")
                    andWriteToFile = false;
                else
                    andWriteToFile = true;
                std::string objNoStr;
                private_nh.getParam ("LearnObjNo", objNoStr);
                learnObjNo = atoi (objNoStr.substr(3,4).c_str());
		subscribeTopics(private_nh);
		advertiseTopics(private_nh);
	}

        virtual void childInit(ros::NodeHandle& nh);

        void rgb_image_msg_cb (const sensor_msgs::ImageConstPtr& msg);

protected:
	/**
	 * Subscribe to the topics of interest. Depending on
	 * some parameters (e.g. use_rois) it can subscribe to different topics.
	 */
	virtual void subscribeTopics(ros::NodeHandle& nh);

	/**
	 *
	 * @param nh
	 */
	virtual void advertiseTopics(ros::NodeHandle& nh);

	virtual void initConfigureService(ros::NodeHandle& nh) {}

	/**
	 * Runs the pose estimation.
	 */
	virtual void runNodelet();

	/**
	 * Publishes detection results
	 */

protected:
	std_msgs::Header header_;

private:
	// publishers
        ros::Publisher detected_objs_;
        ros::Publisher detectionImage_pub_;
        message_filters::Subscriber<sensor_msgs::Image> imgSub;
        typedef message_filters::Subscriber<sensor_msgs::Image> ImageSubscriber;
        boost::shared_ptr<ImageSubscriber> image_sub_;
        bool started;


        int counter;
        std::string codebookFilename;
        std::string loadCodebookFilename;
        std::string outputFilename;

        std::string colorImgSrc;

        int kinectImgSize;
        bool isDepth;
        bool isLearn;
        bool isMask;
        bool isKinect;
        bool isGaze;
        bool isFireWire;
        bool withBG;
        int learnObjNo;
        bool autoStart;
        bool withInterface;
        bool asService;
        bool andPublish;
        bool andWriteToFile;
        bool andPrintToFile;
        int scenarioNo;
        bool isStar;

        MODetector modetector;
};

}

#endif
