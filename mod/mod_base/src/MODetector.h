#ifndef MODETECTOR_H
#define MODETECTOR_H

#include <QThread>
#include <ros/ros.h>

namespace mod_base {
  class MODetector// : public QThread
  {
     public:
       MODetector ();
       void Start (ros::Publisher detectionImage_pub_, ros::Publisher detected_objs_, std::string colorImgSrc, std::string codebookFilename, std::string loadCodebookFilename, std::string outputFilename, int imgSize, int learnObjNo, bool withInterface, bool isStar, bool andPublish, bool andWriteToFile, int scenarioNo);
     protected:
       void run (ros::Publisher detectionImage_pub_, std::string colorImgSrc, int imgSize, int learnObjNo);
       void runDetector(ros::Publisher detectionImage_pub_, std::string colorImgSrc, int imgSize, int learnObjNo, int argc, char **argv);
     private:
       bool run_yet;
       bool withInterface;
       ros::Publisher detected_objs_;
       bool andPublish;
       bool andWriteToFile;
       bool isStar;
       std::string codebookFilename;
       std::string outputFilename;
       std::string loadCodebookFilename;
       int scenarioNo;
  };
}

#endif
