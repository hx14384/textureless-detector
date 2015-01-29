cmake_minimum_required(VERSION 2.4.6)
include($ENV{ROS_ROOT}/core/rosbuild/rosbuild.cmake)

# Set the build type.  Options are:
#  Coverage       : w/ debug symbols, w/o optimization, w/ code-coverage
#  Debug          : w/ debug symbols, w/o optimization
#  Release        : w/o debug symbols, w/ optimization
#  RelWithDebInfo : w/ debug symbols, w/ optimization
#  MinSizeRel     : w/o debug symbols, w/ optimization, stripped binaries
#set(ROS_BUILD_TYPE RelWithDebInfo)

rosbuild_init()
rosbuild_add_boost_directories ()

find_package(OpenCV REQUIRED)

#set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
#set the default path for built libraries to the "lib" directory
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

#set (CMAKE_CXX_FLAGS "-O3 -DNDEBUG -msse -msse2 -msse3 -finline-functions -Wno-write-strings -march=nocona -fno-strict-aliasing" )

#uncomment if you have defined messages
rosbuild_genmsg()
#uncomment if you have defined services
#rosbuild_gensrv()

find_package(Qt4 REQUIRED)
find_package(OpenGL REQUIRED)
SET(QT_USE_QTXML TRUE)
SET(QT_USE_QTOPENGL TRUE)
include(${QT_USE_FILE})
add_definitions(-DQT_NO_KEYWORDS)

include_directories(${CMAKE_CURRENT_BINARY_DIR}) 
file(GLOB QT_FORMS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} src/mainwindow.ui)
file(GLOB QT_MOC RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} src/GLWidget.h src/MainWindow.h src/SaveImageThread.h src/u_mainwindow.h)

QT4_WRAP_UI(QT_FORMS_HPP ${QT_FORMS})
QT4_WRAP_CPP(QT_MOC_HPP ${QT_MOC})

include_directories (. src) 

#common commands for building c++ executables and libraries
#rosbuild_add_library(${PROJECT_NAME} src/example.cpp)
#target_link_libraries(${PROJECT_NAME} another_library)
rosbuild_add_boost_directories()
#rosbuild_link_boost(${PROJECT_NAME} thread)
#rosbuild_add_executable(example examples/example.cpp)
#target_link_libraries(example ${PROJECT_NAME})

file(GLOB QT_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} src/CameraSource.cpp src/CannyEdgeDetector.cpp src/EdgeDetector.cpp src/EdgeObjectDetection.cpp src/GLWidget.cpp src/Homography.cpp src/MainWindow.cpp src/SaveImageThread.cpp src/lsd.c src/EdgeExtraction.cpp)

rosbuild_add_library (mod_nodelets src/MultiObjDetector_nodelet.cpp src/MODetector.cpp ${QT_SOURCES} ${QT_FORMS_HPP} ${QT_MOC_HPP})
target_link_libraries (mod_nodelets ${OUTLIBS} GLEW GLU lapack blas ${QT_LIBRARIES} ${OpenCV_LIBS} ${OPENGL_LIBRARIES})

#rosbuild_add_library (mod_nodelets src/PoseReader_nodelet.cpp)
