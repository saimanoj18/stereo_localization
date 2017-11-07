cmake_minimum_required(VERSION 2.8.3)
project(stereo_localization)

#list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(catkin REQUIRED COMPONENTS 
roscpp rospy std_msgs sensor_msgs
cv_bridge
pcl_conversions
pcl_msgs
pcl_ros
tf
tf_conversions)

find_package(OpenCV REQUIRED)
find_package(Eigen3 REQUIRED)
#find_package(G2O REQUIRED)
#find_package(PCL 1.7 REQUIRED)

## Declare a catkin package
catkin_package(DEPENDS OpenCV Eigen3
               INCLUDE_DIRS include
               LIBRARIES stereo_localization
               CATKIN_DEPENDS roscpp  
               std_msgs sensor_msgs geometry_msgs
               cv_bridge
               pcl_ros pcl_conversions pcl_msgs tf tf_conversions)

# Also add some useful compiler flag
set(CMAKE_CXX_FLAGS
   "${CMAKE_CXX_FLAGS} ${SSE_FLAGS} -march=native  -std=c++0x"
) 

###########
## Build ##
###########

#add_library(stereo_localization
#    src/stereo_localization.cpp
#)
set(CMAKE_BUILD_TYPE Debug)

include_directories(
    include
    ${catkin_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
    ${PCL_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/Thirdparty/Sophus
#    SYSTEM ${G2O_INCLUDE_DIR} ${G2O_INCLUDE_DIRS}
)
#include_directories(SYSTEM BEFORE ${G2O_INCLUDE_DIR} ${G2O_INCLUDE_DIRS})

set (SRC
  src/main.cpp
  src/CamLocalization.cpp
  src/MapPublisher.cpp
  src/DataStructures/Frame.cpp
  src/DataStructures/FramePoseStruct.cpp
  src/DataStructures/FrameMemory.cpp
  src/DepthEstimation/DepthMap.cpp
  src/DepthEstimation/DepthMapPixelHypothesis.cpp
  src/util/globalFuncs.cpp
  src/util/SophusUtil.cpp
  src/util/settings.cpp
  src/ImageDisplay_OpenCV.cpp
  #src/GlobalMapping/FabMap.cpp
  #src/GlobalMapping/KeyFrameGraph.cpp
  #src/GlobalMapping/g2oTypeSim3Sophus.cpp
  #src/GlobalMapping/TrackableKeyFrameSearch.cpp
  src/Tracking/Relocalizer.cpp
  src/Tracking/SE3Tracker.cpp
  src/Tracking/Sim3Tracker.cpp
  src/Tracking/TrackingReference.cpp
)

#link_directories(${G2O_LIBRARY_DIRS})

add_executable(stereo_localization ${SRC})
target_link_libraries(stereo_localization 
                      ${catkin_LIBRARIES} 
                      ${OpenCV_LIBS}
#                      ${G2O_LIBRARIES}
                      ${PROJECT_SOURCE_DIR}/Thirdparty/g2o/lib/libg2o.so)

add_dependencies(stereo_localization ${catkin_EXPORTED_TARGETS})
