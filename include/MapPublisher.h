#ifndef MAPPUBLISHER_H
#define MAPPUBLISHER_H

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry> 
#include <ros/ros.h>

#include <visualization_msgs/Marker.h>
#include <tf/transform_broadcaster.h>
#include <tf/transform_datatypes.h>
#include <tf_conversions/tf_eigen.h>

#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/common/centroid.h>
#include <pcl/correspondence.h>
#include <pcl/common/transforms.h>
#include <pcl_ros/transforms.h> 
#include <opencv2/core/eigen.hpp>

using namespace std;
using namespace Eigen;

class MapPublisher
{
public:
    MapPublisher()
    {
        const char* MAP_FRAME_ID = "/CamLoc/World";
        const char* VELO_POINTS_NAMESPACE = "VeloPoints";
        const char* CAMERA_POINTS_NAMESPACE = "CameraPoints";
        const char* EST_POINTS_NAMESPACE = "EstPoints";
        const char* CAMERA1_NAMESPACE = "OdoCamera";
        const char* CAMERA2_NAMESPACE = "ORBCamera";
        const char* CAMERA3_NAMESPACE = "EstCamera";
        const char* GRAPH1_NAMESPACE = "OdoGraph";
        const char* GRAPH2_NAMESPACE = "ORBGraph";
        const char* GRAPH3_NAMESPACE = "EstGraph";


        //Configure MapPoints
        fPointSize=0.1;
        mVeloPoints.header.frame_id = MAP_FRAME_ID;
        mVeloPoints.ns = VELO_POINTS_NAMESPACE;
        mVeloPoints.id=0;
        mVeloPoints.type = visualization_msgs::Marker::POINTS;
        mVeloPoints.scale.x=fPointSize;
        mVeloPoints.scale.y=fPointSize;
        mVeloPoints.pose.orientation.w=1.0;
        mVeloPoints.action=visualization_msgs::Marker::ADD;
        mVeloPoints.color.r = 0.5f;
        mVeloPoints.color.g = 0.5f;
        mVeloPoints.color.b = 0.5f;
        mVeloPoints.color.a = 1.0;

        fPointSize=0.05;
        mCamPoints.header.frame_id = MAP_FRAME_ID;
        mCamPoints.ns = CAMERA_POINTS_NAMESPACE;
        mCamPoints.id=1;
        mCamPoints.type = visualization_msgs::Marker::POINTS;
        mCamPoints.scale.x=fPointSize;
        mCamPoints.scale.y=fPointSize;
        mCamPoints.pose.orientation.w=1.0;
        mCamPoints.action=visualization_msgs::Marker::ADD;
        mCamPoints.color.g = 1.0f;
        mCamPoints.color.a = 1.0;

        fPointSize=0.05;
        mEstPoints.header.frame_id = MAP_FRAME_ID;
        mEstPoints.ns = EST_POINTS_NAMESPACE;
        mEstPoints.id=2;
        mEstPoints.type = visualization_msgs::Marker::POINTS;
        mEstPoints.scale.x=fPointSize;
        mEstPoints.scale.y=fPointSize;
        mEstPoints.pose.orientation.w=1.0;
        mEstPoints.action=visualization_msgs::Marker::ADD;
        mEstPoints.color.b = 1.0f;
        mEstPoints.color.a = 1.0;

        fCameraSize=0.1;
        //Configure ODO Camera
        mOdoCamera.header.frame_id = MAP_FRAME_ID;
        mOdoCamera.ns = CAMERA1_NAMESPACE;
        mOdoCamera.id=3;
        mOdoCamera.type = visualization_msgs::Marker::LINE_LIST;
        mOdoCamera.scale.x= 0.04;//0.2; 0.03
        mOdoCamera.pose.orientation.w=1.0;
        mOdoCamera.action=visualization_msgs::Marker::ADD;
        mOdoCamera.color.r= 1.0f;
//        mOdoCamera.color.g= 0.5f;
//        mOdoCamera.color.b= 0.5f;
        mOdoCamera.color.a = 1.0;

        //Configure ORB Camera
        mORBCamera.header.frame_id = MAP_FRAME_ID;
        mORBCamera.ns = CAMERA2_NAMESPACE;
        mORBCamera.id=4;
        mORBCamera.type = visualization_msgs::Marker::LINE_LIST;
        mORBCamera.scale.x=0.04;//0.2; 0.03
        mORBCamera.pose.orientation.w=1.0;
        mORBCamera.action=visualization_msgs::Marker::ADD;
        mORBCamera.color.g=1.0f;
        mORBCamera.color.a = 1.0;

        //Configure Est Camera
        mEstCamera.header.frame_id = MAP_FRAME_ID;
        mEstCamera.ns = CAMERA3_NAMESPACE;
        mEstCamera.id=5;
        mEstCamera.type = visualization_msgs::Marker::LINE_LIST;
        mEstCamera.scale.x=0.04;//0.2; 0.03
        mEstCamera.pose.orientation.w=1.0;
        mEstCamera.action=visualization_msgs::Marker::ADD;
        mEstCamera.color.b=1.0f;
        mEstCamera.color.a = 1.0;

        //Configure Publisher
        publisher = nh.advertise<visualization_msgs::Marker>("CamLoc/Map", 10);

        OdoPose = new vector<Matrix4d>();
        OrbPose = new vector<Matrix4d>();
        EstPose = new vector<Matrix4d>();

    }

    //publish
    void PublishPose(Matrix4d pose, int type);
    void PublishMap(const pcl::PointCloud<pcl::PointXYZ>::Ptr& plot_cloud, int type); 

private:
    ros::NodeHandle nh;
    ros::Publisher publisher;

    visualization_msgs::Marker mVeloPoints;
    visualization_msgs::Marker mCamPoints;
    visualization_msgs::Marker mEstPoints;
    visualization_msgs::Marker mOdoCamera;
    visualization_msgs::Marker mORBCamera;
    visualization_msgs::Marker mEstCamera;
    visualization_msgs::Marker mOdoGraph;
    visualization_msgs::Marker mORBGraph;
    visualization_msgs::Marker mEstGraph;

    vector<Matrix4d>* OdoPose;
    vector<Matrix4d>* OrbPose;
    vector<Matrix4d>* EstPose;

    float fCameraSize;
    float fPointSize;


};


#endif // MAPPUBLISHER_H
