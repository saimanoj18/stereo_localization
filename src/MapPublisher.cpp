#include "MapPublisher.h"


void MapPublisher::PublishPose(Matrix4d pose, int type)
{

//    static tf::TransformBroadcaster mTfBr;
//    tf::Transform tfT;
//    tfT.setIdentity();
//    mTfBr.sendTransform(tf::StampedTransform(tfT,ros::Time::now(), "/CamLoc/World", "/CamLoc/Camera"));

    visualization_msgs::Marker mCurrentCamera;
    vector<Matrix4d>* Poses;
    switch (type){
    case 1:
        mCurrentCamera = mOdoCamera;
        Poses = OdoPose;
        break;
    case 2:
        mCurrentCamera = mORBCamera;
        Poses = OrbPose;
        break;
    case 3:
        mCurrentCamera = mEstCamera;
        Poses = EstPose;
        break;
    }

    Poses->push_back(pose);
    geometry_msgs::Point msgs_o_prev;
//    cout<<"poses size: "<<Poses->size()<<endl;
    mCurrentCamera.points.clear();
    float d = fCameraSize;
    //Camera is a pyramid. Define in camera coordinate system
    cv::Mat o = (cv::Mat_<double>(4,1) << 0, 0, 0, 1);
    cv::Mat p1 = (cv::Mat_<double>(4,1) << d, d*0.8, d*0.5, 1);
    cv::Mat p2 = (cv::Mat_<double>(4,1) << d, -d*0.8, d*0.5, 1);
    cv::Mat p3 = (cv::Mat_<double>(4,1) << -d, -d*0.8, d*0.5, 1);
    cv::Mat p4 = (cv::Mat_<double>(4,1) << -d, d*0.8, d*0.5, 1);

    for(vector<Matrix4d>::iterator vit=Poses->begin(), vend=Poses->end(); vit!=vend; vit++){
        cv::Mat Twc;
        cv::eigen2cv((*vit),Twc);
        cv::Mat ow = Twc*o;
        cv::Mat p1w = Twc*p1;
        cv::Mat p2w = Twc*p2;
        cv::Mat p3w = Twc*p3;
        cv::Mat p4w = Twc*p4;

        geometry_msgs::Point msgs_o,msgs_p1, msgs_p2, msgs_p3, msgs_p4;
        msgs_o.x=ow.at<double>(0);
        msgs_o.y=ow.at<double>(1);
        msgs_o.z=ow.at<double>(2);
        msgs_p1.x=p1w.at<double>(0);
        msgs_p1.y=p1w.at<double>(1);
        msgs_p1.z=p1w.at<double>(2);
        msgs_p2.x=p2w.at<double>(0);
        msgs_p2.y=p2w.at<double>(1);
        msgs_p2.z=p2w.at<double>(2);
        msgs_p3.x=p3w.at<double>(0);
        msgs_p3.y=p3w.at<double>(1);
        msgs_p3.z=p3w.at<double>(2);
        msgs_p4.x=p4w.at<double>(0);
        msgs_p4.y=p4w.at<double>(1);
        msgs_p4.z=p4w.at<double>(2);

        mCurrentCamera.points.push_back(msgs_o);
        mCurrentCamera.points.push_back(msgs_p1);
        mCurrentCamera.points.push_back(msgs_o);
        mCurrentCamera.points.push_back(msgs_p2);
        mCurrentCamera.points.push_back(msgs_o);
        mCurrentCamera.points.push_back(msgs_p3);
        mCurrentCamera.points.push_back(msgs_o);
        mCurrentCamera.points.push_back(msgs_p4);
        mCurrentCamera.points.push_back(msgs_p1);
        mCurrentCamera.points.push_back(msgs_p2);
        mCurrentCamera.points.push_back(msgs_p2);
        mCurrentCamera.points.push_back(msgs_p3);
        mCurrentCamera.points.push_back(msgs_p3);
        mCurrentCamera.points.push_back(msgs_p4);
        mCurrentCamera.points.push_back(msgs_p4);
        mCurrentCamera.points.push_back(msgs_p1);
        
        if(vit!=Poses->begin()){
        mCurrentCamera.points.push_back(msgs_o);
        mCurrentCamera.points.push_back(msgs_o_prev);
        }
        msgs_o_prev = msgs_o;
    }

    mCurrentCamera.header.stamp = ros::Time::now();
    publisher.publish(mCurrentCamera);
}

void MapPublisher::PublishMap(const pcl::PointCloud<pcl::PointXYZ>::Ptr& plot_cloud, int type)
{
    visualization_msgs::Marker mPoints;

    switch (type){
    case 1:
        mPoints = mVeloPoints;
        break;
    case 2:
        mPoints = mCamPoints;
        break;
    case 3:
        mPoints = mEstPoints;
        break;
    }

    mPoints.points.clear();
    for (size_t i = 0; i < plot_cloud->points.size(); ++i)
    {
        geometry_msgs::Point p;
        p.x=plot_cloud->points[i].x;
        p.y=plot_cloud->points[i].y;
        p.z=plot_cloud->points[i].z;
        mPoints.points.push_back(p);
    }
    mPoints.header.stamp = ros::Time::now();
    publisher.publish(mPoints);
}
