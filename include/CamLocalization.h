#ifndef CAMLOCALIZATION_H
#define CAMLOCALIZATION_H

#include <iostream>
#include <sys/time.h>
#include <sys/select.h>
#include <numeric>
#include <cmath> 
#include <Eigen/Dense>
#include <Eigen/Geometry> 
#include <Eigen/StdVector>
#include <ros/ros.h>

#include <sensor_msgs/CameraInfo.h>
#include <sensor_msgs/Image.h>
#include <image_transport/image_transport.h>
#include <sensor_msgs/image_encodings.h>
#include <sensor_msgs/PointCloud2.h>
#include <tf/transform_listener.h>
#include <tf_conversions/tf_eigen.h>


#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/common/centroid.h>
#include <pcl/io/ply_io.h>
#include <pcl/filters/voxel_grid.h>

#include <pcl/octree/octree.h>
#include <pcl/octree/octree_pointcloud.h>
#include <pcl/octree/octree_iterator.h>
#include <pcl/octree/octree_impl.h> 

#include <pcl/correspondence.h>
#include <pcl/common/transforms.h>
#include <pcl_ros/transforms.h>

#include <opencv2/highgui/highgui.hpp>
#include <cv_bridge/cv_bridge.h> 
#include <opencv2/calib3d.hpp>

#include "Thirdparty/g2o/g2o/core/block_solver.h"
#include "Thirdparty/g2o/g2o/core/optimization_algorithm_levenberg.h"
#include "Thirdparty/g2o/g2o/solvers/linear_solver_eigen.h"
#include "Thirdparty/g2o/g2o/types/types_six_dof_expmap.h"
#include "Thirdparty/g2o/g2o/core/robust_kernel_impl.h"
#include "Thirdparty/g2o/g2o/solvers/linear_solver_dense.h"
#include "Thirdparty/g2o/g2o/types/types_seven_dof_expmap.h"

#include <liblas/liblas.hpp>
#include <liblas/reader.hpp>

#include "MapPublisher.h"

using namespace std;
using namespace Eigen;

class CamLocalization
{
public:
    CamLocalization():
    velo_raw(new pcl::PointCloud<pcl::PointXYZ>),velo_cloud(new pcl::PointCloud<pcl::PointXYZ>),velo_xyzi(new pcl::PointCloud<pcl::PointXYZI>),velo_global(new pcl::PointCloud<pcl::PointXYZ>),
    fakeTimeStamp(0),frameID(0),
    mode(0),scale(0.42553191),//scale(0.7),
//    mode(1),scale(0.382),//0.472
    Velo_received(false),Left_received(false),Right_received(false), octree(128.0f)
    {
        it = new image_transport::ImageTransport(nh);
        
        //Set Subscriber
        sub_veloptcloud = nh.subscribe("/kitti/velodyne_points", 1, &CamLocalization::VeloPtsCallback, this);
        sub_leftimg = it->subscribeCamera("/kitti/left_image", 10,&CamLocalization::LeftImgCallback, this);
        sub_rightimg = it->subscribeCamera("/kitti/right_image", 10,&CamLocalization::RightImgCallback, this);         
        
//        sub_leftimg = nh.subscribe("/kitti/left_image", 1, &CamLocalization::LeftImgCallback, this);
//        sub_rightimg = nh.subscribe("/kitti/right_image", 1, &CamLocalization::RightImgCallback, this);
//        sub_caminfo = nh.subscribe("/kitti/camera_gray_left/camera_info", 1, &CamLocalization::CamInfoCallback, this);


        EST_pose = Matrix4d::Identity();
        IN_pose = Matrix4d::Identity();
        update_pose = Matrix4d::Identity();
        if(mode ==0 )update_pose(2,3) = 0.8;
        else update_pose(2,3) = 0.3;
        optimized_T = Matrix4d::Identity();
        GT_pose = Matrix4d::Identity();

        base_line = 0.482;

//        data_path_ = "/media/youngji/storagedevice/naver_data/20180125_kitti";
        data_path_ = "/home/irap/data/20180125_kitti";
//        data_path_ = "/home/irap/data/20171120_kitti";

//        read_poses("poses.txt");
//        cout<<"Pose loading is completed"<<endl;
//        cout.precision(20);

    }
    ~CamLocalization(){
        //references.clear();
        delete [] ref_container;
        delete [] igx_container;
        delete [] igy_container;
    
    }
    void CamLocInitialize(cv::Mat image);
    void Refresh();    
    void write_times(std::string fname, float time_diff); 
    
private:



    //for ros subscription
    ros::NodeHandle nh;
    image_transport::ImageTransport *it;
    ros::Subscriber sub_veloptcloud;
    image_transport::CameraSubscriber sub_leftimg;
    image_transport::CameraSubscriber sub_rightimg;
//    ros::Subscriber sub_caminfo;
    tf::TransformListener tlistener;

    //for broadcast    
    tf::TransformBroadcaster mTfBr;

    //for time computation
    int64_t start_time;
    int64_t end_time;    

    //input data
    pcl::PointCloud<pcl::PointXYZ>::Ptr velo_cloud;
    pcl::PointCloud<pcl::PointXYZ>::Ptr velo_raw;
    pcl::PointCloud<pcl::PointXYZ>::Ptr velo_global;
    pcl::PointCloud<pcl::PointXYZI>::Ptr velo_xyzi;
    pcl::octree::OctreePointCloudSearch<pcl::PointXYZ> octree;
    cv::Mat left_image;
    cv::Mat right_image;
    cv::Mat ref_image;
    cv::Mat ref_igx;
    cv::Mat ref_igy;
    cv::Mat ref_depth;
    cv::Mat ref_depth_info;
    float* ref_container;
    float* igx_container;
    float* igy_container;
    double fakeTimeStamp;
    int frameID;
    std::string data_path_;
    cv::Mat igx_image, igy_image;
    cv::Mat left_scaled;
    cv::Mat dgx_image, dgy_image; 
    cv::Mat disp;
    cv::Mat cur_depth_info;
    float* depth;// = new float[width*height]();
    float* image_info;// = new float[width*height]();
    float* depth_gradientX;// = new float[width*height]();
    float* depth_gradientY;// = new float[width*height]();
    float* depth_info;// = new float[width*height]();

    //input transform    
    tf::StampedTransform ctv;
    tf::StampedTransform wtb;
    tf::StampedTransform tfT;
    Matrix4d cTv;

    //input camera info
    Matrix<double,3,4> P0;
    Matrix<double,3,4> P1;
    Matrix3f K;
    double base_line;
    int width;
    int height;
    int ancient_width;
    double scale;
    float matching_thres;
    float d_var;
    float d_limit;

    //result data
    ros::Time ODO_time;
    Matrix4d IN_pose;
    Matrix4d EST_pose;
    Matrix4d GT_pose;
    Matrix4d update_pose;
    Matrix4d optimized_T;
    vector<Matrix4d, Eigen::aligned_allocator<Eigen::Vector4f>> GT_poses;
 

    //Callbacks
    void VeloPtsCallback(const sensor_msgs::PointCloud2::ConstPtr& msg);
//    void LeftImgCallback(const sensor_msgs::Image::ConstPtr& msg);
//    void RightImgCallback(const sensor_msgs::Image::ConstPtr& msg);
    void LeftImgCallback(const sensor_msgs::ImageConstPtr& msg, const sensor_msgs::CameraInfoConstPtr & infomsg);
    void RightImgCallback(const sensor_msgs::ImageConstPtr& msg, const sensor_msgs::CameraInfoConstPtr & infomsg);
    void CamInfoCallback(const sensor_msgs::CameraInfo::ConstPtr& msg);
    bool Velo_received; 
    bool Left_received; 
    bool Right_received;
    int8_t mode;
    void read_poses(std::string fname); 
    void write_poses(std::string fname, Matrix4d saved_pose); 

    //main algorithms
    Matrix4d visual_tracking(const float* ref, const float* r_igx, const float* r_igy, const float* i_var, const float* idepth, cv::Mat cur,Matrix4d init_pose, float thres);
    Matrix4d Optimization(const float* idepth, const float* idepth_var, const float* d_gradientX, const float* d_gradientY, float thres); 
    void depth_propagation(float* idepth, cv::Mat info, Matrix4d pose);    

    void debugImage(cv::Mat& depth_image,cv::Mat& dgx_image,cv::Mat& dgy_image,const float* depth_info);

    int64_t
    timestamp_now (void)
    {
        struct timeval tv;
        gettimeofday (&tv, NULL);
        return (int64_t) tv.tv_sec * 1000000 + tv.tv_usec;
    }

    Vector2d ProjectTo2D(Vector3d v)
    {
      Vector2d res;
      res[0] = v[0]*K(0,0)/v[2] + K(0,2);
      res[1] = v[1]*K(1,1)/v[2] + K(1,2);
      return res;
    }

    Vector3d ReprojectTo3D(double v1, double v2, double v3)
    {
      Vector3d res;
      res[0] = v3/K(0,0)*(v1-K(0,2));
      res[1] = v3/K(1,1)*(v2-K(1,2)); 
      res[2] = v3;
      return res;
    }

    Matrix4d SE3toMat(const g2o::SE3Quat &SE3)
    {
    Eigen::Matrix3d eigR = SE3.rotation().toRotationMatrix();//.cast <float> ();
    Eigen::Vector3d eigt = SE3.translation();//.cast <float> ();
    Matrix4d T = Matrix4d::Identity();
    T.block<3,3>(0,0) = eigR;
    T(0,3) = eigt[0];
    T(1,3) = eigt[1];
    T(2,3) = eigt[2];
    return T;
    }

    Matrix4d Sim3toMat(const g2o::Sim3 &Sim3)
    {
    Eigen::Matrix3d eigR = Sim3.rotation().toRotationMatrix();//.cast <float> ();
    Eigen::Vector3d eigt = Sim3.translation();//.cast <float> ();
    float s = (float) Sim3.scale();
    Matrix4d T = Matrix4d::Identity();
    T.block<3,3>(0,0) = s*eigR;
    T(0,3) = eigt[0];
    T(1,3) = eigt[1];
    T(2,3) = eigt[2];
    return T;
    }
    
    //Degug images
    cv::Vec3b Compute_error_color(float depth_error, float range)
    {
        // rainbow between 0 and 4
	    float r = (0-depth_error) * 255 / range; if(r < 0) r = -r;
	    float g = (range/2-depth_error) * 255 / range; if(g < 0) g = -g;
	    float b = (range-depth_error) * 255 / range; if(b < 0) b = -b;
	    uchar rc = r < 0 ? 0 : (r > 255 ? 255 : r);
	    uchar gc = g < 0 ? 0 : (g > 255 ? 255 : g);
	    uchar bc = b < 0 ? 0 : (b > 255 ? 255 : b);

        
        return cv::Vec3b(255-rc,255-gc,255-bc);
    }

    void save_colormap(cv::Mat& img, std::string img_name, double min, double max)
    {
//        double min;
//        double max;
//        cv::minMaxIdx(img, &min, &max);
        cv::Mat adjMap;
        img.convertTo(adjMap,CV_8UC1, 255 / (max-min), -min); 
        cv::Mat falseColorsMap;
        applyColorMap(adjMap, falseColorsMap, cv::COLORMAP_JET);
        imwrite(img_name, falseColorsMap);
    }

    void save_colormap2d(cv::Mat& img1, cv::Mat& img2, std::string img_name)
    {
        double min;
        double max;

        cv::Mat adjMap1;
        img1.convertTo(adjMap1,CV_32FC1, 90.0 / M_PI, 90.0);


        cv::Mat adjMap2;
        cv::minMaxIdx(img2, &min, &max);
        img2.convertTo(adjMap2,CV_32FC1, 255.0 / (max), 0.0);

        cv::Mat ColorsMap = cv::Mat::zeros(cv::Size(img1.cols, img1.rows), CV_8UC3);
        for(size_t u=0; u<width;u++)
        {
            for(size_t v=0; v<height;v++)
            {
                ColorsMap.at<cv::Vec3b>(v,u)= cv::Vec3b((int)adjMap1.at<float>(v,u),(int)adjMap2.at<float>(v,u),255);
                //if(adjMap1.at<float>(v,u)>127)cout<<adjMap1.at<float>(v,u)<<", "<<adjMap2.at<float>(v,u)<<endl;
            }
        }
        cv::Mat ColorsMap_RGB;
        cvtColor(ColorsMap, ColorsMap_RGB, CV_HSV2RGB);
        imwrite(img_name, ColorsMap_RGB);
    }

    //MapPublisher
    MapPublisher MapPub;
    


};


#endif // CAMLOCALIZATION_H
