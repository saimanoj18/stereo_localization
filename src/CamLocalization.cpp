#include "CamLocalization.h"

void CamLocalization::CamLocInitialize(cv::Mat image)
{
    //Set matrices
//    float fx = 718.856*image.cols/ancient_width;
//    float fy = 718.856*image.cols/ancient_width;        
//    float cx = 607.1928*image.cols/ancient_width;
//    float cy = 185.2157*image.cols/ancient_width;

    float fx = P0(0,0)*image.cols/ancient_width;
    float fy = P0(1,1)*image.cols/ancient_width;        
    float cx = P0(0,2)*image.cols/ancient_width;
    float cy = P0(1,2)*image.cols/ancient_width;
    base_line = -P1(0,3)/P0(0,0);// because R = I;

    cout<<"base_line: "<<base_line<<endl;
    K << fx, 0.0, cx, 0.0, fy, cy, 0.0, 0.0, 1.0;
    cout<<"K_mat: "<<K<<endl;
    cout<<"extrinsic Camera to Velodyne: "<<cTv<<endl;
    cout<<left_image.cols<<", "<<left_image.rows<<endl;
    width = left_image.cols;
    height = left_image.rows;
    
    //Set ref images
    ref_container = new float[width*height];
    igx_container = new float[width*height];
    igy_container = new float[width*height];

    //set initial pose
    EST_pose = GT_pose;

    //set matching thres
    if(mode == 0){   
        d_var = 1.00;
        d_limit = 100.0;
        matching_thres = K(0,0)*base_line*( 1.0/(100.0/16.0) + d_var/((float)(100.0/16.0)*(100.0/16.0)*(100.0/16.0)) );
    }

    if (mode ==1)
    {


        d_var = 1.0;
        d_limit =100.0;
        matching_thres = K(0,0)*base_line*( 1.0/(d_limit/16.0) + d_var/((float)(d_limit/16.0)*(d_limit/16.0)*(d_limit/16.0)) );

        //load velo_raw from .las
        char filename[1000];
        sprintf(filename, "/media/youngji/storagedevice/naver_data/20180125_kitti/sequences/11/sick_pointcloud.las");
        std::ifstream ifs;
        if (!liblas::Open(ifs, filename))
        {
            throw std::runtime_error(std::string("Can not open ") + filename);
        }

        liblas::Reader reader(ifs);
        liblas::Header const& h = reader.GetHeader();        
        velo_raw->width = h.GetPointRecordsCount();
        velo_raw->height = 1;
        velo_raw->points.resize (velo_raw->width * velo_raw->height);
        int count = 0;
        int iter = 0;
//        float ox ,oy, oz;
        Matrix4f init_trans = Matrix4f::Identity();
        while (reader.ReadNextPoint())
        {        
            liblas::Point const& p = reader.GetPoint();
            if(iter == 0)
            {
                //ifstream file("/media/youngji/storagedevice/Initial_pose.txt");
                ifstream file("/media/youngji/storagedevice/naver_data/20180125_kitti/sequences/11/Initial_pose.txt");
                string s;
                getline(file, s, ' ');
                init_trans(0,0) = atof(s.c_str());
                getline(file, s, ' ');
                init_trans(0,1) = atof(s.c_str());
                getline(file, s, ' ');
                init_trans(0,2) = atof(s.c_str());
                getline(file, s, ' ');
                init_trans(0,3) = atof(s.c_str());
                getline(file, s, ' ');
                init_trans(1,0) = atof(s.c_str());
                getline(file, s, ' ');
                init_trans(1,1) = atof(s.c_str());
                getline(file, s, ' ');
                init_trans(1,2) = atof(s.c_str());
                getline(file, s, ' ');
                init_trans(1,3) = atof(s.c_str());
                getline(file, s, ' ');
                init_trans(2,0) = atof(s.c_str());
                getline(file, s, ' ');
                init_trans(2,1) = atof(s.c_str());
                getline(file, s, ' ');
                init_trans(2,2) = atof(s.c_str());
                getline(file, s, ' ');
                init_trans(2,3) = atof(s.c_str());
                file.close();
            }
            if(iter%2 == 0){
            velo_raw->points[count].x = p[0];
            velo_raw->points[count].y = p[1];
            velo_raw->points[count].z = p[2];
            count++;
            }
            iter++;
        }
        cout<<init_trans<<endl;
        cout<<init_trans.inverse()<<endl;
        pcl::transformPointCloud (*velo_raw, *velo_raw, init_trans.inverse());
        pcl::transformPointCloud (*velo_raw, *velo_raw, cTv);

        //filtering
        pcl::VoxelGrid<pcl::PointXYZ> sor;
        sor.setInputCloud (velo_raw);
        sor.setLeafSize (0.01f, 0.01f, 0.01f);
        sor.filter (*velo_raw);

        //set octree
        octree.setInputCloud (velo_raw);
        octree.addPointsFromInputCloud ();
        
        //extract local map
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
        pcl::PointXYZ searchPoint;
        searchPoint.x = EST_pose(0,3);
        searchPoint.y = EST_pose(1,3);
        searchPoint.z = EST_pose(2,3);
        std::vector<int> pointIdxRadiusSearch;
        std::vector<float> pointRadiusSquaredDistance;
        octree.radiusSearch (searchPoint, 100.0f, pointIdxRadiusSearch, pointRadiusSquaredDistance);

        cloud->width = pointIdxRadiusSearch.size ();
        cloud->height = 1;
        cloud->points.resize (cloud->width * cloud->height);
        count = 0;
        for (size_t i = 0; i < pointIdxRadiusSearch.size (); ++i)
        {
            cloud->points[count].x = velo_raw->points[ pointIdxRadiusSearch[i] ].x;
            cloud->points[count].y = velo_raw->points[ pointIdxRadiusSearch[i] ].y;
            cloud->points[count].z = velo_raw->points[ pointIdxRadiusSearch[i] ].z;
            count++;
        }
        MapPub.PublishMap(cloud,1);

    }        

}

void CamLocalization::Refresh()
{
    bool success;
    success = tlistener.waitForTransform("/kitti/World", "/kitti/Velodyne", ros::Time(0), ros::Duration(0.1));
    if (success) {
        tlistener.lookupTransform("/kitti/World", "/kitti/Velodyne", ros::Time(0), ctv);
        pcl_ros::transformAsMatrix (ctv, cTv); 
    }
    
    if(mode ==1) Velo_received = true;

//    if(mode == 1)MapPub.PublishMap(velo_raw,2);//publish velo raw

    if(Velo_received && Left_received && Right_received)
    {
        start_time = timestamp_now ();

        //prepar GT pose and map point clouds        
        if(mode == 0)pcl::transformPointCloud (*velo_cloud, *velo_raw, GT_pose);//transform to world coordinate
//        mTfBr.sendTransform(tf::StampedTransform(wtb,ros::Time::now(), "/CamLoc/World", "/CamLoc/Camera"));

        //initialize 
        if(frameID == 0)CamLocInitialize(right_image);
        
        /////////////////////////disparity map generation/////////////////////////////        
        cv::Mat disp = cv::Mat::zeros(cv::Size(width, height), CV_16S);
        cv::Ptr<cv::StereoSGBM> sbm;
        if(mode == 0)sbm = cv::StereoSGBM::create(0,16*5,7);
        if(mode == 1)sbm = cv::StereoSGBM::create(0,16*5,7);
        sbm->compute(left_image, right_image, disp);
        frameID = frameID+2;

        /////////////////////////depth image generation/////////////////////////////
        float* depth_raw = new float[width*height]();
        float* depth = new float[width*height]();
        //cv::Mat depth_image = cv::Mat(cv::Size(left_image.cols, left_image.rows), CV_32FC1);
        cv::Mat depth_image = cv::Mat::zeros(cv::Size(left_image.cols, left_image.rows), CV_32FC1);
//        float d_var = 0.01;
        for(size_t i=0; i<width*height;i++)
        {
            int u, v;
            u = i%width;
            v = i/width;
            
            //depth
            int16_t d = disp.at<int16_t>(v,u);            
            if(d==0 || d!=d ) d = 0; //
            depth_raw[i] = K(0,0)*base_line*( 1.0/((float)d/16.0) + d_var/((float)(d/16.0)*(d/16.0)*(d/16.0)) );
            if(d==0 || d!=d || d<d_limit) d = 0; //
            depth[i] = K(0,0)*base_line*( 1.0/((float)d/16.0) + d_var/((float)(d/16.0)*(d/16.0)*(d/16.0)) );//base_line*K(0,0)/disp2.at<float>(v,u);


//            if(depth[i]>30.0)depth[i]= -1.0;

            //depth image            
            depth_image.at<float>(v,u) = depth[i];
            
            //reference images
            if(frameID>2){
                ref_container[i] = ref_image.at<float>(v,u);
                igx_container[i] = ref_igx.at<float>(v,u)/32.0f;
                igy_container[i] = ref_igy.at<float>(v,u)/32.0f;
            }

        }

//        cv::imshow("depth_image", depth_image);
//        cv::waitKey(3);
//        cv::moveWindow("depth_image", 50,20);

        /////////////////////////depth gradient generation/////////////////////////////
        cv::Mat dgx_image, dgy_image; // = cv::Mat(cv::Size(left_image.cols, left_image.rows), CV_8UC3);
        cv::Scharr(depth_image, dgx_image, CV_32FC1, 1, 0);
        cv::Scharr(depth_image, dgy_image, CV_32FC1, 0, 1);
        cv::Mat igx_image, igy_image, left_scaled;
        cv::normalize(left_image, left_scaled, 0, 1, CV_MINMAX, CV_32FC1);   
        cv::Scharr(left_scaled, igx_image, CV_32FC1, 1, 0);
        cv::Scharr(left_scaled, igy_image, CV_32FC1, 0, 1);     

        float* depth_gradientX = new float[width*height]();
        float* depth_gradientY = new float[width*height]();
        float* depth_info = new float[width*height]();
        float* image_info = new float[width*height]();
        pcl::PointCloud<pcl::PointXYZ>::Ptr image_cloud (new pcl::PointCloud<pcl::PointXYZ>);
        image_cloud->width    = width;
        image_cloud->height   = height;
        image_cloud->is_dense = false;
        image_cloud->points.resize (image_cloud->width * image_cloud->height);

        
        for(size_t i=0; i<width*height;i++)
        {

            int u, v;
            u = i%width;
            v = i/width;

            //depth gradient
            depth_gradientX[i] = dgx_image.at<float>(v,u)/32.0f;
            depth_gradientY[i] = dgy_image.at<float>(v,u)/32.0f;

            //depth info
            float info_denom = sqrt(depth_gradientX[i]*depth_gradientX[i]+depth_gradientY[i]*depth_gradientY[i]);
            if (!isfinite(info_denom)) depth_info[i] = 0;
            else if (info_denom<0.001) depth_info[i] = 0;
            else depth_info[i] = 10.0/info_denom;
            float igx = igx_image.at<float>(v,u)/32.0f;
            float igy = igy_image.at<float>(v,u)/32.0f;
            if(i%2==0)image_info[i] = 1000.0f*sqrt(igx*igx+igy*igy);
            else image_info[i] = 0;
//            if(image_info[i]>500.0)image_info[i] = 0;
            
            //cloud plot
            if(depth_info[i]>0 && isfinite(depth[i])){
                image_cloud->points[i].x = depth[i]/K(0,0)*(u-K(0,2));
                image_cloud->points[i].y = depth[i]/K(1,1)*(v-K(1,2)); 
                image_cloud->points[i].z = depth[i];
            }    
              
        }

        pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
        if(frameID>2){
            //tracking
            update_pose = visual_tracking(ref_container,igx_container,igy_container,image_info,depth_raw,left_scaled,update_pose);
            cout<<update_pose<<endl;

            //prepare EST_pose, velo_cloud
            EST_pose = EST_pose*update_pose;
            if (mode == 0)pcl::transformPointCloud (*velo_raw, *velo_cloud, EST_pose.inverse());
            else
            {
                //extract local map
                pcl::PointXYZ searchPoint;
                searchPoint.x = EST_pose(0,3);
                searchPoint.y = EST_pose(1,3);
                searchPoint.z = EST_pose(2,3);
                std::vector<int> pointIdxRadiusSearch;
                std::vector<float> pointRadiusSquaredDistance;
                octree.radiusSearch (searchPoint, 30.0f, pointIdxRadiusSearch, pointRadiusSquaredDistance);

                cloud->width = pointIdxRadiusSearch.size ();
                cloud->height = 1;
                cloud->points.resize (cloud->width * cloud->height);
                int count = 0;
                for (size_t i = 0; i < pointIdxRadiusSearch.size (); ++i)
                {
                    cloud->points[count].x = velo_raw->points[ pointIdxRadiusSearch[i] ].x;
                    cloud->points[count].y = velo_raw->points[ pointIdxRadiusSearch[i] ].y;
                    cloud->points[count].z = velo_raw->points[ pointIdxRadiusSearch[i] ].z;
                    count++;
                }
//                //filtering
//                pcl::VoxelGrid<pcl::PointXYZ> sor;
//                sor.setInputCloud (cloud);
//                sor.setLeafSize (0.02f, 0.02f, 0.02f);
//                sor.filter (*cloud);

                pcl::transformPointCloud (*cloud, *velo_cloud, EST_pose.inverse());
            }
           
            //localization
            optimized_T = Matrix4f::Identity();
            optimized_T = Optimization(depth,depth_info,depth_gradientX,depth_gradientY);
            cout<<optimized_T<<endl;
            EST_pose = EST_pose*optimized_T.inverse();

//            if(mode == 0)debugImage(depth_image,dgx_image,dgy_image,depth_info);//save debug images
//            if(mode == 1)save_colormap(depth_image, "image_depth2.jpg",0,30);

            if(mode == 1) MapPub.PublishMap(cloud,1);//publish velo raw

        }

        
        if(mode == 0)MapPub.PublishMap(velo_raw,1);//publish velo raw
        if(mode == 0)MapPub.PublishPose(GT_pose,1);//publish GT pose
        MapPub.PublishPose(EST_pose,3);
        pcl::transformPointCloud (*image_cloud, *image_cloud, EST_pose);
        MapPub.PublishMap(image_cloud,3);
        
        //prepare reference images
        left_scaled.copyTo(ref_image);
        igx_image.copyTo(ref_igx);
        igy_image.copyTo(ref_igy);

        //save poses 
        write_poses("EST_poses.txt", EST_pose);
        if(mode == 0)write_poses("GT_poses.txt", GT_pose);

        if(mode == 1){
            Eigen::Affine3d e;
            e.matrix() = EST_pose.matrix().cast<double>();
            tf::transformEigenToTF(e, wtb);
        }     
        mTfBr.sendTransform(tf::StampedTransform(wtb,ros::Time::now(), "/CamLoc/World", "/CamLoc/Camera"));

        if(mode == 0)Velo_received = false;
        Left_received = false;
        Right_received = false;
        delete [] depth_gradientX;
        delete [] depth_gradientY;
        delete [] depth;
        delete [] depth_info;

        end_time = timestamp_now ();
        float time_diff = (end_time - start_time)/1000000.0;
        write_times("Elapsed_times.txt", time_diff);
        cout<<"Elapsed time: %"<<time_diff<<" secs\n"<<endl;

    }    



}

void CamLocalization::read_poses(std::string fname)
{
    int count = 0;
    GT_poses.clear();
//    GT_poses.reserve(4500);
    Matrix4f tmp = Matrix4f::Identity();
    ifstream file("poses.txt");
     while(!file.eof()){
        string s;
        for(int i = 0; i <3; ++i){
            getline(file, s, ' ');
            tmp(i,0) = atof(s.c_str());
            getline(file, s, ' ');
            tmp(i,1) = atof(s.c_str());
            getline(file, s, ' ');
            tmp(i,2) = atof(s.c_str());
            getline(file, s, '\n');
            tmp(i,3) = atof(s.c_str());                
        }
        GT_poses.push_back(tmp);
        tmp = Matrix4f::Identity();
    }   
    file.close();
}

void CamLocalization::write_poses(std::string fname, Matrix4f saved_pose)
{
    //write poses
    ofstream poses_file(fname, std::ios::app);
    if (poses_file.is_open()){
        poses_file << saved_pose(0,0) << ' ';
        poses_file << saved_pose(0,1) << ' ';
        poses_file << saved_pose(0,2) << ' ';
        poses_file << saved_pose(0,3) << '\n';
        poses_file << saved_pose(1,0) << ' ';
        poses_file << saved_pose(1,1) << ' ';
        poses_file << saved_pose(1,2) << ' ';
        poses_file << saved_pose(1,3) << '\n';
        poses_file << saved_pose(2,0) << ' ';
        poses_file << saved_pose(2,1) << ' ';
        poses_file << saved_pose(2,2) << ' ';
        poses_file << saved_pose(2,3) << '\n';
    }
    poses_file.close();
}

void CamLocalization::write_times(std::string fname, float time_diff)
{
    //write poses
    ofstream times_file(fname, std::ios::app);
    if (times_file.is_open()){
        times_file << time_diff << '\n';
    }
    times_file.close();
}

void CamLocalization::VeloPtsCallback(const sensor_msgs::PointCloud2::ConstPtr& msg)
{

    if(Velo_received==false)
    {
        pcl::PCLPointCloud2 pcl_pc2;
        pcl_conversions::toPCL(*msg,pcl_pc2);
        pcl::fromPCLPointCloud2(pcl_pc2,*velo_cloud);
        pcl::fromPCLPointCloud2(pcl_pc2,*velo_xyzi);       
        if (velo_cloud->points.size()>0){

            bool success;
            success = tlistener.waitForTransform("/kitti/World", "/kitti/Current", ros::Time(0), ros::Duration(0.1));
            if (success) {
                tlistener.lookupTransform("/kitti/World", "/kitti/Current", ros::Time(0), wtb);
                if(ODO_time!=wtb.stamp_){
                    Matrix4f prev_pose = GT_pose;
                    pcl_ros::transformAsMatrix (wtb, GT_pose);
//                    update_pose = prev_pose.inverse()*GT_pose;     
                    ODO_time=wtb.stamp_;

                    cout<<"Velodyne input: "<<velo_cloud->points.size()<<endl;
                    pcl::transformPointCloud (*velo_cloud, *velo_cloud, cTv);//transform to camera keyframe coordinate
                    Velo_received = true;
                }
            }

        }
    }
}

void CamLocalization::LeftImgCallback(const sensor_msgs::ImageConstPtr& msg, const sensor_msgs::CameraInfoConstPtr & infomsg)
{
    if(Left_received==false)
    {
        //image processing
        cv_bridge::CvImagePtr cv_ptr;
        cv_ptr = cv_bridge::toCvCopy(msg, sensor_msgs::image_encodings::MONO8);
        left_image = cv_ptr->image;       
        ancient_width = left_image.cols;
        cv::resize(left_image, left_image, cv::Size(), scale, scale);
        Left_received = true; 
        cv::Ptr<cv::CLAHE> clahe = cv::createCLAHE();
        clahe->setClipLimit(5);
        clahe->apply(left_image,left_image);
    
        //camera info
        P0 = Map<const MatrixXd>(&infomsg->P[0], 3, 4);         
    }
}

void CamLocalization::RightImgCallback(const sensor_msgs::ImageConstPtr& msg, const sensor_msgs::CameraInfoConstPtr & infomsg)
{
    if(Right_received==false)
    {
        //image processing
        cv_bridge::CvImagePtr cv_ptr;
        cv_ptr = cv_bridge::toCvCopy(msg, sensor_msgs::image_encodings::MONO8);
        right_image = cv_ptr->image;
        cv::resize(right_image, right_image, cv::Size(), scale, scale);
        Right_received = true; 
        cv::Ptr<cv::CLAHE> clahe = cv::createCLAHE();
        clahe->setClipLimit(5);
        clahe->apply(right_image,right_image);

        //camera info
        P1 = Map<const MatrixXd>(&infomsg->P[0], 3, 4);   
    }
}

//void CamLocalization::CamInfoCallback(const sensor_msgs::CameraInfo::ConstPtr& msg)
//{
//    //K= Map<const Matrix3d>(&msg->K[0], 3, 3);
//    R_rect = Map<const Matrix3d>(&msg->R[0], 3, 3);
//    P_rect = Map<const MatrixXd>(&msg->P[0], 3, 4);

////    cout<<"K Matrix: "<<K<<endl;
//}

Matrix4f CamLocalization::visual_tracking(const float* ref, const float* r_igx, const float* r_igy, const float* i_var, const float* idepth, cv::Mat cur, Matrix4f init_pose)
{
    cout<<"tracking start"<<endl;
    const float deltaHuber = sqrt(10);//10 may be the best choice
    //solver initialization
    g2o::SparseOptimizer optimizer;
    g2o::BlockSolverX::LinearSolverType * linearSolver;
    linearSolver = new g2o::LinearSolverDense<g2o::BlockSolverX::PoseMatrixType>();
    g2o::BlockSolverX * solver_ptr = new g2o::BlockSolverX(linearSolver);
    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
    optimizer.setAlgorithm(solver);

    // SET SIMILARITY VERTEX
    g2o::VertexSim3Expmap * vSim3 = new g2o::VertexSim3Expmap();
    vSim3->_fix_scale= false;
    Matrix4d Rt = init_pose.cast<double> ();
    Matrix3d R = Rt.block<3,3>(0,0);
    Vector3d t(Rt(0,3),Rt(1,3),Rt(2,3));
    const double s = 1;
    g2o::Sim3 g2oS_init(R,t,s);
    vSim3->setEstimate(g2oS_init);
    vSim3->setId(0);
    vSim3->setFixed(false);
    vSim3->_principle_point1[0] = K(0,2);
    vSim3->_principle_point1[1] = K(1,2);
    vSim3->_focal_length1[0] = K(0,0);
    vSim3->_focal_length1[1] = K(1,1);
    vSim3->_width = width;
    vSim3->_height = height;
    vSim3->ImageD = ref;
    vSim3->ImageGx = r_igx;
    vSim3->ImageGy = r_igy;
    vSim3->ImageInfo = i_var;
    optimizer.addVertex(vSim3);

    Matrix<double, 1, 1> info;
    info << 0.0f;

    int index = 1;
    for(size_t i=0; i<width*height;i++)
    {

        int u, v;
        u = i%width;
        v = i/width;

        // Set points
        Matrix<double,3,1> pts;
        pts <<idepth[i]/K(0,0)*(u-K(0,2)), idepth[i]/K(1,1)*(v-K(1,2)), idepth[i] ;

        Vector2d Ipos( vSim3->cam_map(vSim3->estimate().map(pts)) );
        int i_idx = ((int)Ipos[1])*vSim3->_width+((int)Ipos[0]);

        if ( pts[2]>0.0f && isfinite(pts[2]) && pts[2]<matching_thres){ //pts[2]<16*K(0,0)*base_line/100){ //pts[2]<30.0){//
            if (Ipos[0]<vSim3->_width && Ipos[0]>=0 && Ipos[1]<vSim3->_height && Ipos[1]>=0 && i_var[i_idx]>100)
            {
                // SET PointXYZ VERTEX
                g2o::VertexSBAPointXYZ* vPoint = new g2o::VertexSBAPointXYZ();
                vPoint->setEstimate(pts);
                vPoint->setId(index);
                vPoint->setFixed(true);
                optimizer.addVertex(vPoint);
                
                // Set Edges
                g2o::EdgeSim3ProjectXYZ* e01 = new g2o::EdgeSim3ProjectXYZ();
                e01->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(index)));
                e01->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));
                e01->setMeasurement(cur.at<float>(v,u));//(e1-e2);//(pts[i][2]);//
                info << i_var[i_idx];
//                info << 1000.0f;
                e01->setInformation(info);
                g2o::RobustKernelHuber* rk1 = new g2o::RobustKernelHuber;
                rk1->setDelta(deltaHuber);
                e01->setRobustKernel(rk1);

                optimizer.addEdge(e01);

                index++;
            }
        }

    }
    
    optimizer.initializeOptimization();
    optimizer.computeActiveErrors();
//    optimizer.setVerbose(true);
    int g2oresult = optimizer.optimize(10);

    // Recover optimized Sim3
    g2o::VertexSim3Expmap* vSim3_recov = static_cast<g2o::VertexSim3Expmap*>(optimizer.vertex(0));
    g2o::Sim3 g2oS12 = vSim3_recov->estimate();      

    Matrix4f result_mat = Sim3toMat(g2oS12);
//    Matrix4f result_mat;
    return result_mat;
}

Matrix4f CamLocalization::Optimization(const float* idepth, const float* idepth_var, const float* d_gradientX, const float* d_gradientY)
{

    //g2o optimization 
//    cout<<"g2o Optimization"<<endl;
    const float deltaHuber = sqrt(10);//10 may be the best choice

    //solver initialization
    g2o::SparseOptimizer optimizer;
    g2o::BlockSolverX::LinearSolverType * linearSolver;
    linearSolver = new g2o::LinearSolverDense<g2o::BlockSolverX::PoseMatrixType>();
    g2o::BlockSolverX * solver_ptr = new g2o::BlockSolverX(linearSolver);
    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
//    solver->setUserLambdaInit(100.0f);
    optimizer.setAlgorithm(solver);

    // SET SIMILARITY VERTEX
    g2o::VertexSim3Expmap * vSim3 = new g2o::VertexSim3Expmap();
    vSim3->_fix_scale= false;
    Matrix3d R = Matrix3d::Identity();
    Vector3d t(0,0,0);
    const double s = 1;
    g2o::Sim3 g2oS_init(R,t,s);
    vSim3->setEstimate(g2oS_init);
    vSim3->setId(0);
    vSim3->setFixed(false);
    vSim3->_principle_point1[0] = K(0,2);
    vSim3->_principle_point1[1] = K(1,2);
    vSim3->_focal_length1[0] = K(0,0);
    vSim3->_focal_length1[1] = K(1,1);
    vSim3->_width = width;
    vSim3->_height = height;
    vSim3->ImageD = idepth;
    vSim3->ImageGx = d_gradientX;
    vSim3->ImageGy = d_gradientY;
    vSim3->ImageInfo = idepth_var;
    optimizer.addVertex(vSim3);

    //Set map point vertices
    std::vector<Matrix<double,3,1> > pts;
    int numpts = velo_cloud->points.size();
    pts.reserve(numpts);

    Matrix<double, 1, 1> info;
    info << 0.0f;

    int index = 1;
    for(size_t i=0; i<numpts;i++)
    {
        // Set map points
        pts[i] <<velo_cloud->points[i].x, velo_cloud->points[i].y, velo_cloud->points[i].z ;

        Vector2d Ipos( vSim3->cam_map(vSim3->estimate().map(pts[i])) );
        int i_idx = ((int)Ipos[1])*vSim3->_width+((int)Ipos[0]);
        
        
        if ( pts[i][2]>0.0f && pts[i][2]<matching_thres){ //pts[i][2]<16*K(0,0)*base_line/100){//pts[i][2]<30.0){//
                if (Ipos[0]<vSim3->_width && Ipos[0]>=0 && Ipos[1]<vSim3->_height && Ipos[1]>=0 && idepth_var[i_idx]>0)
                {
                    // SET PointXYZ VERTEX
//                    pts[i][2] = 11.0f;
                    g2o::VertexSBAPointXYZ* vPoint = new g2o::VertexSBAPointXYZ();
                    vPoint->setEstimate(pts[i]);
                    vPoint->setId(index);
                    vPoint->setFixed(true);
                    optimizer.addVertex(vPoint);
                    
                    // Set Edges
                    g2o::EdgeSim3ProjectXYZD* e01 = new g2o::EdgeSim3ProjectXYZD();
                    e01->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(index)));
                    e01->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));
                    e01->setMeasurement(1.0f);//(e1-e2);//(pts[i][2]);//
                    info << idepth_var[i_idx];
//                    info << 1000.0f;
                    e01->setInformation(info);
                    g2o::RobustKernelHuber* rk1 = new g2o::RobustKernelHuber;
                    rk1->setDelta(deltaHuber);
                    e01->setRobustKernel(rk1);

                    optimizer.addEdge(e01);

                    index++;
                }
        }

    }
    
    cout<<index<<endl;

    optimizer.initializeOptimization();
    optimizer.computeActiveErrors();

//    optimizer.setVerbose(true);
    int g2oresult = optimizer.optimize(100);
    
    cout<<g2oresult<<endl;

    // Recover optimized Sim3
    g2o::VertexSim3Expmap* vSim3_recov = static_cast<g2o::VertexSim3Expmap*>(optimizer.vertex(0));
    g2o::Sim3 g2oS12 = vSim3_recov->estimate();      

    Matrix4f result_mat = Sim3toMat(g2oS12);

    return result_mat;

    
}

void CamLocalization::debugImage(cv::Mat& depth_image,cv::Mat& dgx_image,cv::Mat& dgy_image, const float* depth_info)
{
    //generate lidar intensity image 
    int numpts = velo_cloud->points.size();    
    cv::Mat frame = cv::Mat::zeros(cv::Size(left_image.cols, left_image.rows), CV_8UC1);
    cv::Mat depth = cv::Mat::zeros(cv::Size(left_image.cols, left_image.rows), CV_32FC1);
    for(size_t i=0; i<numpts;i++)
    {
        if (velo_cloud->points[i].z>-1.0){
            // Set frame and depth image
            Vector3d point(velo_cloud->points[i].x, velo_cloud->points[i].y, velo_cloud->points[i].z);
            Vector2d uv = ProjectTo2D(point);
            int u = (int) uv(0);
            int v = (int) uv(1);
            
            if (u<left_image.cols && u>=0 && v<left_image.rows && v>=0 ){
                if (depth.at<float>(v,u)>0){
                    if(depth.at<float>(v,u)>velo_cloud->points[i].z){
                        depth.at<float>(v,u) = velo_cloud->points[i].z;
                        frame.at<int8_t>(v,u) = (int)255.0f*velo_xyzi->points[i].intensity;
                    }
                }
                else{
                    depth.at<float>(v,u) = velo_cloud->points[i].z;
                    frame.at<int8_t>(v,u) = (int)255.0f*velo_xyzi->points[i].intensity;
                }
                
                
            }
            
        }
        
    }


    cv::Mat mag = cv::Mat::zeros(cv::Size(left_image.cols, left_image.rows), CV_32FC1);
    cv::Mat dir = cv::Mat::zeros(cv::Size(left_image.cols, left_image.rows), CV_32FC1);
    cv::Mat info = cv::Mat::zeros(cv::Size(left_image.cols, left_image.rows), CV_32FC1);
    float delta = 10.0;
    float max_info = 0;
    for(size_t i=0; i<width*height;i++)
    {
        int u, v;
        u = i%width;
        v = i/width;

        float v1 = 0 ,v2 = 0;
        if(isfinite(dgx_image.at<float>(v,u)))
        {
            if(dgx_image.at<float>(v,u)>delta)v1 = delta;
            else if(dgx_image.at<float>(v,u)<-delta)v1 = -delta;
            else v1 = dgx_image.at<float>(v,u);
        }
        if(isfinite(dgy_image.at<float>(v,u)))
        {
            if(dgy_image.at<float>(v,u)>delta)v2 = delta;
            else if(dgy_image.at<float>(v,u)<-delta)v2 = -delta;
            else v2 = dgy_image.at<float>(v,u);
        }
        mag.at<float>(v,u) = sqrt(v1*v1+v2*v2);
        dir.at<float>(v,u) = atan2(v2,v1);
        float tmp = dir.at<float>(v,u); 
        if(tmp<0)dir.at<float>(v,u) = tmp+2*M_PI;

        info.at<float>(v,u) = depth_info[i];
        if(depth_info[i]>max_info)max_info = depth_info[i];
    }
    

    cv::Mat intensity_residual = cv::Mat::zeros(cv::Size(left_image.cols, left_image.rows), CV_32FC1);
    cv::Mat depth_residual = cv::Mat::zeros(cv::Size(left_image.cols, left_image.rows), CV_32FC1);

    for(size_t i=0; i<numpts;i++)
    {
        if (velo_cloud->points[i].z>-1.0){
            Vector3d point(velo_cloud->points[i].x, velo_cloud->points[i].y, velo_cloud->points[i].z);
            Vector2d uv = ProjectTo2D(point);
            int u = (int) uv(0);
            int v = (int) uv(1);
            if (u<left_image.cols && u>=0 && v<left_image.rows && v>=0 ){
                float d_float;
                d_float = depth_image.at<float>(v,u);
                if(d_float>0 && d_float==d_float && d_float<100)
                {
                }
                else
                {
                    d_float = 0;
                }
                depth_residual.at<float>(v,u) = abs(depth.at<float>(v,u) - d_float);
                intensity_residual.at<float>(v,u) = abs((float)frame.at<int8_t>(v,u) - (float)left_image.at<int8_t>(v,u));
            }
        }
    }     

    cout<<max_info<<endl;
    
    save_colormap(frame,"lidar_intensity.jpg",0,255);
    save_colormap(left_image,"image_intensity.jpg",0,255);
    save_colormap(depth,"lidar_depth.jpg",0,30);
    save_colormap(depth_image, "image_depth.jpg",0,30);
    save_colormap(intensity_residual, "intensity_residual.jpg",0,255);
    save_colormap(depth_residual, "depth_residual.jpg",0,30);
    save_colormap2d(dir,mag, "gradient.jpg");
    save_colormap(info, "information.jpg",0,100);
   
}


