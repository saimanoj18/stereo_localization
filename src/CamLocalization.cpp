#include "CamLocalization.h"

void CamLocalization::CamLocInitialize(cv::Mat image)
{
    //Set matrices
    float fx = 718.856*image.cols/ancient_width;
    float fy = 718.856*image.cols/ancient_width;        
    float cx = 607.1928*image.cols/ancient_width;
    float cy = 185.2157*image.cols/ancient_width;
    K << fx, 0.0, cx, 0.0, fy, cy, 0.0, 0.0, 1.0;
    cout<<K<<endl;
    cout<<left_image.cols<<", "<<left_image.rows<<endl;
    width = left_image.cols;
    height = left_image.rows;
    
    //Set ref images
    ref_container = new float[width*height];
    igx_container = new float[width*height];
    igy_container = new float[width*height];        

}

void CamLocalization::Refresh()
{
    bool success;
    success = tlistener.waitForTransform("/kitti/World", "/kitti/Velodyne", ros::Time(0), ros::Duration(0.1));
    if (success) {
        tlistener.lookupTransform("/kitti/World", "/kitti/Velodyne", ros::Time(0), ctv);
        pcl_ros::transformAsMatrix (ctv, cTv); 
    }
    

    if(Velo_received && Left_received && Right_received)
    {
        start_time = timestamp_now ();

        //prepar GT pose and map point clouds        
//        int frame_half = (int) frameID/2;
//        cout<<frame_half<<endl;
//        GT_pose = GT_poses[frame_half];//ODO_pose;//
        pcl::transformPointCloud (*velo_cloud, *velo_raw, GT_pose);//transform to world coordinate
        mTfBr.sendTransform(tf::StampedTransform(wtb,ros::Time::now(), "/CamLoc/World", "/CamLoc/Camera"));
        MapPub.PublishMap(velo_raw,1);//publish velo raw
        MapPub.PublishPose(GT_pose,1);//publish GT pose

        //initialize 
        if(frameID == 0)CamLocInitialize(right_image);
        
        /////////////////////////disparity map generation/////////////////////////////        
        cv::Mat disp = cv::Mat::zeros(cv::Size(width, height), CV_16S);
        cv::Ptr<cv::StereoSGBM> sbm = cv::StereoSGBM::create(0,16*2,9);
        sbm->compute(left_image, right_image, disp);
        frameID = frameID+2;

        /////////////////////////depth image generation/////////////////////////////
        float* depth = new float[width*height]();
        cv::Mat depth_image = cv::Mat(cv::Size(left_image.cols, left_image.rows), CV_32FC1);
        for(size_t i=0; i<width*height;i++)
        {
            int u, v;
            u = i%width;
            v = i/width;
            
            //depth
            int16_t d = disp.at<int16_t>(v,u);
            if(d==0 || d!=d || d<100) d = 0;
            depth[i] = 16*K(0,0)*0.54/((float)d);//0.54*K(0,0)/disp2.at<float>(v,u);

            //depth image            
            depth_image.at<float>(v,u) = depth[i];
            
            //reference images
            if(frameID>2){
                ref_container[i] = ref_image.at<char>(v,u); 
                igx_container[i] = ref_igx.at<float>(v,u);
                igy_container[i] = ref_igy.at<float>(v,u);
            }

        }

        /////////////////////////depth gradient generation/////////////////////////////
        cv::Mat dgx_image, dgy_image; // = cv::Mat(cv::Size(left_image.cols, left_image.rows), CV_8UC3);
        cv::Scharr(depth_image, dgx_image, CV_32FC1, 1, 0);
        cv::Scharr(depth_image, dgy_image, CV_32FC1, 0, 1);
        cv::Mat igx_image, igy_image;
        cv::Scharr(left_image, igx_image, CV_32FC1, 1, 0);
        cv::Scharr(left_image, igy_image, CV_32FC1, 0, 1);        

        float* depth_gradientX = new float[width*height]();
        float* depth_gradientY = new float[width*height]();
        float* depth_info = new float[width*height]();
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
            float info_denom = depth_gradientX[i]*depth_gradientX[i]+depth_gradientY[i]*depth_gradientY[i];
            if (!isfinite(info_denom)) depth_info[i] = 0;
            else if (info_denom<0.001) depth_info[i] = 1000.0;
            else depth_info[i] = 10.0/(info_denom);


            //cloud plot
            if(depth_info[i]>0 && isfinite(depth[i])){
                image_cloud->points[i].x = depth[i]/K(0,0)*(u-K(0,2));
                image_cloud->points[i].y = depth[i]/K(1,1)*(v-K(1,2)); 
                image_cloud->points[i].z = depth[i];
            }    
              
        }


        if(frameID>2){
//            //tracking
//            update_pose = Matrix4f::Identity();
//            update_pose = visual_tracking(ref_container,igx_container,igy_container,depth,depth_info,left_image);
//            cout<<update_pose<<endl;

            //prepare EST_pose, velo_cloud
            EST_pose = EST_pose*update_pose;
            pcl::transformPointCloud (*velo_raw, *velo_cloud, EST_pose.inverse());

            //localization
            optimized_T = Matrix4f::Identity();
            optimized_T = Optimization(depth,depth_info,depth_gradientX,depth_gradientY);
            pcl::transformPointCloud (*velo_cloud, *velo_cloud, optimized_T);
            EST_pose = EST_pose*optimized_T.inverse();

        }

        MapPub.PublishPose(EST_pose,3);
        pcl::transformPointCloud (*image_cloud, *image_cloud, EST_pose);
        MapPub.PublishMap(image_cloud,3);
        
        //prepare reference images
        left_image.copyTo(ref_image);
        igx_image.copyTo(ref_igx);
        igy_image.copyTo(ref_igy);

        //save poses 
        write_poses("EST_poses.txt", EST_pose);
        write_poses("GT_poses.txt", GT_pose);
//        write_poses("ODO_poses.txt", ODO_pose);

        
        mTfBr.sendTransform(tf::StampedTransform(wtb,ros::Time::now(), "/CamLoc/World", "/CamLoc/Camera"));

        Velo_received = false;
        Left_received = false;
        Right_received = false;
        delete [] depth_gradientX;
        delete [] depth_gradientY;
        delete [] depth;
        delete [] depth_info;

        end_time = timestamp_now ();

        cout<<"Elapsed time: %"<<(end_time - start_time)/1000000.0<<" secs\n"<<endl;

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

void CamLocalization::VeloPtsCallback(const sensor_msgs::PointCloud2::ConstPtr& msg)
{

    if(Velo_received==false)
    {
        pcl::PCLPointCloud2 pcl_pc2;
        pcl_conversions::toPCL(*msg,pcl_pc2);
        pcl::fromPCLPointCloud2(pcl_pc2,*velo_cloud);       
        if (velo_cloud->points.size()>0){

            bool success;
            success = tlistener.waitForTransform("/kitti/World", "/kitti/Current", ros::Time(0), ros::Duration(0.1));
            if (success) {
                tlistener.lookupTransform("/kitti/World", "/kitti/Current", ros::Time(0), wtb);
                if(ODO_time!=wtb.stamp_){
                    Matrix4f prev_pose = GT_pose;
                    pcl_ros::transformAsMatrix (wtb, GT_pose);
                    update_pose = prev_pose.inverse()*GT_pose;     
                    ODO_time=wtb.stamp_;

                    cout<<"Velodyne input: "<<velo_cloud->points.size()<<endl;
                    pcl::transformPointCloud (*velo_cloud, *velo_cloud, cTv);//transform to camera keyframe coordinate
                    Velo_received = true;
                }
            }

        }
    }
}

void CamLocalization::LeftImgCallback(const sensor_msgs::Image::ConstPtr& msg)
{
    if(Left_received==false)
    {
        cv_bridge::CvImagePtr cv_ptr;
        cv_ptr = cv_bridge::toCvCopy(msg, sensor_msgs::image_encodings::MONO8);
        left_image = cv_ptr->image;       
        ancient_width = left_image.cols;
        cv::resize(left_image, left_image, cv::Size(), scale, scale);
        Left_received = true; 
        cv::Ptr<cv::CLAHE> clahe = cv::createCLAHE();
        clahe->setClipLimit(5);
        clahe->apply(left_image,left_image);         
    }
}

void CamLocalization::RightImgCallback(const sensor_msgs::Image::ConstPtr& msg)
{
    if(Right_received==false)
    {
        cv_bridge::CvImagePtr cv_ptr;
        cv_ptr = cv_bridge::toCvCopy(msg, sensor_msgs::image_encodings::MONO8);
        right_image = cv_ptr->image;
        cv::resize(right_image, right_image, cv::Size(), scale, scale);
        Right_received = true; 
        cv::Ptr<cv::CLAHE> clahe = cv::createCLAHE();
        clahe->setClipLimit(5);
        clahe->apply(right_image,right_image);  
    }
}

void CamLocalization::CamInfoCallback(const sensor_msgs::CameraInfo::ConstPtr& msg)
{
    //K= Map<const Matrix3d>(&msg->K[0], 3, 3);
    R_rect = Map<const Matrix3d>(&msg->R[0], 3, 3);
    P_rect = Map<const MatrixXd>(&msg->P[0], 3, 4);

//    cout<<"K Matrix: "<<K<<endl;
}

Matrix4f CamLocalization::visual_tracking(const float* ref, const float* r_igx, const float* r_igy, const float* i_var, const float* idepth, cv::Mat cur)
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

        if (Ipos[0]<vSim3->_width && Ipos[0]>=0 && Ipos[1]<vSim3->_height && Ipos[1]>=0 && i_var[i_idx]>0)
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
            e01->setMeasurement(cur.at<char>(v,u));//(e1-e2);//(pts[i][2]);//
            info << i_var[i_idx];
//            info << 1000.0f;
            e01->setInformation(info);
            g2o::RobustKernelHuber* rk1 = new g2o::RobustKernelHuber;
            rk1->setDelta(deltaHuber);
            e01->setRobustKernel(rk1);

            optimizer.addEdge(e01);

            index++;
        }

    }
    
    optimizer.initializeOptimization();
    optimizer.computeActiveErrors();
    optimizer.setVerbose(true);
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
        
        
        if ( pts[i][2]>0.0f && pts[i][2]<16*K(0,0)*0.54/100){
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
    
    optimizer.initializeOptimization();
    optimizer.computeActiveErrors();

//    optimizer.setVerbose(true);
    int g2oresult = optimizer.optimize(10);

    // Recover optimized Sim3
    g2o::VertexSim3Expmap* vSim3_recov = static_cast<g2o::VertexSim3Expmap*>(optimizer.vertex(0));
    g2o::Sim3 g2oS12 = vSim3_recov->estimate();      

    Matrix4f result_mat = Sim3toMat(g2oS12);

    return result_mat;

    
}
