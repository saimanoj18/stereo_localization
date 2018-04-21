#include "CamLocalization.h"

void CamLocalization::CamLocInitialize(cv::Mat image)
{
    //Set matrices
    float fx = P0(0,0)*image.cols/ancient_width;
    float fy = P0(1,1)*image.cols/ancient_width;        
    float cx = P0(0,2)*image.cols/ancient_width;
    float cy = P0(1,2)*image.cols/ancient_width;
    base_line = -P1(0,3)/P0(0,0);// because R = I;

    cout<<"width: "<<ancient_width<<endl;
    cout<<"base_line: "<<base_line<<endl;
    K << fx, 0.0, cx, 0.0, fy, cy, 0.0, 0.0, 1.0;
    cout<<"K_mat: "<<K<<endl;
    cout<<"extrinsic Camera to Velodyne: "<<cTv<<endl;
    cout<<left_image.cols<<", "<<left_image.rows<<endl;
    width = left_image.cols;
    height = left_image.rows;
    
    //Set ref images
    ref_container = new float[width*height];
    src_container = new float[width*height];
    image_gradientX = new float[width*height];
    image_gradientY = new float[width*height];
    image_info= new float[width*height]();

    //Set ref depth
    ref_depth = cv::Mat::zeros(cv::Size(left_image.cols, left_image.rows), CV_32FC1);
    ref_depth_info = cv::Mat::zeros(cv::Size(left_image.cols, left_image.rows), CV_32FC1);

    //Set informations
    depth = new float[width*height]();
    depth_gradientX = new float[width*height]();
    depth_gradientY = new float[width*height]();
    depth_info = new float[width*height]();

    //set matching thres
    if(mode == 0){
        //set initial pose
        EST_pose = GT_pose;   
        d_var = 0.01;
        d_limit = 0.0;
        matching_thres = K(0,0)*base_line*( 1.0/(d_limit/16.0) + d_var/((float)(d_limit/16.0)*(d_limit/16.0)*(d_limit/16.0)) );
    }

    if (mode ==1)
    {
    
        //set initial pose       
        IN_pose = GT_pose*cTv.inverse();
        Matrix4d IN_pose_inv = IN_pose.inverse();
        EST_pose = Matrix4d::Identity();

        d_var = 0.01;
        d_limit = 50.0;
        matching_thres = K(0,0)*base_line*( 1.0/(d_limit/16.0) + d_var/((float)(d_limit/16.0)*(d_limit/16.0)*(d_limit/16.0)) );

        //load velo_global from .las
        std:string filename;
        filename = data_path_+"/sequences/00/sick_pointcloud.las";
        std::ifstream ifs;
        if (!liblas::Open(ifs, filename))
        {
            throw std::runtime_error(std::string("Can not open ") + filename);
        }

        liblas::Reader reader(ifs);
        liblas::Header const& h = reader.GetHeader();        
        velo_global->width = h.GetPointRecordsCount();
        velo_global->height = 1;
        velo_global->points.resize (velo_global->width * velo_global->height);
        int count = 0;
        int iter = 0;
        while (reader.ReadNextPoint())
        {        
            liblas::Point const& p = reader.GetPoint();

            velo_global->points[count].x = p[0];
            velo_global->points[count].y = p[1];
            velo_global->points[count].z = p[2];

            count++;
            iter++;
        }
        
        pcl::transformPointCloud (*velo_global, *velo_global, IN_pose.inverse().matrix().cast <float> ());

        //set octree
        octree.setInputCloud (velo_global);
        octree.addPointsFromInputCloud ();
        
        
    }        

}

void CamLocalization::Refresh()
{
    bool success;
    success = tlistener.waitForTransform("/kitti/World", "/kitti/Velodyne", ros::Time(0), ros::Duration(0.1));
    if (success) {
        tlistener.lookupTransform("/kitti/World", "/kitti/Velodyne", ros::Time(0), ctv);
        Eigen::Affine3d e_temp;
        tf::transformTFToEigen (ctv, e_temp);
        cTv.matrix() = e_temp.matrix();       
    }
    
    if(mode ==1) Velo_received = true;


    if(Velo_received && Left_received && Right_received)
    {
          
        int u, v;//for loops

        start_time = timestamp_now ();

        //prepare GT pose and map point clouds 
        bool success_pose;
        success_pose = tlistener.waitForTransform("/kitti/World", "/kitti/Current", ros::Time(0), ros::Duration(0.1));
        if (success_pose) {
            tlistener.lookupTransform("/kitti/World", "/kitti/Current", ros::Time(0), wtb);
            Eigen::Affine3d e_temp;
            tf::transformTFToEigen(wtb, e_temp);
            GT_pose.matrix() = e_temp.matrix();  
        }

        //initialize 
        if(frameID == 0)CamLocInitialize(right_image);            

        //prepare image gradients & depth gradients
        cv::normalize(left_image, left_scaled, 0, 1, CV_MINMAX, CV_32FC1);   
        cv::Scharr(left_scaled, igx_image, CV_32FC1, 1, 0);
        cv::Scharr(left_scaled, igy_image, CV_32FC1, 0, 1);   

        /////////////////////////disparity map generation/////////////////////////////        
        disp = cv::Mat::zeros(cv::Size(width, height), CV_16S);
        cv::Ptr<cv::StereoSGBM> sbm;
        if(mode == 0)sbm = cv::StereoSGBM::create(0,16*5,7);
        if(mode == 1)sbm = cv::StereoSGBM::create(0,16*5,7);
        sbm->compute(left_image, right_image, disp);
        frameID = frameID+1;

        /////////////////////////depth image generation/////////////////////////////
        pcl::PointCloud<pcl::PointXYZ>::Ptr image_cloud (new pcl::PointCloud<pcl::PointXYZ>);
        image_cloud->width    = width;
        image_cloud->height   = height;
        image_cloud->is_dense = false;
        image_cloud->points.resize (image_cloud->width * image_cloud->height);
        for(size_t i=0; i<width*height;i++)
        {
            u = i%width;
            v = i/width;
            
            int16_t d = disp.at<int16_t>(v,u);            
            if(d==0 || d!=d || d<d_limit) d = 0; //
            depth[i] = K(0,0)*base_line*( 1.0/((float)d/16.0) + d_var/((float)(d/16.0)*(d/16.0)*(d/16.0)) );

            //depth image            
            ref_depth.at<float>(v,u) = depth[i];

            //image gradient
            float igx = igx_image.at<float>(v,u)/32.0f;
            float igy = igy_image.at<float>(v,u)/32.0f;
            float info_nom = sqrt(igx*igx+igy*igy);
            if(!isfinite(info_nom))image_info[i] = 0;
            else image_info[i] = 10.0f*sqrt(igx*igx+igy*igy);       

            //reference images
            if(frameID>1){
                ref_container[i] = ref_image.at<float>(v,u);
                src_container[i] = left_scaled.at<float>(v,u);
                image_gradientX[i] = igx_image.at<float>(v,u)/32.0f;
                image_gradientY[i] = igy_image.at<float>(v,u)/32.0f;
            } 

        }

//        std::string depth_file;
//        depth_file =  "./depth/"+std::to_string(frameID-1)+".jpg";
//        save_colormap(depth_image, depth_file,0,30);

//        cv::imshow("depth_image", depth_image);
//        cv::waitKey(3);
//        cv::moveWindow("depth_image", 50,20);



        if(mode == 1){
            GT_pose = IN_pose.inverse()*GT_pose*cTv.inverse();
        }
        if(frameID>1){                   

            bool success_imu;
            success_imu = tlistener.waitForTransform("/kitti/World", "/kitti/IMU", ros::Time(0), ros::Duration(0.1));
            if (success_imu) {
                tlistener.lookupTransform("/kitti/World", "/kitti/IMU", ros::Time(0), wtb);
                Eigen::Affine3d e_temp;
                tf::transformTFToEigen(wtb, e_temp);
                update_pose.matrix() = e_temp.matrix(); 
            }            
            update_pose = cTv*update_pose*cTv.inverse();


            /////////////////////////depth gradient generation/////////////////////////////
            //depth gradient
            cv::Scharr(ref_depth, dgx_image, CV_32FC1, 1, 0);
            cv::Scharr(ref_depth, dgy_image, CV_32FC1, 0, 1);
            int count_gradient = 0; 
            for(size_t i=0; i<width*height;i++)
            {
                u = i%width;
                v = i/width;

                //depth gradient
                depth_gradientX[i] = dgx_image.at<float>(v,u)/32.0f;
                depth_gradientY[i] = dgy_image.at<float>(v,u)/32.0f;

                //depth info
                float info_denom = sqrt(depth_gradientX[i]*depth_gradientX[i]+depth_gradientY[i]*depth_gradientY[i]);
                if (!isfinite(info_denom)) depth_info[i] = 0;
                else if (info_denom<0.01) depth_info[i] = 0;
                else depth_info[i] = 100.0/info_denom;

                //cloud plot
                if(isfinite(depth[i])){
                    image_cloud->points[i].x = depth[i]/K(0,0)*(u-K(0,2));
                    image_cloud->points[i].y = depth[i]/K(1,1)*(v-K(1,2)); 
                    image_cloud->points[i].z = depth[i];
                }   
            }

            //prepare velo_raw
            EST_pose = EST_pose*update_pose;
            if(mode == 0)pcl::transformPointCloud (*velo_cloud, *velo_raw, GT_pose.matrix().cast <float> ());//transform to world coordinate
            else
            {
                //extract local map
                pcl::PointXYZ searchPoint;
                Eigen::Matrix4f pose_f = EST_pose.matrix().cast<float> ();
                searchPoint.x = pose_f(0,3);
                searchPoint.y = pose_f(1,3);
                searchPoint.z = pose_f(2,3);
                std::vector<int> pointIdxRadiusSearch;
                std::vector<float> pointRadiusSquaredDistance;
                octree.radiusSearch (searchPoint, 30.0f, pointIdxRadiusSearch, pointRadiusSquaredDistance);

                velo_raw->clear();
                velo_raw->width = pointIdxRadiusSearch.size()/10+1;
                velo_raw->height = 1;
                velo_raw->points.resize (velo_raw->width * velo_raw->height);
                int count = 0;
                for (size_t i = 0; i < pointIdxRadiusSearch.size (); ++i)
                {
                    if(i%10==0){
                    velo_raw->points[count].x = velo_global->points[ pointIdxRadiusSearch[i] ].x;
                    velo_raw->points[count].y = velo_global->points[ pointIdxRadiusSearch[i] ].y;
                    velo_raw->points[count].z = velo_global->points[ pointIdxRadiusSearch[i] ].z;
                    count++;
                    }
                }
//                //filtering
//                pcl::VoxelGrid<pcl::PointXYZ> sor;
//                sor.setInputCloud (velo_raw);
//                sor.setLeafSize (0.01f, 0.01f, 0.01f);
//                sor.filter (*velo_raw);

            }

            //prepare velo_cloud
            pcl::transformPointCloud (*velo_raw, *velo_cloud, EST_pose.inverse().matrix().cast <float> ());
            
            //localization
            optimized_T = Matrix4d::Identity();
            optimized_T = Optimization_combined(ref_container, src_container, image_info, image_gradientX, image_gradientY, depth, depth_info, depth_gradientX, depth_gradientY, update_pose);
            cout<<optimized_T<<endl;
//            update_pose = update_pose*optimized_T.inverse();
            EST_pose = EST_pose*optimized_T.inverse();
            
//            if(mode == 0)debugImage(depth_image,dgx_image,dgy_image,depth_info);//save debug images
//            if(mode == 1)save_colormap(depth_image, "image_depth2.jpg",0,30);


        }
        
        //publish map and pose
        MapPub.PublishMap(velo_raw,1);//publish velo raw
        MapPub.PublishPose(GT_pose,1);//publish GT pose
        MapPub.PublishPose(EST_pose,3);
        pcl::transformPointCloud (*image_cloud, *image_cloud, EST_pose.matrix().cast <float> ());
        MapPub.PublishMap(image_cloud,3);
        
        //prepare reference images
        left_scaled.copyTo(ref_image);

        //save poses 
        write_poses("EST_poses.txt", EST_pose);
        write_poses("GT_poses.txt", GT_pose);
        
        //broadcast
        Eigen::Affine3d e;
        e.matrix() = EST_pose.matrix();
        tf::transformEigenToTF(e, wtb);
        mTfBr.sendTransform(tf::StampedTransform(wtb,ros::Time::now(), "/CamLoc/World", "/CamLoc/Camera"));

        if(mode == 0)Velo_received = false;
        Left_received = false;
        Right_received = false;

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
    Matrix4d tmp = Matrix4d::Identity();
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
        tmp = Matrix4d::Identity();
    }   
    file.close();
}

void CamLocalization::write_poses(std::string fname, Matrix4d saved_pose)
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

        cout<<"Velodyne input: "<<velo_cloud->points.size()<<endl;
        pcl::transformPointCloud (*velo_cloud, *velo_cloud, cTv.matrix().cast <float> ());//transform to camera keyframe coordinate
        Velo_received = true;

        }
    }
}

void CamLocalization::EncoderCallback(const irp_sen_msgs::encoder::ConstPtr& msg)
{
    int64_t cur_enc_left = msg->left_count;
    int64_t cur_enc_right = msg->right_count;

    if(frameID >0 ){
        int64_t L_diff = cur_enc_left - prev_enc_left;
        int64_t R_diff = cur_enc_right - prev_enc_right;
        dL = ((double)L_diff/((double)ENCODER_RESOLUTION)) * LEFT_DIAMETER * M_PI;
        dR = ((double)R_diff/((double)ENCODER_RESOLUTION)) * RIGHT_DIAMETER * M_PI;
//        double dist = 0.5*(dL+dR);        
//        double dth = (dR-dL)/ VEHICLE_THREAD;
//        double dx = dist * cos(dth);
//        double dy = dist * sin(dth); 
//        update_pose = Matrix4d::Identity();
//        update_pose(0,3) = dx;
//        update_pose(1,3) = dy;
//        AngleAxisd heading(dth,Eigen::Vector3d::UnitZ());
//        update_pose.block<3,3>(0,0) = heading.toRotationMatrix();
//        cout<<update_pose<<endl;
//        update_pose = cTv*update_pose*cTv.inverse();

    }
    prev_enc_left = cur_enc_left;
    prev_enc_right = cur_enc_right;
        
}

void CamLocalization::FogCallback(const irp_sen_msgs::fog_3axis::ConstPtr& msg)
{
    fog_angles<< msg->d_roll, msg->d_pitch, msg->d_yaw;
//    if(frameID >0 ){
//        AngleAxisd roll(cur_imu(0), Vector3d::UnitX());
//        AngleAxisd pitch(cur_imu(1), Vector3d::UnitY());
//        AngleAxisd yaw(0, Vector3d::UnitZ());
//        Matrix3d cur_mat = (yaw*pitch*roll).matrix();
//    
//        AngleAxisd p_roll(prev_imu(0), Vector3d::UnitX());
//        AngleAxisd p_pitch(prev_imu(1), Vector3d::UnitY());
//        AngleAxisd p_yaw(0, Vector3d::UnitZ());
//        Matrix3d prev_mat = (p_yaw*p_pitch*p_roll).matrix();
//        
//        update_angular_pose = Matrix4d::Identity();
//        update_angular_pose.block<3,3>(0,0) = prev_mat.inverse()*cur_mat;

//    }
//    prev_imu = cur_imu;
    
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


Matrix4d CamLocalization::Optimization_combined(const float* ref, const float* image, const float* image_var, float* i_gradientX, const float* i_gradientY, const float* idepth, const float* idepth_var, const float* d_gradientX, const float* d_gradientY, Matrix4d init_pose)
{
    //g2o optimization 
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
    vSim3->Image = image;
    vSim3->ImageGx = i_gradientX;
    vSim3->ImageGy = i_gradientY;
    vSim3->ImageInfo = image_var;
    vSim3->Depth = idepth;
    vSim3->DepthGx = d_gradientX;
    vSim3->DepthGy = d_gradientY;
    vSim3->DepthInfo = idepth_var;
    float* occ_container = new float[width*height]();
    vSim3->occ_image = occ_container;
    int* occ_idx = new int[width*height]();
    vSim3->occ_idx = occ_idx;
    optimizer.addVertex(vSim3);

    g2o::VertexSim3Expmap * vSim3_fake = new g2o::VertexSim3Expmap();
    R = init_pose.block<3,3>(0,0);
    t<<init_pose(0,3),init_pose(1,3),init_pose(2,3);
    g2o::Sim3 g2oS_init2(R,t,s);
    vSim3_fake->setEstimate(g2oS_init2);
    vSim3_fake->_principle_point1[0] = K(0,2);
    vSim3_fake->_principle_point1[1] = K(1,2);
    vSim3_fake->_focal_length1[0] = K(0,0);
    vSim3_fake->_focal_length1[1] = K(1,1);
    vSim3_fake->_width = width;
    vSim3_fake->_height = height;

    //Set map point vertices
    Matrix<double,3,1>pts;
    int numpts = velo_cloud->points.size();
    Matrix<double, 1, 1> info;
    info << 0.0f;

    int index = 1;
    for(size_t i=0; i<numpts;i++)
    {
        // map points for image n
        pts <<velo_cloud->points[i].x, velo_cloud->points[i].y, velo_cloud->points[i].z ;
        Vector2d Ipos( vSim3->cam_map(vSim3->estimate().map(pts)) );
        int i_idx = ((int)Ipos[1])*vSim3->_width+((int)Ipos[0]);

        // map points for image n-1
        Vector2d Ipos2( vSim3->cam_map(vSim3_fake->estimate().map(pts)) );
        int i_idx2 = ((int)Ipos2[1])*vSim3->_width+((int)Ipos2[0]);         
        
        
        if ( pts[2]>0.0f && pts[2]<matching_thres){ 
                if (Ipos[0]<vSim3->_width && Ipos[0]>=0 && Ipos[1]<vSim3->_height && Ipos[1]>=0 && idepth_var[i_idx]>100.0) // && image_var[i_idx]>100.0)
                {
                    // SET PointXYZ VERTEX
                    g2o::VertexSBAPointXYZ* vPoint = new g2o::VertexSBAPointXYZ();
                    vPoint->setEstimate(pts);
                    vPoint->setId(index);
                    vPoint->setFixed(true);
                    optimizer.addVertex(vPoint);
                    
                    // Set Image Edge
                    g2o::EdgeSim3ProjectXYZ* e_image = new g2o::EdgeSim3ProjectXYZ();
                    e_image->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(index)));
                    e_image->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));
                    if (Ipos2[0]<vSim3->_width && Ipos2[0]>=0 && Ipos2[1]<vSim3->_height && Ipos2[1]>=0)e_image->setMeasurement(ref[i_idx2]);
                    else e_image->setMeasurement(-1);
                    info << image_var[i_idx];
                    e_image->setInformation(info);
                    g2o::RobustKernelHuber* rk1 = new g2o::RobustKernelHuber;
                    rk1->setDelta(deltaHuber);
                    e_image->setRobustKernel(rk1);
                    optimizer.addEdge(e_image);
                    
                    // Set Depth Edge
                    g2o::EdgeSim3ProjectXYZD* e_depth = new g2o::EdgeSim3ProjectXYZD();
                    e_depth->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(index)));
                    e_depth->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));
                    e_depth->setMeasurement(1.0f);
                    info << idepth_var[i_idx];
                    e_depth->setInformation(info);
                    g2o::RobustKernelHuber* rk2 = new g2o::RobustKernelHuber;
                    rk2->setDelta(deltaHuber);
                    e_depth->setRobustKernel(rk2);
                    optimizer.addEdge(e_depth);
                    
                    index++;

                }
        }

    }     

    cout<<index<<endl;
    if(index>10)
    {
    optimizer.initializeOptimization();
    optimizer.computeActiveErrors();
    int g2oresult = optimizer.optimize(100);

    g2o::VertexSim3Expmap* vSim3_recov = static_cast<g2o::VertexSim3Expmap*>(optimizer.vertex(0));
    g2o::Sim3 g2oS12 = vSim3_recov->estimate();      

    Matrix4d result_mat = Sim3toMat(g2oS12);
    return result_mat;
    }
    else return Matrix4d::Identity();

}

void CamLocalization::debugImage(cv::Mat& depth_image,cv::Mat& dgx_image,cv::Mat& dgy_image, const float* depth_info)
{
    //generate lidar intensity image 
    int numpts = velo_cloud->points.size();    
    cv::Mat frame = cv::Mat::zeros(cv::Size(left_image.cols, left_image.rows), CV_8UC1);
    cv::Mat depth = cv::Mat::zeros(cv::Size(left_image.cols, left_image.rows), CV_32FC1);
    int u, v;
    for(size_t i=0; i<numpts;i++)
    {
        if (velo_cloud->points[i].z>-1.0){
            // Set frame and depth image
            Vector3d point(velo_cloud->points[i].x, velo_cloud->points[i].y, velo_cloud->points[i].z);
            Vector2d uv = ProjectTo2D(point);
            u = (int) uv(0);
            v = (int) uv(1);
            
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
            u = (int) uv(0);
            v = (int) uv(1);
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


