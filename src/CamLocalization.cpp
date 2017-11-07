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


        int frame_half = (int) frameID/2;
        cout<<frame_half<<endl;
        GT_pose = GT_poses[frame_half];//ODO_pose;//
        cout<<"GT_pose:"<<GT_pose<<endl;
        pcl::transformPointCloud (*velo_cloud, *velo_raw, GT_pose);//transform to world coordinate
        mTfBr.sendTransform(tf::StampedTransform(wtb,ros::Time::now(), "/CamLoc/World", "/CamLoc/Camera"));
        MapPub.PublishMap(velo_raw,1);//publish velo raw
        MapPub.PublishPose(GT_pose,1);//publish GT pose


        //initialize 
        if(frameID == 0)CamLocInitialize(right_image);
        
        /////////////////////////disparity map generation/////////////////////////////        
        cv::Mat disp = cv::Mat::zeros(cv::Size(width, height), CV_16S);
        //disparity image estimate
        cv::Ptr<cv::StereoSGBM> sbm = cv::StereoSGBM::create(0,16*2,5);
        sbm->compute(left_image, right_image, disp);
        frameID = frameID+2;

        cout<<width<<", "<<height<<endl;
        /////////////////////////prepare optimization/////////////////////////////
        float* depth_gradientX = new float[width*height]();
        float* depth_gradientY = new float[width*height]();
        float* depth = new float[width*height]();
        float* depth_info = new float[width*height]();
        cv::Mat depth_image = cv::Mat(cv::Size(left_image.cols, left_image.rows), CV_32FC1);
        cv::Mat info_image = cv::Mat::zeros(cv::Size(width, height), CV_8UC3);

        cout<<"K: "<<K(0,0)<<endl;
        for(size_t i=0; i<width*height;i++)
        {
            int u, v;
            u = i%width;
            v = i/width;
            
            //depth
            int16_t d = disp.at<int16_t>(v,u);
            if(d==0 || d!=d || d<0.001) depth[i] = 1000;
            else depth[i] = 16*K(0,0)*0.54/((float)d);//0.54*K(0,0)/disp2.at<float>(v,u);
//            depth[i] = 10;
            depth_info[i] = 1/depth[i];

            if(depth[i]>1000)cout<<"error1: "<<depth[i]<<endl;
            if(depth[i]!=depth[i])cout<<"error2: "<<depth[i]<<endl;
 
            //depth image            
            depth_image.at<float>(v,u) = depth[i];

            //info image            
            info_image.at<cv::Vec3b>(v,u) = Compute_error_color(depth_info[i],1000);

        }
        cv::Mat dgx_image, dgy_image; // = cv::Mat(cv::Size(left_image.cols, left_image.rows), CV_8UC3);
        cv::Scharr(depth_image, dgx_image, CV_32FC1, 1, 0);
        cv::Scharr(depth_image, dgy_image, CV_32FC1, 0, 1);



        cv::Mat dgx_plot = cv::Mat::zeros(cv::Size(left_image.cols, left_image.rows), CV_8UC3);
        cv::Mat dgy_plot = cv::Mat::zeros(cv::Size(left_image.cols, left_image.rows), CV_8UC3); 
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

            if(depth_gradientX[i]!=depth_gradientX[i])depth_gradientX[i] = 0;
            if(depth_gradientY[i]!=depth_gradientY[i])depth_gradientY[i] = 0;

            if(depth_gradientX[i]>1||depth_gradientX[i]<-1)depth_gradientX[i] = 0;
            if(depth_gradientY[i]>1||depth_gradientY[i]<-1)depth_gradientY[i] = 0;


            
            //for visualization
            //image plot
            if(depth_gradientX[i]>0)dgx_plot.at<cv::Vec3b>(v,u)=cv::Vec3b(255*depth_gradientX[i],0,0);
            else if(depth_gradientX[i]==0)dgx_plot.at<cv::Vec3b>(v,u)=cv::Vec3b(0,0,0);
            else dgx_plot.at<cv::Vec3b>(v,u)=cv::Vec3b(0,0,-255*depth_gradientX[i]);
              
            if(depth_gradientY[i]>0)dgy_plot.at<cv::Vec3b>(v,u)=cv::Vec3b(255*depth_gradientY[i],0,0);
            else if(depth_gradientY[i]==0)dgy_plot.at<cv::Vec3b>(v,u)=cv::Vec3b(0,0,0);
            else dgy_plot.at<cv::Vec3b>(v,u)=cv::Vec3b(0,0,-255*depth_gradientY[i]);

            //cloud plot
            if(depth_info[i]>0){
                image_cloud->points[i].x = depth[i]/K(0,0)*(u-K(0,2));
                image_cloud->points[i].y = depth[i]/K(1,1)*(v-K(1,2)); 
                image_cloud->points[i].z = depth[i];
            }    
              
        }


        cv::imshow("dgx_plot", dgx_plot);
        cv::waitKey(3);
        cv::moveWindow("dgx_plot", 50,20);
        cv::imshow("dgy_plot", dgy_plot);
        cv::waitKey(3);
        cv::moveWindow("dgy_plot", 600,20);  
        cv::imshow("info_image", info_image);
        cv::waitKey(3);
        cv::moveWindow("info_image", 1000,20);
        cv::imshow("depth_image", depth_image);
        cv::waitKey(3);
        cv::moveWindow("depth_image", 1000,70);

        
        //Initialize EST_pose, velo_cloud
        EST_pose = EST_pose*update_pose;
        pcl::transformPointCloud (*velo_raw, *velo_cloud, EST_pose.inverse());
        cout<<"EST_pose"<<EST_pose<<endl;
        MapPub.PublishPose(ODO_pose,2);
        pcl::PointCloud<pcl::PointXYZ>::Ptr image_cloud_pre (new pcl::PointCloud<pcl::PointXYZ>); //publish previously estimated status
        pcl::transformPointCloud (*image_cloud, *image_cloud_pre, EST_pose);


        if(frameID>2){


            optimized_T = Matrix4f::Identity();
            
            //optimization
            optimized_T = Optimization(depth,depth_info,depth_gradientX,depth_gradientY);
            pcl::transformPointCloud (*velo_cloud, *velo_cloud, optimized_T);
            EST_pose = EST_pose*optimized_T.inverse();

        }


        cout<<"EST_pose"<<EST_pose<<endl;
//        cout<<"update_pose: "<<update_pose<<endl;
        MapPub.PublishPose(EST_pose,3);
        pcl::transformPointCloud (*image_cloud, *image_cloud, EST_pose);
        MapPub.PublishMap(image_cloud_pre,2); //image cloud before optimization 
        MapPub.PublishMap(image_cloud,3);

        //save poses 
        write_poses("EST_poses.txt", EST_pose);
        write_poses("GT_poses.txt", GT_pose);
        write_poses("ODO_poses.txt", ODO_pose);

        
        mTfBr.sendTransform(tf::StampedTransform(wtb,ros::Time::now(), "/CamLoc/World", "/CamLoc/Camera"));

        Velo_received = false;
        Left_received = false;
        Right_received = false;
        delete [] depth_gradientX;
        delete [] depth_gradientY;
        delete [] depth;
        delete [] depth_info;

//        cout<<"Elapsed time: %"<<end_time - start_time<<" usecs\n"<<endl;
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
                    Matrix4f prev_pose = ODO_pose;
                    pcl_ros::transformAsMatrix (wtb, ODO_pose);
                    update_pose = prev_pose.inverse()*ODO_pose;
//                    update_pose(2,3) = update_pose(2,3) -0.1f;         
                    ODO_time=wtb.stamp_;

                    cout<<"Velodyne input: "<<velo_cloud->points.size()<<endl;
                    pcl::transformPointCloud (*velo_cloud, *velo_cloud, cTv);//transform to camera keyframe coordinate
                    //pcl::transformPointCloud (*velo_cloud, *velo_cloud, ODO_pose);//transform to world coordinate
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

    cout<<"K Matrix: "<<K<<endl;
}


Matrix4f CamLocalization::Optimization(const float* idepth)
{
    //ICP optimization 
    Matrix<double,3,1> Sum_x;
    Matrix<double,3,1> Sum_x2;
    int numpts = velo_cloud->points.size();
    int index = 0;
    cout<<index<<endl;
    for(size_t i=0; i<numpts;i++)
    {
        // point sums
        Matrix<double,3,1> pts;
        pts <<velo_cloud->points[i].x, velo_cloud->points[i].y, velo_cloud->points[i].z ;

        Vector2d Ipos = ProjectTo2D(pts);
        int i_idx = ((int)Ipos[1])*width+((int)Ipos[0]);

        if ( pts[2]>1.0f && pts[2]<10.0f && pts[1]<0.0f)
        { 
            Matrix<double,3,1> pts2 = ReprojectTo3D(Ipos[0],Ipos[1],idepth[i_idx]);
            double dif = pts[2]-pts2[2];
            if (Ipos[0]<width && Ipos[0]>=0 && Ipos[1]<height && Ipos[1]>=0)// && dif<0.5f && dif>-0.5f)
            {
//                cout<<dif<<endl;
                Sum_x+=pts;
                Sum_x2+=pts2;
                index++;
                
            }
        }

    }
    cout<<index<<endl;
//    //solve Ax=b
//    Matrix<double,3,6> A;
//    Matrix<double,3,1> b;
//    Matrix3d A1, A2;
//    A1 <<0,-Sum_x2[2],Sum_x2[1],Sum_x2[2],0,-Sum_x2[0],-Sum_x2[1],Sum_x2[0],0;
//    A2 = ((double)index)*Matrix3d::Identity();
//    A.block<3,3>(0,0) = A1;
//    A.block<3,3>(0,3) = A2;
//    b = Sum_x-Sum_x2;
//    Matrix<double,6,1> sol;
//    sol = (A.transpose()*A).inverse()*A.transpose()*b;

    //build result_mat
    Matrix4f result_mat = Matrix4f::Identity();
//    result_mat << 1,0,0,sol[3],0,1,0,sol[4],0,0,1,sol[5],0,0,0,1;
//    cout<<"result mat"<<sol<<endl;
//    if(abs(sol[0])>10)return result_mat;
//    else return Matrix4f::Identity();
    if(index>0){
    Matrix<double,3,1> sol;
    sol = Sum_x-Sum_x2;
    result_mat(0,3) = 0;//sol[0]/((double)index);
    result_mat(1,3) = 0;//sol[1]/((double)index);
    result_mat(2,3) = -sol[2]/((double)index);
    }
    return result_mat;
    
}


Matrix4f CamLocalization::Optimization(const float* idepth, const float* idepth_var, const float* d_gradientX, const float* d_gradientY)
{
    //g2o optimization 
    cout<<"g2o Optimization"<<endl;
    const float deltaHuber = sqrt(1000);//1000 may be the best choice

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
    float* occ_container = new float[width*height]();
    vSim3->occ_image = occ_container;
    int* occ_idx = new int[width*height]();
    vSim3->occ_idx = occ_idx;
    optimizer.addVertex(vSim3);

    //Set map point vertices
    std::vector<Matrix<double,3,1> > pts;
    int numpts = velo_cloud->points.size();
    pts.reserve(numpts);

    vector<g2o::EdgeSim3ProjectXYZD*> vpEdges12;
    vector<size_t> vnIndexEdge;
    vnIndexEdge.reserve(2*numpts);
    vpEdges12.reserve(2*numpts);

    Matrix<double, 1, 1> info;
    info << 0.0f;

    int index = 1;
    for(size_t i=0; i<numpts;i++)
    {
        // Set map points
        pts[i] <<velo_cloud->points[i].x, velo_cloud->points[i].y, velo_cloud->points[i].z ;

        Vector2d Ipos( vSim3->cam_map(vSim3->estimate().map(pts[i])) );
        int i_idx = ((int)Ipos[1])*vSim3->_width+((int)Ipos[0]);
        
        
        if ( pts[i][2]>0.0f ){// && pts[i][2]<10.0f){
//            if(Ipos[0]<vSim3->_width && Ipos[0]>=0 && Ipos[1]<vSim3->_height && Ipos[1]>=0)
//            {
                
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

                    
                    vpEdges12.push_back(e01);
                    vnIndexEdge.push_back(index);

                    index++;
                }
//                else
//                {
//                    // SET PointXYZ VERTEX
//                    g2o::VertexSBAPointXYZ* vPoint = new g2o::VertexSBAPointXYZ();
//                    vPoint->setEstimate(pts[i]);
//                    vPoint->setId(index);
//                    vPoint->setFixed(true);
//                    optimizer.addVertex(vPoint);
//                    
//                    // Set Edges
//                    g2o::EdgeSim3ProjectXYZD* e01 = new g2o::EdgeSim3ProjectXYZD();
//                    e01->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(index)));
//                    e01->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));
//                    e01->setMeasurement(0.0f);//(e1-e2);//(pts[i][2]);//
//                    e01->setInformation(info);
//                    g2o::RobustKernelHuber* rk1 = new g2o::RobustKernelHuber;
//                    rk1->setDelta(deltaHuber);
//                    e01->setRobustKernel(rk1);

//                    optimizer.addEdge(e01);

//                    vpEdges12.push_back(e01);
//                    vnIndexEdge.push_back(index);

//                    index++;
//                }
//            }

        }

    }
    
    optimizer.initializeOptimization();
    optimizer.computeActiveErrors();
    cout << "Inlier nums: " << index << endl;

    optimizer.setVerbose(true);
    int g2oresult = optimizer.optimize(100);

    // Recover optimized Sim3
    g2o::VertexSim3Expmap* vSim3_recov = static_cast<g2o::VertexSim3Expmap*>(optimizer.vertex(0));
    g2o::Sim3 g2oS12 = vSim3_recov->estimate();      

    //debug image
    cv::Mat error_image = cv::Mat::zeros(cv::Size(width, height), CV_8UC3);
    cv::Mat error_image_after = cv::Mat::zeros(cv::Size(width, height), CV_8UC3);
    float error1 = 0;
    int sum1 = 0;
    float error2 = 0;
    int sum2 = 0;
    for(size_t i=0; i<numpts;i++)
    {

        Vector2d Ipos( vSim3->cam_map(pts[i]) );
        int i_idx = ((int)Ipos[1])*vSim3->_width+((int)Ipos[0]);

        Matrix<double,3,1> pts_transformed = vSim3->estimate().map(pts[i]); 
        Vector2d Ipos_trans( vSim3->cam_map(pts_transformed) );
        int i_idx_trans = ((int)Ipos_trans[1])*vSim3->_width+((int)Ipos_trans[0]);

        if (Ipos[0]<vSim3->_width && Ipos[0]>=0 && Ipos[1]<vSim3->_height && Ipos[1]>=0)
        {
            if ( pts[i][2]>0 && idepth_var[i_idx]>0) //&& idepth_var[i_idx]>0)// && pts[i][2]<20 )//&& idepth_var[i_idx]>0)
            {
                float derr = pts[i][2]-idepth[i_idx];
                
                if(derr<1.0f && derr>-1.0f)
                {
                cv::Vec3b residual_color;
                if(derr>0)residual_color = cv::Vec3b(0,0,derr*255);//Red
                else residual_color = cv::Vec3b(-derr*255,0,0);//Blue
                error_image.at<cv::Vec3b>(Ipos[1],Ipos[0]) = residual_color;//Compute_error_color(derr,10.0f);
                error1 += abs(derr);
                sum1++;
                }
            }
        }
        if (Ipos_trans[0]<vSim3->_width && Ipos_trans[0]>=0 && Ipos_trans[1]<vSim3->_height && Ipos_trans[1]>=0)
        {
            if ( pts[i][2]>0 && idepth_var[i_idx_trans]>0) // && idepth_var[i_idx]>0)// && pts[i][2]<20 )//&& idepth_var[i_idx]>0)
            {
                float derr = pts_transformed[2]-idepth[i_idx_trans];
                error_image_after.at<cv::Vec3b>(Ipos_trans[1],Ipos_trans[0]) = Compute_error_color(derr,10.0f);
                error2 += abs(derr);
                sum2++;
            }
        }
    }
    cv::imshow("error_image", error_image);
    cv::waitKey(3);
    cv::moveWindow("error_image", 50,300);
    cv::imshow("error_image_after", error_image_after);
    cv::waitKey(3);
    cv::moveWindow("error_image_after", 600,300); 

    cout<<"errors: "<<error1/sum1<<", "<<error2/sum2<<endl;
    cout<<"nums: "<<sum1<<", "<<sum2<<endl;


    // Check inliers
    int nBad=0;
    for(size_t i=0; i<vpEdges12.size();i++)
    {
        g2o::EdgeSim3ProjectXYZD* e12 = vpEdges12[i];
        if(!e12)
            continue;

        if(e12->chi2()>100)
        {
            size_t idx = vnIndexEdge[i];
            optimizer.removeEdge(e12);
            vpEdges12[i]=NULL;
            nBad++;
        }
    }
//    return Sim3toMat(g2oS12);

//    cout<<"Number of nBad is: "<<nBad<<", "<<vpEdges12.size()<<endl;
    cout<<"g2o success?: "<<g2oresult<<endl;
    Matrix4f result_mat = Sim3toMat(g2oS12);
//    result_mat(2,3) = -result_mat(2,3);

    delete [] occ_container;
    delete [] occ_idx;

    float res_x, res_y, res_z;
    res_x = result_mat(0,3)>0?result_mat(0,3):-result_mat(0,3);
    res_y = result_mat(1,3)>0?result_mat(1,3):-result_mat(1,3);
    res_z = result_mat(2,3)>0?result_mat(2,3):-result_mat(2,3);
//    result_mat(2,3) = -result_mat(2,3);
    if (res_x>0.3) result_mat(0,3) = 0;//result_mat(0,3)>0?0.1:-0.1;
    if (res_y>0.3) result_mat(1,3) = 0;//result_mat(1,3)>0?0.1:-0.1;
    if (res_z>1.0) result_mat(2,3) = 0;//result_mat(2,3)>0?0.1:-0.1;// || g2oresult<4 
    else return result_mat;


//    int nMoreIterations;
//    if(nBad>0){
//        nMoreIterations=10;
//        if(((float)vpEdges12.size()/nBad)<1.5f)//((max_size-nBad*2)<50)//((max_size/nBad)<2.1)//(nBad>500)//
//        {
//            cout<<"Number of nBad is large"<<endl;
//            nMoreIterations=-1;
//        }
//    }
//    else
//        nMoreIterations=5;



//    if(nMoreIterations>0){
//        cout<<nMoreIterations<<" more optimizations"<<endl;
//        // Optimize again only with inliers
//        optimizer.initializeOptimization();
//        optimizer.optimize(nMoreIterations);
//        // Recover optimized Sim3
//        g2o::VertexSim3Expmap* vSim3_recov = static_cast<g2o::VertexSim3Expmap*>(optimizer.vertex(0));
//        g2o::Sim3 g2oS12= vSim3_recov->estimate();
//        cout<<"optimized_T2:"<<Sim3toMat(g2oS12)<<endl;      
//        return Sim3toMat(g2oS12);
//    }
//    else{
//        return Matrix4f::Identity();
//    }
   
    

    
}
