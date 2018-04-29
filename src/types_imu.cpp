#include "types_imu.h"

namespace g2o {

  VertexImu::VertexImu() : BaseVertex<15, ImuState>()
  {
  }

  EdgeImu::EdgeImu() : BaseBinaryEdge<9, ImuState, VertexImu, VertexImu>()
  {
  }

  EdgeImuProjectXYZD::EdgeImuProjectXYZD() : BaseMultiEdge<1, double>()
  {
    resize(3);
  }

  EdgeImuProjectXYZ::EdgeImuProjectXYZ() : BaseMultiEdge<1, double>()
  {
    resize(3);
  }


  bool VertexImu::read(std::istream& is)
  {
    return true;
  }
  bool VertexImu::write(std::ostream& os) const
  {
    return true;
  }

  bool EdgeImu::write(std::ostream& os) const
  {
    return true;
  }

  bool EdgeImu::read(std::istream& is)
  {
    return true;
  }

  bool EdgeImuProjectXYZ::write(std::ostream& os) const
  {
    return true;
  }

  bool EdgeImuProjectXYZ::read(std::istream& is)
  {
    return true;
  }

  bool EdgeImuProjectXYZD::write(std::ostream& os) const
  {
    return true;
  }

  bool EdgeImuProjectXYZD::read(std::istream& is)
  {
    return true;
  }

  void EdgeImu::linearizeOplus()
  {
    VertexImu* vi = static_cast<VertexImu*>(_vertices[0]);
    ImuState Imu_i = vi->estimate();
    VertexImu* vj = static_cast<VertexImu*>(_vertices[1]);
    ImuState Imu_j = vj->estimate();

    _jacobianOplusXi.block<3,15>(0,0) = Imu_j.Jacobian_Ri(_measurement, Imu_i);
    _jacobianOplusXi.block<3,15>(3,0) = Imu_j.Jacobian_vi(Imu_i);
    _jacobianOplusXi.block<3,15>(6,0) = Imu_j.Jacobian_ti(Imu_i);

    _jacobianOplusXj.block<3,15>(0,0) = Imu_j.Jacobian_Rj(_measurement, Imu_i);
    _jacobianOplusXj.block<3,15>(3,0) = Imu_j.Jacobian_vj(Imu_i);
    _jacobianOplusXj.block<3,15>(6,0) = Imu_j.Jacobian_tj(Imu_i);

//    cout<<_jacobianOplusXj<<endl;
  }


  void EdgeImuProjectXYZD::linearizeOplus()
  {
    _jacobianOplus[0].resize(1,15);
    _jacobianOplus[1].resize(1,15);
    _jacobianOplus[2].resize(1,3);

    VertexImu* vj = static_cast<VertexImu*>(_vertices[1]);
    ImuState Imu_j = vj->estimate();
    VertexSBAPointXYZ* v_point = static_cast<VertexSBAPointXYZ*>(_vertices[2]);
    Vector3d xyz = v_point->estimate();

    Matrix3d cTv_R = vj->cTv.block<3,3>(0,0);
    Vector3d cTv_t = vj->cTv.block<3,1>(0,3);

    Vector3d xyz_trans = Imu_j.map_inv(xyz);
    Vector3d xyz_trans2 = cTv_R*Imu_j.map_inv(xyz)+cTv_t;
    Vector2d Ipos(vj->cam_map(xyz_trans2));
    int idx = (int)(((int)Ipos[1])*vj->_width+((int)Ipos[0]));


    if (_measurement<0.1)
    {
        _jacobianOplus[0].setZero();
        _jacobianOplus[1].setZero();
        _jacobianOplus[2].setZero();
    }
    else
    {
//        cout<<_error<<endl;
        Matrix<double,1,2> D_u;
        Matrix<double,2,4> K_p;

        D_u(0,0) = vj->DepthGx[idx];
        D_u(0,1) = vj->DepthGy[idx];


        K_p(0,0) = vj->_focal_length1[0]/xyz_trans2[2];
        K_p(0,1) = 0;
        K_p(0,2) = -vj->_focal_length1[0]*xyz_trans2[0]/(xyz_trans2[2]*xyz_trans2[2]);
        K_p(0,3) = 0;//vj->_principle_point1[0];//
        K_p(1,0) = 0;
        K_p(1,1) = vj->_focal_length1[1]/xyz_trans2[2];
        K_p(1,2) = -vj->_focal_length1[1]*xyz_trans2[1]/(xyz_trans2[2]*xyz_trans2[2]);
        K_p(1,3) = 0;//vj->_principle_point1[1];//


        Matrix<double,4,6> Tp_note;
        Tp_note.block<3,3>(0,0) = -skew(xyz_trans);
        Tp_note.block<3,3>(0,3) = Matrix3d::Identity();
        Tp_note(3,0) = 0;
        Tp_note(3,1) = 0;
        Tp_note(3,2) = 0;
        Tp_note(3,3) = 0;
        Tp_note(3,4) = 0;
        Tp_note(3,5) = 0;

        Matrix<double,1,4> dm;
        dm<<0.0,0.0,1.0,0.0;
        Matrix<double,1,6> delta = -(dm-D_u*K_p)*vj->cTv*Tp_note;
        _jacobianOplus[0].setZero();
        _jacobianOplus[1].setZero();
        _jacobianOplus[1].block<1,3>(0,0) = delta.block<1,3>(0,0);        
        _jacobianOplus[1].block<1,3>(0,6) = delta.block<1,3>(0,3);
        _jacobianOplus[2].setZero();

//        cout<<_jacobianOplus[1]<<endl;
         
    }



  }



  void EdgeImuProjectXYZ::linearizeOplus()
  {
    _jacobianOplus[0].resize(1,15);
    _jacobianOplus[1].resize(1,15);
    _jacobianOplus[2].resize(1,3);

    VertexImu* vi = static_cast<VertexImu*>(_vertices[0]);
    ImuState Imu_i = vi->estimate();
    VertexImu* vj = static_cast<VertexImu*>(_vertices[1]);
    ImuState Imu_j = vj->estimate();
    VertexSBAPointXYZ* v_point = static_cast<VertexSBAPointXYZ*>(_vertices[2]);
    Vector3d xyz = v_point->estimate();

    Matrix3d cTv_R = vj->cTv.block<3,3>(0,0);
    Vector3d cTv_t = vj->cTv.block<3,1>(0,3);

    Vector2d Ipos(vi->cam_map(cTv_R*Imu_i.map_inv(xyz)+cTv_t));
    int i_idx = (int)(((int)Ipos[1])*vi->_width+((int)Ipos[0]));

    Vector2d Jpos(vj->cam_map(cTv_R*Imu_j.map_inv(xyz)+cTv_t));
    int j_idx = (int)(((int)Jpos[1])*vj->_width+((int)Jpos[0]));


    if(_measurement<0.1)
    {
        _jacobianOplus[0].setZero();
        _jacobianOplus[1].setZero();
        _jacobianOplus[2].setZero();
    }
    else
    {
//        cout<<_error<<endl;
        Matrix<double,1,2> imagei_u;        
        Matrix<double,1,2> imagej_u;
        Matrix<double,2,4> K_pi;
        Matrix<double,2,4> K_pj;
        Matrix<double,4,6> p_notei;
        Matrix<double,4,6> p_notej;

        imagei_u(0,0) = vi->ImageGx[i_idx];
        imagei_u(0,1) = vi->ImageGy[i_idx];
        imagej_u(0,0) = vj->ImageGx[j_idx];
        imagej_u(0,1) = vj->ImageGy[j_idx];

        Vector3d xyz_trans = Imu_i.map_inv(xyz);
        Vector3d xyz_trans2 = cTv_R*xyz_trans+cTv_t;        
        K_pi(0,0) = vi->_focal_length1[0]/xyz_trans2[2];
        K_pi(0,1) = 0;
        K_pi(0,2) = -vi->_focal_length1[0]*xyz_trans2[0]/(xyz_trans2[2]*xyz_trans2[2]);
        K_pi(0,3) = 0;
        K_pi(1,0) = 0;
        K_pi(1,1) = vi->_focal_length1[1]/xyz_trans2[2];
        K_pi(1,2) = -vi->_focal_length1[1]*xyz_trans2[1]/(xyz_trans2[2]*xyz_trans2[2]);
        K_pi(1,3) = 0;

        Matrix<double,4,6> Tp_notei;
        Tp_notei.block<3,3>(0,0) = -1*skew(xyz_trans);
        Tp_notei.block<3,3>(0,3) = Matrix3d::Identity();
        Tp_notei(3,0) = 0;
        Tp_notei(3,1) = 0;
        Tp_notei(3,2) = 0;
        Tp_notei(3,3) = 0;
        Tp_notei(3,4) = 0;
        Tp_notei(3,5) = 0;

        xyz_trans = Imu_j.map_inv(xyz);
        xyz_trans2 = cTv_R*xyz_trans+cTv_t;  
        K_pj(0,0) = vj->_focal_length1[0]/xyz_trans2[2];
        K_pj(0,1) = 0;
        K_pj(0,2) = -vj->_focal_length1[0]*xyz_trans2[0]/(xyz_trans2[2]*xyz_trans2[2]);
        K_pj(0,3) = 0;
        K_pj(1,0) = 0;
        K_pj(1,1) = vj->_focal_length1[1]/xyz_trans2[2];
        K_pj(1,2) = -vj->_focal_length1[1]*xyz_trans2[1]/(xyz_trans2[2]*xyz_trans2[2]);
        K_pj(1,3) = 0;

        Matrix<double,4,6> Tp_notej;
        Tp_notej.block<3,3>(0,0) = -1*skew(xyz_trans);
        Tp_notej.block<3,3>(0,3) = Matrix3d::Identity();
        Tp_notej(3,0) = 0;
        Tp_notej(3,1) = 0;
        Tp_notej(3,2) = 0;
        Tp_notej(3,3) = 0;
        Tp_notej(3,4) = 0;
        Tp_notej(3,5) = 0;

        Matrix<double,1,6> delta = -imagei_u*K_pi*vi->cTv*Tp_notei;
        _jacobianOplus[0].setZero();
        _jacobianOplus[0].block<1,3>(0,0) = delta.block<1,3>(0,0);        
        _jacobianOplus[0].block<1,3>(0,6) = delta.block<1,3>(0,3);
        delta = imagej_u*K_pj*vj->cTv*Tp_notej;
        _jacobianOplus[1].setZero();
        _jacobianOplus[1].block<1,3>(0,0) = delta.block<1,3>(0,0);        
        _jacobianOplus[1].block<1,3>(0,6) = delta.block<1,3>(0,3);
        _jacobianOplus[2].setZero();

//        cout<<_jacobianOplus[1]<<endl;
         
    }

  }


} // end namespace
