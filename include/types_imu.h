#ifndef G2O_TYPES_IMU_H
#define G2O_TYPES_IMU_H

#include "Thirdparty/g2o/g2o/core/base_vertex.h"
#include "Thirdparty/g2o/g2o/core/base_multi_edge.h"
#include "Thirdparty/g2o/g2o/core/base_binary_edge.h"
#include "Thirdparty/g2o/g2o/core/base_unary_edge.h"
#include "Thirdparty/g2o/g2o/types/types_sba.h"
#include "ImuState.h"

namespace g2o {

  using namespace Eigen;
  using namespace std;

  class VertexImu : public BaseVertex<15, ImuState>
  {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    VertexImu();

    virtual bool read(std::istream& is);
    virtual bool write(std::ostream& os) const;

    virtual void setToOriginImpl() {
      _estimate = ImuState();
      _estimate_prev = ImuState();
    }

    virtual void oplusImpl(const double* update_)
    {
      Eigen::Map<Vector15d> update(const_cast<double*>(update_));

      setEstimatePrev(estimate());
      ImuState s(update);
//      s.check_print();
      setEstimate(estimate()*s);
    }

    Vector2d _principle_point;
    Vector2d _focal_length;
    double _width, _height;

    Vector2d cam_map(const Vector3d & v) const
    {
      Vector2d res;
      res[0] = v[0]*_focal_length[0]/v[2] + _principle_point[0];
      res[1] = v[1]*_focal_length[1]/v[2] + _principle_point[1];
      return res;
    }

    Vector3d cam_map_inv(const Vector2d & v, const double d) const
    {
      Vector3d res;
      res[0] = d/_focal_length[0]*(v[0]-_principle_point[0]);
      res[1] = d/_focal_length[1]*(v[1]-_principle_point[1]);
      res[2] = d;
      return res;
    }

    Matrix4d cTv;

    const float* Image;
    const float* ImageGx;
    const float* ImageGy;
    const float* ImageInfo;

    const float* Depth;
    const float* DepthGx;
    const float* DepthGy;
    const float* DepthInfo;

    float* occ_image;
    int* occ_idx;


  protected:
  };

  class EdgeImu : public BaseBinaryEdge<9, ImuState, VertexImu, VertexImu>
  {
    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      EdgeImu();

      void initOcclusionimg()
      {
      }
      void clearMeasurement()
      {
      }
      void computeError()
      {
      }
      int computeError2(int& return_idx)
      {
        const VertexImu* v1 = static_cast<const VertexImu*>(_vertices[0]);
        const VertexImu* v2 = static_cast<const VertexImu*>(_vertices[1]);
        ImuState currrent = v2->estimate();
        _error = currrent.compute_error(_measurement, v1->estimate());
//        cout<<_error<<endl;
//        _information<< 1000; 
        return_idx = -1; 
        return 1;
      }

      virtual void linearizeOplus();
      virtual bool read(std::istream& is);
      virtual bool write(std::ostream& os) const;

  };

  class EdgeImuProjectXYZD : public BaseMultiEdge<1, double>
  {
    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      EdgeImuProjectXYZD();

      void initOcclusionimg()
      {
        const VertexImu* v1 = static_cast<const VertexImu*>(_vertices[1]);
        for(size_t i=0; i<v1->_width*v1->_height;i++)
        {
            v1->occ_image[i] = -1.0f;
            v1->occ_idx[i] = -1;
        }
      }
      void clearMeasurement()
      {
        _error<< 0.0f;
        _measurement = 0.0f;
      }
      void computeError()
      {

      }

      int computeError2(int& return_idx)
      {
        const VertexImu* v0 = static_cast<const VertexImu*>(_vertices[0]);
        const VertexImu* v1 = static_cast<const VertexImu*>(_vertices[1]);
        const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[2]);
        Matrix3d cTv_R = v1->cTv.block<3,3>(0,0);
        Vector3d cTv_t = v1->cTv.block<3,1>(0,3);

        Vector3d map_trans = cTv_R*v1->estimate().map_inv(v2->estimate())+cTv_t;
        Vector2d Ipos( v1->cam_map(map_trans) );
        int idx = (int)(((int)Ipos[1])*v1->_width+((int)Ipos[0]));

        if (Ipos[0]>=v1->_width || Ipos[0]<0 || Ipos[1]>=v1->_height || Ipos[1]<0 )
        {
          _error<< 0.0f;
          _measurement = 0.0f;
          return_idx = -1;
          return 0;
        }
        else if(!std::isfinite(v1->Depth[idx]))
        {
          _error<< 0.0f;
          _measurement = 0.0f;
          return_idx = -1;
          return 0;
        }
        else if(!std::isfinite(v1->DepthGx[idx]) || !std::isfinite(v1->DepthGy[idx]))
        {
          _error<< 0.0f;
          _measurement = 0.0f;
          return_idx = -1;
          return 0;

        }
        else
        {
          Matrix<double, 1, 1> e1(v1->Depth[idx]);
          Matrix<double, 1, 1> obsz(map_trans[2]);
          _information<< v1->DepthInfo[idx];
          _error = obsz-e1;


          float error_abs = _error[0]>0.0?_error[0]:-_error[0]; 
          if(v1->occ_image[idx]>0.0f){
            if(error_abs>v1->occ_image[idx]){
              _error<< 0.0f;
              _measurement = 0.0f;
              return_idx = -1;
              return 0;            
            }
            else{
              v1->occ_image[idx] = error_abs;
              int in_idx = return_idx;
              return_idx = v1->occ_idx[idx];
              v1->occ_idx[idx] = in_idx;
              _measurement = 1.0f;
              return 0;
            }  
          }
          else{
              if(_error[0]>3.0f || _error[0]<-3.0f){// && _measurement == 0.0f)
                  _error<< 0.0f;
                  _measurement = 0.0f;
                  return_idx = -1; 
                  return 0;
              }
              else{
                v1->occ_image[idx] = error_abs;
                v1->occ_idx[idx] = return_idx;
                _measurement = 1.0f;
                return_idx = -1;   
                return 1;
              } 
          }
            
        }

      }


      virtual void linearizeOplus();
      virtual bool read(std::istream& is);
      virtual bool write(std::ostream& os) const;

  };

  class EdgeImuProjectXYZ : public BaseMultiEdge<1, double>
  {
    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      EdgeImuProjectXYZ();

      virtual bool read(std::istream& is);
      virtual bool write(std::ostream& os) const;

      void initOcclusionimg()
      {
      }
      void clearMeasurement()
      {
      }

      void computeError()
      {
      }
    
      int computeError2(int& return_idx)
      {
        const VertexImu* v0 = static_cast<const VertexImu*>(_vertices[0]);
        const VertexImu* v1 = static_cast<const VertexImu*>(_vertices[1]);
        const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[2]);
        Matrix3d cTv_R = v1->cTv.block<3,3>(0,0);
        Vector3d cTv_t = v1->cTv.block<3,1>(0,3);

        Vector2d Ipos( v0->cam_map( cTv_R*v0->estimate().map_inv(v2->estimate())+cTv_t) );
        int i_idx = (int)(((int)Ipos[1])*v0->_width+((int)Ipos[0]));

        Vector2d Jpos( v1->cam_map( cTv_R*v1->estimate().map_inv(v2->estimate())+cTv_t) );
        int j_idx = (int)(((int)Jpos[1])*v1->_width+((int)Jpos[0]));

        ImuState imu_i = v0->estimate();
        if(!std::isfinite(Ipos[0])||!std::isfinite(Ipos[1])||!std::isfinite(Jpos[0])||!std::isfinite(Jpos[1]))
        {
          _error<< 0.0f;
          _measurement = 0.0f;
          return_idx = -1;
          return 0;            
        }
        else if (Ipos[0]>=v0->_width || Ipos[0]<0 || Ipos[1]>=v0->_height || Ipos[1]<0 )
        {
          _error<< 0.0f;
          _measurement = 0.0f;
          return_idx = -1;
          return 0;
        }
        else if (Jpos[0]>=v1->_width || Jpos[0]<0 || Jpos[1]>=v1->_height || Jpos[1]<0 )
        {
          _error<< 0.0f;
          _measurement = 0.0f;
          return_idx = -1;
          return 0;
        }
        else if(!std::isfinite(v0->Image[i_idx]))
        {
          _error<< 0.0f;
          _measurement = 0.0f;
          return_idx = -1;
          return 0;
        }
        else if(!std::isfinite(v1->Image[j_idx]))
        {
          _error<< 0.0f;
          _measurement = 0.0f;
          return_idx = -1;
          return 0;
        }
        else if(!std::isfinite(v0->ImageGx[i_idx]) || !std::isfinite(v0->ImageGy[i_idx]))
        {
          _error<< 0.0f;
          _measurement = 0.0f;
          return_idx = -1;
          return 0;
        }
        else if(!std::isfinite(v1->ImageGx[j_idx]) || !std::isfinite(v1->ImageGy[j_idx]))
        {
          _error<< 0.0f;
          _measurement = 0.0f;
          return_idx = -1;
          return 0;
        }
        else if(_measurement<0)
        {
          _error<< 0.0f;
          _measurement = 0.0f;
          return_idx = -1;
          return 0;
        }
        else
        {
          Matrix<double, 1, 1> e1(v0->Image[i_idx]);
          Matrix<double, 1, 1> obsz(v1->Image[j_idx]);
          _information<< v0->ImageInfo[i_idx];// + v1->ImageInfo[j_idx];//1000;// 
          _error = e1-obsz;
          _measurement = 1.0f;
          return_idx = -1;
          return 1;
        }

      }

      virtual void linearizeOplus();

  };

  class EdgeImuPhotometric : public BaseBinaryEdge<1, Vector2d, VertexImu, VertexImu>
  {
    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      EdgeImuPhotometric();

      virtual bool read(std::istream& is);
      virtual bool write(std::ostream& os) const;

      void initOcclusionimg()
      {
      }
      void clearMeasurement()
      {
      }

      void computeError()
      {
      }
    
      int computeError2(int& return_idx)
      {
        const VertexImu* v0 = static_cast<const VertexImu*>(_vertices[0]);
        const VertexImu* v1 = static_cast<const VertexImu*>(_vertices[1]);
        Matrix3d cTv_R = v1->cTv.block<3,3>(0,0);
        Vector3d cTv_t = v1->cTv.block<3,1>(0,3);

        Vector2d Ipos = _measurement;
        int i_idx = (int)(((int)Ipos[1])*v0->_width+((int)Ipos[0]));

        Vector3d xyz_i = v0->cam_map_inv(Ipos,v0->Depth[i_idx]);
        xyz_i = cTv_R.inverse()*(xyz_i-cTv_t);
        Vector3d xyz_j = v1->estimate().map_inv(v0->estimate().map(xyz_i));
        xyz_j = cTv_R*xyz_j+cTv_t;

        Vector2d Jpos( v1->cam_map(xyz_j) );
        int j_idx = (int)(((int)Jpos[1])*v1->_width+((int)Jpos[0]));

        ImuState imu_i = v0->estimate();
        if(!std::isfinite(Ipos[0])||!std::isfinite(Ipos[1])||!std::isfinite(Jpos[0])||!std::isfinite(Jpos[1]))
        {
          _error<< 0.0f;
          return_idx = -1;
          return 0;            
        }
        else if (Ipos[0]>=v0->_width || Ipos[0]<0 || Ipos[1]>=v0->_height || Ipos[1]<0 )
        {
          _error<< 0.0f;
          return_idx = -1;
          return 0;
        }
        else if (Jpos[0]>=v1->_width || Jpos[0]<0 || Jpos[1]>=v1->_height || Jpos[1]<0 )
        {
          _error<< 0.0f;
          return_idx = -1;
          return 0;
        }
        else if(!std::isfinite(v0->Image[i_idx]))
        {
          _error<< 0.0f;
          return_idx = -1;
          return 0;
        }
        else if(!std::isfinite(v1->Image[j_idx]))
        {
          _error<< 0.0f;
          return_idx = -1;
          return 0;
        }
        else if(!std::isfinite(v0->ImageGx[i_idx]) || !std::isfinite(v0->ImageGy[i_idx]))
        {
          _error<< 0.0f;
          return_idx = -1;
          return 0;
        }
        else if(!std::isfinite(v1->ImageGx[j_idx]) || !std::isfinite(v1->ImageGy[j_idx]))
        {
          _error<< 0.0f;
          return_idx = -1;
          return 0;
        }
        else
        {
          Matrix<double, 1, 1> e1(v0->Image[i_idx]);
          Matrix<double, 1, 1> obsz(v1->Image[j_idx]);
          _information<< v0->ImageInfo[i_idx];// + v1->ImageInfo[j_idx];//1000;// 
          _error = e1-obsz;
          return_idx = -1;
          return 1;
        }

      }

      virtual void linearizeOplus();

  };


}//end namespace



#endif // G2O_TYPES_IMU_H
