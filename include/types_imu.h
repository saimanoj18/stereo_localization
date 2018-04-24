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
      setEstimate(estimate()*s);
    }

    Vector2d _principle_point1, _principle_point2;
    Vector2d _focal_length1, _focal_length2;
    double _width, _height;

    Vector2d cam_map(const Vector3d & v) const
    {
      Vector2d res;
      res[0] = v[0]*_focal_length1[0]/v[2] + _principle_point1[0];
      res[1] = v[1]*_focal_length1[1]/v[2] + _principle_point1[1];
      return res;
    }

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

      void computeError()
      {
      }
      int computeError2(int& return_idx)
      {
        const VertexImu* v1 = static_cast<const VertexImu*>(_vertices[0]);
        const VertexImu* v2 = static_cast<const VertexImu*>(_vertices[1]);
//        Vector3d updated_g = v1->estimate().bias_g_ - _measurement.bias_g_;
//        Vector3d e1 = Log(Exp(_measurement.J_Rg*updated_g).transpose()*_measurement.R_.transpose()*v1->estimate().R_.inverse()*v2->estimate.R_);
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

      void initOcclusionimg()
      {
        const VertexImu* v2 = static_cast<const VertexImu*>(_vertices[1]);
        for(size_t i=0; i<v2->_width*v2->_height;i++)
        {
            v2->occ_image[i] = -1.0f;
            v2->occ_idx[i] = -1;
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

        Vector2d Ipos( v1->cam_map(v1->estimate().map(v2->estimate())) );
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
          Matrix<double, 1, 1> obsz(v1->estimate().map(v2->estimate())[2]);
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
              if(_error[0]>10.0f || _error[0]<-10.0f){// && _measurement == 0.0f)
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



}//end namespace



#endif // G2O_TYPES_IMU_H
