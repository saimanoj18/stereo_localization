// g2o - General Graph Optimization
// Copyright (C) 2011 H. Strasdat
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// Modified by Raúl Mur Artal (2014)
// - Added EdgeInverseSim3ProjectXYZ 
// - Modified VertexSim3Expmap to represent relative transformation between two cameras. Includes calibration of both cameras.

#ifndef G2O_SEVEN_DOF_EXPMAP_TYPES
#define G2O_SEVEN_DOF_EXPMAP_TYPES

#include "../core/base_vertex.h"
#include "../core/base_binary_edge.h"
#include "../core/base_unary_edge.h"
#include "types_six_dof_expmap.h"
#include "sim3.h"

namespace g2o {

  using namespace Eigen;

  /**
 * \brief Sim3 Vertex, (x,y,z,qw,qx,qy,qz)
 * the parameterization for the increments constructed is a 7d vector
 * (x,y,z,qx,qy,qz) (note that we leave out the w part of the quaternion.
 */
  class VertexSim3Expmap : public BaseVertex<7, Sim3>
  {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    VertexSim3Expmap();
    virtual bool read(std::istream& is);
    virtual bool write(std::ostream& os) const;

    virtual void setToOriginImpl() {
      _estimate = Sim3();
      _estimate_prev = Sim3();
    }

    virtual void oplusImpl(const double* update_)
    {
      Eigen::Map<Vector7d> update(const_cast<double*>(update_));

      if (_fix_scale)
        update[6] = 0;

//      std::cout<<update<<std::endl;
      setEstimatePrev(estimate());

      Sim3 s(update);
      setEstimate(s*estimate());
      //setEstimate(estimate()*s.inverse());
//      std::cout<<estimate()<<std::endl;
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

    Vector2d cam_map1(const Vector2d & v) const
    {
      Vector2d res;
      res[0] = v[0]*_focal_length1[0] + _principle_point1[0];
      res[1] = v[1]*_focal_length1[1] + _principle_point1[1];
      return res;
    }

    Vector2d cam_map2(const Vector2d & v) const
    {
      Vector2d res;
      res[0] = v[0]*_focal_length2[0] + _principle_point2[0];
      res[1] = v[1]*_focal_length2[1] + _principle_point2[1];
      return res;
    }
    const float* ImageD;
    const float* ImageGx;
    const float* ImageGy;
    const float* ImageInfo;

//    float* occ_image;
//    int* occ_idx;

    bool _fix_scale;


  protected:
  };

  /**
 * \brief 7D edge between two Vertex7
 */
  class EdgeSim3 : public BaseBinaryEdge<7, Sim3, VertexSim3Expmap, VertexSim3Expmap>
  {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    EdgeSim3();
    virtual bool read(std::istream& is);
    virtual bool write(std::ostream& os) const;
    void computeError()
    {
      const VertexSim3Expmap* v1 = static_cast<const VertexSim3Expmap*>(_vertices[0]);
      const VertexSim3Expmap* v2 = static_cast<const VertexSim3Expmap*>(_vertices[1]);

      Sim3 C(_measurement);
      Sim3 error_=C*v1->estimate()*v2->estimate().inverse();
      _error = error_.log();
    }

    virtual double initialEstimatePossible(const OptimizableGraph::VertexSet& , OptimizableGraph::Vertex* ) { return 1.;}
    virtual void initialEstimate(const OptimizableGraph::VertexSet& from, OptimizableGraph::Vertex* /*to*/)
    {
      VertexSim3Expmap* v1 = static_cast<VertexSim3Expmap*>(_vertices[0]);
      VertexSim3Expmap* v2 = static_cast<VertexSim3Expmap*>(_vertices[1]);
      if (from.count(v1) > 0)
  v2->setEstimate(measurement()*v1->estimate());
      else
  v1->setEstimate(measurement().inverse()*v2->estimate());
    }
  };

 /**
 * \brief 7D edge between pointxyz and sim3exp (YJKim)
 */
class EdgeSim3ProjectXYZD : public  BaseBinaryEdge<1, double, VertexSBAPointXYZ, VertexSim3Expmap>
{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    EdgeSim3ProjectXYZD();
    virtual bool read(std::istream& is);
    virtual bool write(std::ostream& os) const;
//    void initOcclusionimg()
//    {
//        const VertexSim3Expmap* v1 = static_cast<const VertexSim3Expmap*>(_vertices[1]);
////        const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[0]);
//        for(size_t i=0; i<v1->_width*v1->_height;i++)
//        {
//            v1->occ_image[i] = 0.0f;
//            v1->occ_idx[i] = -1;
//        }
//    }
//    void clearMeasurement()
//    {
//        _error<< 0.0f;
//        _measurement = 0.0f;
//    }
    void computeError()
    {
      const VertexSim3Expmap* v1 = static_cast<const VertexSim3Expmap*>(_vertices[1]);
      const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[0]);

      Vector2d Ipos( v1->cam_map(v1->estimate().map(v2->estimate())) );
      int idx = (int)(((int)Ipos[1])*v1->_width+((int)Ipos[0]));


      if (Ipos[0]>=v1->_width || Ipos[0]<0 || Ipos[1]>=v1->_height || Ipos[1]<0 )
      {
          _error<< 0.0f;
          _measurement = 0.0f;
      }
      else if(v1->ImageD[idx]<=0 || v1->ImageInfo[idx]<=0)
      {
          _error<< 0.0f;
          _measurement = 0.0f;
      }
      else
      {
          Matrix<double, 1, 1> e1(v1->ImageD[idx]);
          Matrix<double, 1, 1> obsz(v1->estimate().map(v2->estimate())[2]);
          _information<< v1->ImageInfo[idx];
          _error = (obsz-e1);
          _measurement = 1.0f;
      }

    }


    int computeError2(int& return_idx)
    {
      const VertexSim3Expmap* v1 = static_cast<const VertexSim3Expmap*>(_vertices[1]);
      const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[0]);

      Vector2d Ipos( v1->cam_map(v1->estimate().map(v2->estimate())) );
      int idx = (int)(((int)Ipos[1])*v1->_width+((int)Ipos[0]));

      if (Ipos[0]>=v1->_width || Ipos[0]<0 || Ipos[1]>=v1->_height || Ipos[1]<0 )
      {
          _error<< 0.0f;
          _measurement = 0.0f;
          return_idx = -1;
          return 0;
      }
      else if(!std::isfinite(v1->ImageD[idx]))
      {
          _error<< 0.0f;
          _measurement = 0.0f;
          return_idx = -1;
          return 0;
      }
     else if(!std::isfinite(v1->ImageGx[idx]) || !std::isfinite(v1->ImageGy[idx]))
     {
          _error<< 0.0f;
          _measurement = 0.0f;
          return_idx = -1;
          return 0;

     }
      else
      {
          Matrix<double, 1, 1> e1(v1->ImageD[idx]);
          Matrix<double, 1, 1> obsz(v1->estimate().map(v2->estimate())[2]);
          _information<< v1->ImageInfo[idx];
          _error = obsz-e1;

          _measurement = 1.0f;
          return_idx = -1; 
          return 1;

//          if(v1->occ_image[idx]>0.0f){
//            float error_abs = _error[0]>0?_error[0]:-_error[0]; 
//            if(error_abs>v1->occ_image[idx]){
//              _error<< 0.0f;
//              _measurement = 0.0f;
//              return_idx = -1;
//              return 0;            
//            }
//            else{
//              v1->occ_image[idx] = error_abs;
//              int in_idx = return_idx;
//              return_idx = v1->occ_idx[idx];
//              v1->occ_idx[idx] = in_idx;
//              _measurement = 1.0f;
//              return 0;
//            }
//          }  
            
      }

    }


    virtual void linearizeOplus();

};


/**/
class EdgeSim3ProjectXYZ : public  BaseBinaryEdge<1, double,  VertexSBAPointXYZ, VertexSim3Expmap>
{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    EdgeSim3ProjectXYZ();
    virtual bool read(std::istream& is);
    virtual bool write(std::ostream& os) const;

    void computeError()
    {
      const VertexSim3Expmap* v1 = static_cast<const VertexSim3Expmap*>(_vertices[1]);
      const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[0]);

//      Vector3d obs(_measurement);
//      _error = obs-v1->estimate().map(v2->estimate());
    }
    
    int computeError2(int& return_idx)
    {
      const VertexSim3Expmap* v1 = static_cast<const VertexSim3Expmap*>(_vertices[1]);
      const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[0]);

      const Vector3d v3 = v1->estimate().map(v2->estimate());
      
     
      Vector2d Ipos( v1->cam_map(v1->estimate().map(v2->estimate())) );
      int idx = (int)(((int)Ipos[1])*v1->_width+((int)Ipos[0]));

//      std::cout<<Ipos<<std::endl;
      if(!std::isfinite(Ipos[0])||!std::isfinite(Ipos[1]))
      {
          _error<< 0.0f;
//          _measurement = 0.0f;
          return_idx = -1;
          return 0;            
      }
      else if (Ipos[0]>=v1->_width || Ipos[0]<0 || Ipos[1]>=v1->_height || Ipos[1]<0 )
      {
          _error<< 0.0f;
//          _measurement = 0.0f;
          return_idx = -1;
          return 0;
      }
      else if(!std::isfinite(v1->ImageD[idx]))
      {
          _error<< 0.0f;
//          _measurement = 0.0f;
          return_idx = -1;
          return 0;
      }
     else if(!std::isfinite(v1->ImageGx[idx]) || !std::isfinite(v1->ImageGy[idx]))
     {
          _error<< 0.0f;
//          _measurement = 0.0f;
          return_idx = -1;
          return 0;

     }
      else
      {
          Matrix<double, 1, 1> e1(_measurement);
          Matrix<double, 1, 1> obsz(v1->ImageD[idx]);
          _information<< v1->ImageInfo[idx];//1000;// 
          _error = obsz-e1;
          
//          if(!std::isfinite(_information[0]))std::cout<<_error<<std::endl;
        

//          std::cout<<"obsz: "<<obsz<<std::endl;
//          std::cout<<"e1: "<<e1<<std::endl;  
//          std::cout<<_error<<std::endl;
          return_idx = -1; 
          return 1;

      }

    }

   virtual void linearizeOplus();

};

/**/
class EdgeInverseSim3ProjectXYZ : public  BaseBinaryEdge<3, Vector3d,  VertexSBAPointXYZ, VertexSim3Expmap>
{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    EdgeInverseSim3ProjectXYZ();
    virtual bool read(std::istream& is);
    virtual bool write(std::ostream& os) const;

    void computeError()
    {
      const VertexSim3Expmap* v1 = static_cast<const VertexSim3Expmap*>(_vertices[1]);
      const VertexSBAPointXYZ* v2 = static_cast<const VertexSBAPointXYZ*>(_vertices[0]);

      Vector3d obs(_measurement);
      _error = obs-v1->estimate().inverse().map(v2->estimate());
    }

   // virtual void linearizeOplus();

};

} // end namespace

#endif

