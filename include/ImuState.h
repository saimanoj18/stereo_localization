#ifndef G2O_IMUSTATE_H
#define G2O_IMUSTATE_H

#include <Eigen/Core>

#include "Thirdparty/g2o/g2o/types/se3_ops.h"

#ifndef GA
#define GA -9.8
#endif


namespace g2o {

    using namespace Eigen;

    typedef  Eigen::Matrix <double, 15, 1> Vector15d;

    class ImuState
    {
    public:
        ImuState():time_(0.0),R_(Matrix3d::Identity()),t_(Vector3d::Zero()),v_(Vector3d::Zero()),bias_g_(Vector3d::Zero()),bias_a_(Vector3d::Zero())
        {
        }

        ImuState(const Vector15d & update)
        {
            Vector3d omega;
            for (int i=0; i<3; i++){
              omega[i]=update[i];
              v_[i]=update[i+3];
              t_[i]=update[i+6];
              bias_g_[i]=update[i+9];
              bias_a_[i]=update[i+12]; 
            }

            R_ = Exp(omega);
            
        }
        
        void set_all(double time, Matrix3d R, Vector3d t, Vector3d v, Vector3d bias_g, Vector3d bias_a)
        {
            time_ = time;
            R_ = R;
            t_ = t;
            v_ = v;
            bias_g_ = bias_g; 
            bias_a_ = bias_a;
        }
        
        void set_time(double time){time_ = time;}
        void set_rotation(Matrix3d R){R_ = R;}
        void set_translation(Vector3d t){t_ = t;}
        void set_velocity(Vector3d v){v_ = v;}
        void set_biases(Vector3d bias_g, Vector3d bias_a){bias_g_ = bias_g; bias_a_ = bias_a;}

        void update_all(Vector3d w, Vector3d a, double delta_t)
        {
            Matrix3d R_i = R_;
            Vector3d v_i = v_;
            Vector3d gravity;
            gravity<<0,0,GA;

            double d2 = delta_t * delta_t;

            R_ = R_i * Exp(delta_t*(w-bias_g_));
            v_ = v_ + gravity*delta_t + R_i*(a-bias_a_)*delta_t;
            t_ = t_ + v_i*delta_t + 0.5*gravity*d2 + 0.5*R_i*(a-bias_a_)*d2;

        }

        ImuState operator *(const ImuState& other) const {
        ImuState ret;
        ret.R_ = R_*other.R_;
        ret.t_= t_ + (R_*other.t_);
        ret.v_ = v_ + other.v_;
        ret.bias_g_ = bias_g_ + other.bias_g_;
        ret.bias_a_ = bias_a_ + other.bias_a_;
        return ret;
        }
    //    void update_R(Matrix3d r){R_ = R_ * r;}
    //    void update_t(Vector3d t){t_ = t_ + R_*t;}
    //    void update_v(Vector3d v){v_ = v_ + v;}
    //    void update_bias_g(Vector3d bias_g){bias_g_ = bias_g_ + bias_g;}
    //    void update_bias_a(Vector3d bias_a){bias_a_ = bias_a_ + bias_a;}
        Vector3d map (const Vector3d& xyz) const {
            return (R_*xyz) + t_;
        }
     
    private:
        double time_;//sec

        Matrix3d R_;
        Vector3d t_;
        Vector3d v_;
        Vector3d bias_a_;
        Vector3d bias_g_;

        Matrix3d R_cov_;
        Matrix3d t_cov_;
        Matrix3d v_cov_;
        Matrix3d bias_a_cov_;
        Matrix3d bias_g_cov_;

     

    };
}
#endif // G2O_IMUSTATE_H
