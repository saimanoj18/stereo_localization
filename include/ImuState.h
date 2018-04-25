#ifndef G2O_IMUSTATE_H
#define G2O_IMUSTATE_H

#include <Eigen/Core>

#include "Thirdparty/g2o/g2o/types/se3_ops.h"


namespace g2o {

    using namespace Eigen;

    typedef  Eigen::Matrix <double, 15, 1> Vector15d;

    class ImuState
    {
    public:
        ImuState():time_(0.0),R_(Matrix3d::Identity()),t_(Vector3d::Zero()),v_(Vector3d::Zero()),bias_g_(Vector3d::Zero()),bias_a_(Vector3d::Zero()),
        delta_ba_(Vector3d::Zero()),delta_bg_(Vector3d::Zero())
        {
            gravity<<0,0,-9.8;
        }

        ImuState(const Vector15d & update)
        {
            Vector3d omega;
            for (int i=0; i<3; i++){
              omega[i]=update[i];
              v_[i]=update[i+3];
              t_[i]=update[i+6];
              delta_bg_[i]=update[i+9];
              delta_ba_[i]=update[i+12]; 
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

        void update_all(Vector3d w, Vector3d a)
        {
            w_ = w;
            a_ = a;
        }

        void update_all(Vector3d w, Vector3d a, double delta_t, bool restart)
        {
            Matrix3d R_i = R_;
            Vector3d v_i = v_;
            double d2 = delta_t * delta_t;
            Matrix3d R_update = Exp(delta_t*(w_-bias_g_));

            //update time_
            time_ = time_ + delta_t;

            //update R, v, t
            R_ = R_i * R_update;
            v_ = v_ + gravity*delta_t + R_i*(a_-bias_a_)*delta_t;
            t_ = t_ + v_i*delta_t + 0.5*gravity*d2 + 0.5*R_i*(a_-bias_a_)*d2;

            //save R for bias jacobians
            for (std::vector<Matrix3d>::iterator it = R_for_b.begin() ; it != R_for_b.end(); ++it){
                *it = *it * R_update;
            }
            R_for_b.push_back(R_update);
            Jr_for_b.push_back(right_Jacobian(R_update)*delta_t);
            R_ij = R_ij * R_update;

            //update bias jacobians
            R_for_b.push_back(Matrix3d::Identity());
            R_bg = Matrix3d::Zero();
            for (std::pair<std::vector<Matrix3d>::iterator, std::vector<Matrix3d>::iterator> i(R_for_b.begin(), Jr_for_b.begin());i.first != R_for_b.end();++i.first, ++i.second){
                R_bg = R_bg - (*i.first) * (*i.second);
            }
            R_for_b.pop_back();
            v_ba = v_ba - R_ij*delta_t;
            v_bg = v_bg - R_ij*skew(a-bias_a_)*R_bg*delta_t;
            t_ba = t_ba + v_ba*delta_t - 0.5*R_ij*delta_t*delta_t;
            t_bg = t_bg + v_bg*delta_t - 0.5*R_ij*skew(a-bias_a_)*R_bg*delta_t*delta_t; 

            if(restart){

                //erase R for bias jacobians
                R_for_b.clear();
                Jr_for_b.clear();
                Jr_for_b.push_back(right_Jacobian(R_update)*delta_t);
                R_ij = R_update;
                delta_ba_ = Vector3d::Zero();
                delta_bg_ = Vector3d::Zero();
            }
    
            w_ = w;
            a_ = a;

        }

        ImuState operator *(const ImuState& other) const {
            ImuState ret;
            ret.R_ = R_*other.R_;
            ret.v_ = v_ + other.v_;
            ret.t_ = t_ + (R_*other.t_);
            ret.bias_g_ = bias_g_ + other.delta_bg_;
            ret.bias_a_ = bias_a_ + other.delta_ba_;
            ret.delta_bg_ = delta_bg_ + other.delta_bg_;
            ret.delta_ba_ = delta_ba_ + other.delta_ba_;

            ret.R_bg = R_bg;
            ret.v_ba = v_ba;        
            ret.v_bg = v_bg;
            ret.t_ba = t_ba;        
            ret.t_bg = t_bg;
            return ret;
        }

        //Rj - Ri
        ImuState operator -(const ImuState& other) const {
            ImuState ret;
            double time_ij = time_-other.time_;
            ret.R_ = other.R_.inverse()*R_;
            ret.v_ = other.R_*(v_-other.v_-gravity*time_ij);
            ret.t_ = other.R_*(t_-other.t_-other.v_*time_ij-0.5*gravity*time_ij*time_ij);
            return ret;
        }

        Matrix<double, 9, 1> compute_error(const ImuState& measurement, const ImuState& stateA)
        {
            Matrix<double, 9, 1> ret;
            double time_ij = time_-stateA.time_;
            Matrix3d exp1 = Exp(R_bg*delta_bg_);
            Matrix3d R_ab = stateA.R_.inverse()*R_;
            ret.block<3,1>(0,0) = Log(exp1.inverse()*measurement.R_.inverse()*R_ab);

            Vector3d v_ab = stateA.R_*(v_ - stateA.v_ - gravity*time_ij);
            Vector3d v_bgba = v_bg*delta_bg_ + v_ba*delta_ba_;
            ret.block<3,1>(3,0) = v_ab - measurement.v_ - v_bgba;

            Vector3d t_ab = stateA.R_*(t_ - stateA.t_ - stateA.v_*time_ij - 0.5*gravity*time_ij*time_ij);
            Vector3d t_bgba = t_bg*delta_bg_ + t_ba*delta_ba_;
            ret.block<3,1>(6,0) = t_ab - measurement.t_ - t_bgba;

            return ret;
        }

        Matrix<double,3,15> Jacobian_Ri(const ImuState& measurement, const ImuState& stateA)
        {
            Matrix<double,3,15> ret;
            ret.setZero();
            double time_ij = time_-stateA.time_;
            Matrix3d exp1 = Exp(R_bg*delta_bg_);
            Matrix3d R_ab = stateA.R_.inverse()*R_;
            Matrix3d res = exp1.inverse()*measurement.R_.inverse()*R_ab;
            
            ret.block<3,3>(0,0) = -left_Jacobian(res)*R_.inverse()*stateA.R_;
            return ret;
        }
        Matrix<double,3,15> Jacobian_Rj(const ImuState& measurement, const ImuState& stateA)
        {
            Matrix<double,3,15> ret;
            ret.setZero();
            double time_ij = time_-stateA.time_;
            Matrix3d exp1 = Exp(R_bg*delta_bg_);
            Matrix3d R_ab = stateA.R_.inverse()*R_;
            Matrix3d res = exp1.inverse()*measurement.R_.inverse()*R_ab;
            
            ret.block<3,3>(0,0) = left_Jacobian(res);
            ret.block<3,3>(0,9) = -left_Jacobian(res)*res.inverse()*right_Jacobian(Exp(R_bg*delta_bg_))*R_bg;
            return ret;
        }
        Matrix<double,3,15> Jacobian_vi(const ImuState& stateA)
        {
            Matrix<double,3,15> ret;
            ret.setZero();
            double time_ij = time_-stateA.time_;
            
            ret.block<3,3>(0,0) = skew(stateA.R_.inverse()*(v_-stateA.v_-gravity*time_ij));
            ret.block<3,3>(0,3) = -stateA.R_.inverse();
            return ret;
        }
        Matrix<double,3,15> Jacobian_vj(const ImuState& stateA)
        {
            Matrix<double,3,15> ret;
            ret.setZero();
            
            ret.block<3,3>(0,3) = stateA.R_.inverse();
            ret.block<3,3>(0,9) = -v_bg;
            ret.block<3,3>(0,12) = -v_ba;
            return ret;
        }
        Matrix<double,3,15> Jacobian_ti(const ImuState& stateA)
        {
            Matrix<double,3,15> ret;
            ret.setZero();
            double time_ij = time_-stateA.time_;
            
            ret.block<3,3>(0,0) = skew(stateA.R_.inverse()*(t_-stateA.t_-stateA.v_*time_ij-0.5*gravity*time_ij*time_ij));
            ret.block<3,3>(0,3) = -stateA.R_.inverse()*time_ij;
            ret.block<3,3>(0,6) = -Matrix3d::Identity();

            return ret;
        }
        Matrix<double,3,15> Jacobian_tj(const ImuState& stateA)
        {
            Matrix<double,3,15> ret;
            ret.setZero();
            double time_ij = time_-stateA.time_;
            
            ret.block<3,3>(0,6) = stateA.R_.inverse()*R_;
            ret.block<3,3>(0,9) = -t_bg;
            ret.block<3,3>(0,12) = -t_ba;
            return ret;
        }

        Vector3d map (const Vector3d& xyz) const {
            return (R_*xyz) + t_;
        }

        Vector3d map_inv (const Vector3d& xyz) const {
            return R_.inverse()*(xyz - t_);
        }
     
    private:
        double time_;//sec

        Matrix3d R_;
        Vector3d t_;
        Vector3d v_;
        Vector3d bias_a_;
        Vector3d bias_g_;
        
        Vector3d delta_ba_;
        Vector3d delta_bg_;

        Matrix3d R_bg;
        Matrix3d v_ba;        
        Matrix3d v_bg;
        Matrix3d t_ba;        
        Matrix3d t_bg;

        Matrix3d J_r;

        Vector3d w_;
        Vector3d a_;

        std::vector<Matrix3d> R_for_b;
        std::vector<Matrix3d> Jr_for_b;
        Matrix3d R_ij;
        

        Vector3d gravity;
        

//        Matrix3d R_cov_;
//        Matrix3d t_cov_;
//        Matrix3d v_cov_;
//        Matrix3d bias_a_cov_;
//        Matrix3d bias_g_cov_;

     

    };
}
#endif // G2O_IMUSTATE_H
