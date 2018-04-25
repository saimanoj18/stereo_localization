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

  Matrix3d skew(const Vector3d&v)
  {
    Matrix3d m;
    m.fill(0.);
    m(0,1)  = -v(2);
    m(0,2)  =  v(1);
    m(1,2)  = -v(0);
    m(1,0)  =  v(2);
    m(2,0) = -v(1);
    m(2,1) = v(0);
    return m;
  }

  Vector3d deltaR(const Matrix3d& R)
  {
    Vector3d v;
    v(0)=R(2,1)-R(1,2);
    v(1)=R(0,2)-R(2,0);
    v(2)=R(1,0)-R(0,1);
    return v;
  }

  Vector2d project(const Vector3d& v)
  {
    Vector2d res;
    res(0) = v(0)/v(2);
    res(1) = v(1)/v(2);
    return res;
  }

  Vector3d project(const Vector4d& v)
  {
    Vector3d res;
    res(0) = v(0)/v(3);
    res(1) = v(1)/v(3);
    res(2) = v(2)/v(3);
    return res;
  }

  Vector3d unproject(const Vector2d& v)
  {
    Vector3d res;
    res(0) = v(0);
    res(1) = v(1);
    res(2) = 1;
    return res;
  }

  Vector4d unproject(const Vector3d& v)
  {
    Vector4d res;
    res(0) = v(0);
    res(1) = v(1);
    res(2) = v(2);
    res(3) = 1;
    return res;
  }

  Matrix3d Exp(const Vector3d& omega)
  {
    double theta = omega.norm();
    Vector3d a = omega/theta;
    Matrix3d aa;
    aa <<a[0]*a[0],a[0]*a[1],a[0]*a[2],
         a[1]*a[0],a[1]*a[1],a[1]*a[2],
         a[2]*a[0],a[2]*a[1],a[2]*a[2];
    Matrix3d Omega = skew(omega);
    Matrix3d Omega2 = Omega*Omega;
    Matrix3d I = Matrix3d::Identity();
    Matrix3d R;

    double eps = 0.000001;
    if (theta<eps)
    {
        R = (I + Omega + 0.5*Omega2);
    }
    else
    {
        R = cos(theta)*I + sin(theta)*skew(a) + (1-cos(theta))*aa;
    }
    return R;
  }

  Vector3d Log(const Matrix3d& R)
  {
    Vector3d omega;
    double d = 0.5*(R(0,0)+R(1,1)+R(2,2)-1);
    Matrix3d Omega;

    double eps = 0.00001;
    Matrix3d I = Matrix3d::Identity();

    if (d>1-eps)
    {
        omega=0.5*deltaR(R);
    }
    else
    {
        double theta = acos(d);
        omega = theta/(2*sqrt(1-d*d))*deltaR(R);
    }
    return omega;        

  }

  Matrix3d right_Jacobian(const Matrix3d& R)
  {
    Vector3d xi = Log(R);
    Matrix3d Jr;
    Matrix3d I = Matrix3d::Identity();
    double theta = xi.norm();
    double theta2 = theta*theta;
    double theta3 = theta*theta*theta;
    Jr = I - (1-cos(theta))/theta2*skew(xi) + (theta - sin(theta))/theta3*skew(xi)*skew(xi);
    return Jr;
  }

  Matrix3d left_Jacobian(const Matrix3d& R)
  {
    Vector3d xi = Log(R);
    Matrix3d Jl;
    Matrix3d I = Matrix3d::Identity();
    double theta = xi.norm();
    double theta2 = theta*theta;
    Jl = I + 0.5*skew(xi) + (1.0/theta2 + (1+cos(theta))/(2*theta*sin(theta)) )*skew(xi)*skew(xi);
    return Jl;    

  }

