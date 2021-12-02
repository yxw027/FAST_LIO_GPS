// This is an advanced implementation of the algorithm described in the following paper:
//   J. Zhang and S. Singh. LOAM: Lidar Odometry and Mapping in Real-time.
//     Robotics: Science and Systems Conference (RSS). Berkeley, CA, July 2014. 

// Modifier: Tong Qin               qintonguav@gmail.com
// 	         Shaozu Cao 		    saozu.cao@connect.ust.hk


// Copyright 2013, Ji Zhang, Carnegie Mellon University
// Further contributions copyright (c) 2016, Southwest Research Institute
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from this
//    software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#pragma once

#include <cmath>

#include <pcl/point_types.h>

typedef pcl::PointXYZI PointType;
struct Pose6D {
  double x;
  double y;
  double z;
  double roll;
  double pitch;
  double yaw;
};
inline double rad2deg(double radians)
{
  return radians * 180.0 / M_PI;
}

inline double deg2rad(double degrees)
{
  return degrees * M_PI / 180.0;
}

inline Eigen::Matrix3d rotX(const double &degree)
{
  Eigen::Matrix3d R;
  double a = deg2rad(degree); 
  R<< 1,0,0, 
      0,cos(a),-sin(a), 
      0,sin(a),cos(a);  
  return R;
}
inline Eigen::Matrix3d rotY(const double &degree)
{
  Eigen::Matrix3d R;
  double a = deg2rad(degree); 
  R<< cos(a), 0, sin(a),
       0 , 1 ,0, 
      -sin(a) ,0 ,cos(a);  
  return R;
}

inline Eigen::Matrix3d rotZ(const double &degree)
{
  Eigen::Matrix3d R;
  double a = deg2rad(degree); 
  R<< cos(a), -sin(a), 0,
       sin(a),  cos(a), 0, 
       0,    0 ,  1;  
  return R;
}


inline void att2q(const double pitch, const double roll, const double yaw,
                   double *w, double *x, double *y, double *z)
{
            double roll2 = roll / 2;
            double pitch2 = pitch / 2;
            double yaw2 = yaw / 2;
            double sin_roll2 = sin(roll2);
            double sin_pitch2 = sin(pitch2);
            double sin_yaw2 = sin(yaw2);
            double cos_roll2 = cos(roll2);
            double cos_pitch2 = cos(pitch2);
            double cos_yaw2 = cos(yaw2);
            double sp = sin_pitch2, sr = sin_roll2, sy = sin_yaw2;
            double cp = cos_pitch2, cr = cos_roll2, cy = cos_yaw2;

            *w = cp * cr * cy - sp * sr * sy;
            *x = sp * cr * cy - cp * sr * sy;
            *y = cp * sr * cy + sp * cr * sy;
            *z = cp * cr * sy + sp * sr * cy;
}

inline Pose6D transform(Pose6D poseFrom,const Eigen::Quaterniond &QImu)
{
  Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
  Eigen::Matrix4d T_from = Eigen::Matrix4d::Identity();
  Eigen::Matrix3d Rwi = QImu.matrix();
  Eigen::Vector3d twi(0.82, 0, 1.60);   //lidar ->gps外参数
 
  double W = 1, X = 0, Y = 0, Z = 0;
  Eigen::Quaterniond Qtemp;
  // printf("poseFrom--pitch,roll,yaw= %f,%f,%f\n",poseFrom.pitch*180/M_PI,poseFrom.roll*180/M_PI,poseFrom.yaw*180/M_PI);
  att2q(poseFrom.pitch, poseFrom.roll, poseFrom.yaw, &W,&X, &Y, &Z);
  Qtemp.w() =W; Qtemp.x() =X;Qtemp.y() =Y;Qtemp.z() =Z;
  Eigen::Vector3d t_from(poseFrom.x, poseFrom.y, poseFrom.z);
  Eigen::Matrix3d r_from = Qtemp.matrix();

  T.block<3, 3>(0, 0) = Rwi;
  T.block<3, 1>(0, 3) = twi;
  T_from.block<3, 3>(0, 0) = r_from;
  T_from.block<3, 1>(0, 3) = t_from;

  Eigen::Matrix4d T_to = T*T_from;
  // Eigen::Matrix3d Rt = Ro * Rwi;
  Eigen::Quaterniond Qt(T_to.block<3, 3>(0, 0));
  double roll, pitch, yaw;
  tf2::Matrix3x3(tf2::Quaternion(Qt.x(), Qt.y(), Qt.z(), Qt.w())).getRPY(roll, pitch, yaw);
  Pose6D poseTo;
  poseTo.x = T_to(0,3);
  poseTo.y = T_to(1,3);
  poseTo.z = T_to(2,3);
  poseTo.roll = roll;
  poseTo.pitch = pitch;
  poseTo.yaw = yaw;
  // printf("poseTo---------->pitch,roll,yaw= %f,%f,%f\n",poseTo.pitch*180/M_PI,poseTo.roll*180/M_PI,poseTo.yaw*180/M_PI);
  return poseTo;

}



inline sensor_msgs::Imu transformRotZ(sensor_msgs::Imu  rawIMU,const double degree)
{       
        sensor_msgs::Imu  IMU_out;
        Eigen::Matrix3d R = rotZ(degree);
 
        Eigen::Quaterniond Qtemp(rawIMU.orientation.w,rawIMU.orientation.x,rawIMU.orientation.y,rawIMU.orientation.z);
        Eigen::Matrix3d R_from = Qtemp.matrix();
        Eigen::Matrix3d R_to =  R*R_from;

        Eigen::Quaterniond Q_to(R_to);
        IMU_out = rawIMU;
        IMU_out.orientation.w = Q_to.w();
        IMU_out.orientation.x = Q_to.x();
        IMU_out.orientation.y = Q_to.y();
        IMU_out.orientation.z = Q_to.z();

        return IMU_out;

}


inline Eigen::Vector3d  transform_gps2lidar( Eigen::Vector3d xyz,const sensor_msgs::Imu  curr_IMU )
{

  Eigen::Quaterniond  QImu;
  QImu.w() = curr_IMU.orientation.w;
  QImu.x() = curr_IMU.orientation.x;
  QImu.y() = curr_IMU.orientation.y;
  QImu.z() = curr_IMU.orientation.z;

  Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
  Eigen::Matrix4d T_from = Eigen::Matrix4d::Identity();
  Eigen::Matrix3d Rwi = QImu.matrix();
  Eigen::Vector3d twi(0.82, 0, 1.60);   //lidar ->gps外参数
 

  Eigen::Vector3d t_from(xyz[0], xyz[1], xyz[2]);
  Eigen::Matrix3d r_from = QImu.matrix();

  T.block<3, 1>(0, 3) = twi;
  T_from.block<3, 3>(0, 0) = r_from;
  T_from.block<3, 1>(0, 3) = t_from;

  Eigen::Matrix4d T_to = T*T_from;
  Eigen::Vector3d xyz_to;
  xyz_to[0] = T_to(0,3);
  xyz_to[1] = T_to(1,3);
  xyz_to[2] = T_to(2,3);
  
  return xyz_to;

}

