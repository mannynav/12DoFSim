#pragma once

#include <Eigen/Dense>
#include <map>
#include <cmath>
#include <numbers>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

Eigen::VectorXd EulersKinematicalEq(double& p_b, double& q_b, double& r_b, double& phi_rad, double& theta_rad, double& psi_rad)
{
	Eigen::VectorXd result(3, 1);
	result.setZero();

	//Kinematic equations
	result[0] = p_b + sin(phi_rad) * tan(theta_rad) * q_b + cos(phi_rad) * tan(theta_rad) * r_b;

	result[1] = cos(phi_rad) * q_b - sin(phi_rad) * r_b;

	result[2] = sin(phi_rad) / cos(theta_rad) * q_b + cos(phi_rad) / cos(theta_rad) * r_b;

	return result;

}
Eigen::MatrixXd PoissonsKinematicalEq(double& p_b, double& q_b, double& r_b, Eigen::MatrixXd& DCM)
{

	Eigen::VectorXd omega(3, 3);
	omega.setZero();

	omega << 0, -r_b, q_b,
		r_b, 0, -p_b,
		-q_b, p_b, 0;

	Eigen::MatrixXd dDCM = omega * DCM;

	return dDCM;

}


Eigen::MatrixXd DCM(double phi, double theta, double psi)
{
	//Trig values from Euler angles
	double c_phi = cos(phi);
	double c_theta = cos(theta);
	double c_psi = cos(psi);
	double s_phi = sin(phi);
	double s_theta = sin(theta);
	double s_psi = sin(psi);
	double t_theta = tan(theta);

	Eigen::MatrixXd DCM(3, 3);
	DCM.setZero();

	//DCM coefficients
	double c00 = c_theta * c_psi;
	double c01 = c_theta * s_phi;
	double c02 = -s_theta;

	double c10 = c_psi * s_theta * s_phi - c_phi * s_psi;
	double c11 = c_phi * c_psi + s_theta * s_phi * s_psi;
	double c12 = c_theta * s_phi;

	double c20 = c_phi * c_psi * s_theta + s_phi * s_psi;
	double c21 = c_phi * s_theta * s_psi - c_psi * s_phi;
	double c22 = c_theta * c_phi;

	DCM << c00, c01, c02,
		c10, c11, c12,
		c20, c21, c22;

	return DCM;
}
Eigen::VectorXd quaternion(double p_b, double q_b, double r_b, Eigen::MatrixXd& DCM) {

	Eigen::VectorXd q(4, 1);
	q.setZero();

	Eigen::MatrixXd bigOmega(4, 4);
	bigOmega.setZero();

	bigOmega << 0, -p_b, -q_b, -r_b,
		p_b, 0, q_b, r_b,
		q_b, -r_b, 0, p_b,
		r_b, q_b, -p_b, 0;

	//std::cout << bigOmega << std::endl;

	double c11 = DCM(0, 0);
	double c22 = DCM(1, 1);
	double c33 = DCM(2, 2);

	double sTilda = sqrt(0.25 * (1 + c11 + c22 + c33));
	double qxTilda = sqrt(0.25 * (1 + c11 - c22 - c33));
	double qyTilda = sqrt(0.25 * (1 - c11 + c22 - c33));
	double qzTilda = sqrt(0.25 * (1 - c11 - c22 - c33));

	Eigen::VectorXd quatTilda(4, 1);
	quatTilda[0] = sTilda;
	quatTilda[1] = qxTilda;
	quatTilda[2] = qyTilda;
	quatTilda[3] = qzTilda;

	double qTildaMax = quatTilda.maxCoeff();

	if (qTildaMax == sTilda) {

		q[0] = sTilda;
		q[1] = ((DCM(1, 2) - DCM(2, 1)) / 4 * sTilda);
		q[2] = ((DCM(2, 0) - DCM(0, 2)) / 4 * sTilda);
		q[3] = ((DCM(0, 1) - DCM(1, 0)) / 4 * sTilda);

	}

	else if (qTildaMax == qxTilda) {

		q[0] = ((DCM(1, 2) - DCM(2, 1)) / 4 * qxTilda);
		q[1] = qxTilda;
		q[2] = ((DCM(0, 1) - DCM(1, 0)) / 4 * qxTilda);
		q[3] = ((DCM(2, 0) - DCM(0, 2)) / 4 * qxTilda);

	}

	else if (qTildaMax == qyTilda) {

		q[0] = ((DCM(2, 0) - DCM(0, 2)) / 4 * qyTilda);
		q[1] = ((DCM(0, 1) + DCM(1, 0)) / 4 * qyTilda);
		q[2] = qyTilda;
		q[3] = ((DCM(1, 2) - DCM(2, 1)) / 4 * qyTilda);

	}

	else if (qTildaMax == qzTilda) {

		q[0] = ((DCM(0, 1) - DCM(1, 0)) / 4 * qzTilda);
		q[1] = ((DCM(2, 0) + DCM(0, 2)) / 4 * qzTilda);
		q[2] = ((DCM(1, 2) + DCM(2, 1)) / 4 * qzTilda);
		q[3] = qzTilda;

	}


	Eigen::VectorXd result(4, 1);
	result = 0.5 * bigOmega * q;

	return result;

}
Eigen::VectorXd Quaternion2EulerAngle(Eigen::VectorXd& q) {

	//Compute DCM using quat values
	Eigen::MatrixXd DCM(3, 3);
	DCM.setZero();

	DCM(0, 0) = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
	DCM(0, 1) = 2 * (q[1] * q[2] + q[3] * q[0]);
	DCM(0, 2) = 2 * (q[1] * q[3] - q[2] * q[0]);


	DCM(1, 0) = 2 * (q[1] * q[2] - q[3] * q[0]);
	DCM(1, 1) = q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3];
	DCM(1, 2) = 2 * (q[2] * q[3] + q[1] * q[0]);

	DCM(2, 0) = 2 * (q[1] * q[3] + q[2] * q[0]);
	DCM(2, 1) = 2 * (q[2] * q[3] - q[1] * q[0]);
	DCM(2, 2) = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3];



	//Compute Euler angles from DCM
	double phi = atan2(DCM(1, 2), DCM(2, 1));
	double theta = asin(DCM(0, 2));
	double psi = atan2(DCM(0, 1), DCM(1, 0));

	Eigen::VectorXd result(3, 1);
	result[0] = phi;
	result[1] = theta;
	result[2] = psi;

	return result;

}
Eigen::VectorXd QuaternionKinematicalEqs(double& phi_rad, double& theta_rad, double& psi_rad, double& p_b, double& q_b, double& r_b)
{
	Eigen::MatrixXd DirCosineMatrix = DCM(phi_rad, theta_rad, psi_rad);
	Eigen::VectorXd quat = quaternion(p_b, q_b, r_b, DirCosineMatrix);

	return Quaternion2EulerAngle(quat);

}
Eigen::VectorXd DCM2angle(Eigen::MatrixXd& DCM)
{
	Eigen::VectorXd result(3, 1);
	result.setZero();

	//Compute Euler angles from DCM
	double phi = atan2(DCM(1, 2), DCM(2, 1));
	double theta = asin(DCM(0, 2));
	double psi = atan2(DCM(0, 1), DCM(1, 0));

	result[0] = phi;
	result[1] = theta;
	result[2] = psi;

	return result;

}