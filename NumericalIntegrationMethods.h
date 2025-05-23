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

//Eigen::MatrixXd ForwardEuler(std::vector<double>& time_vector, Eigen::MatrixXd solution_matrix, std::map<std::string, double>& vehicle_model, double step_size, std::map<std::string, std::vector<double> >& atmosphere_mod) {
//
//
//	for (int i = 1; i < time_vector.size(); i++) {
//
//		solution_matrix.col(i) = solution_matrix.col(i - 1) + step_size * flat_earth_EOM(time_vector[i - 1], solution_matrix.col(i - 1), vehicle_model, atmosphere_mod);
//	}
//
//
//	return solution_matrix;
//}

//Eigen::MatrixXd AdamsBashforth2(std::vector<double>& time_vector, Eigen::MatrixXd solution_matrix, std::map<std::string, double>& vehicle_model, double step_size, std::map<std::string, std::vector<double> >& atmosphere_mod)
//{
//	for (int i = 1; i < time_vector.size(); i++)
//	{
//		if (i == 1) {
//			solution_matrix.col(i) = solution_matrix.col(i - 1) + step_size * flat_earth_EOM(time_vector[i - 1], solution_matrix.col(i - 1), vehicle_model, atmosphere_mod);
//		}
//		else {
//
//			Eigen::VectorXd fim1 = flat_earth_EOM(time_vector[i - 1], solution_matrix.col(i - 1), vehicle_model, atmosphere_mod);
//			Eigen::VectorXd fim2 = flat_earth_EOM(time_vector[i - 1], solution_matrix.col(i - 2), vehicle_model, atmosphere_mod);
//			solution_matrix.col(i) = solution_matrix.col(i - 1) + 1.5 * step_size * fim1 - 0.5 * step_size * fim2;
//
//
//		}
//	}
//	return solution_matrix;
//
//}

//Eigen::MatrixXd FourOrderRungeKutta(std::vector<double>& time_vector, Eigen::MatrixXd solution_matrix, std::map<std::string, double>& vehicle_model, double step_size, std::map<std::string, std::vector<double> >& atmosphere_mod)
//{
//	for (int i = 1; i < time_vector.size(); i++)
//	{
//
//		Eigen::VectorXd k1 = step_size * flat_earth_EOM(time_vector[i - 1], solution_matrix.col(i - 1), vehicle_model, atmosphere_mod);
//		Eigen::VectorXd k2 = step_size * flat_earth_EOM(time_vector[i - 1] + 0.5 * step_size, solution_matrix.col(i - 1) + 0.5 * k1, vehicle_model, atmosphere_mod);
//		Eigen::VectorXd k3 = step_size * flat_earth_EOM(time_vector[i - 1] + 0.5 * step_size, solution_matrix.col(i - 1) + 0.5 * k2, vehicle_model, atmosphere_mod);
//		Eigen::VectorXd k4 = step_size * flat_earth_EOM(time_vector[i - 1] + step_size, solution_matrix.col(i - 1) + k3, vehicle_model, atmosphere_mod);
//
//		Eigen::VectorXd increment = (1.0 / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
//
//		solution_matrix.col(i) = solution_matrix.col(i - 1) + increment;
//
//	}
//	return solution_matrix;
//}