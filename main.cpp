

#include <Eigen/Dense>
#include <map>
#include <cmath>
#include <numbers>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "Input.h"
#include "Models.h"
#include "NumericalIntegrationMethods.h"
#include "MethodsForEulerAngles.h"
#include "Utility.h"
#include "Output.h"


Eigen::VectorXd TwelveDofSimulation(double time, Eigen::VectorXd state_vector_at_t, std::map<std::string, double> vehicle_model, std::map<std::string, std::vector<double> > atmosphere_mod)
{
	//return vector
	Eigen::VectorXd dx(12);

	double u = state_vector_at_t[0];
	double v = state_vector_at_t[1];
	double w = state_vector_at_t[2];
	double p = state_vector_at_t[3];
	double q = state_vector_at_t[4];
	double r = state_vector_at_t[5];
	double phi = state_vector_at_t[6];
	double theta = state_vector_at_t[7];
	double psi = state_vector_at_t[8];
	double p10_n = state_vector_at_t[9];
	double p20_n = state_vector_at_t[10];
	double p30_n = state_vector_at_t[11];

	//Trig values from Euler angles
	double c_phi = cos(phi);
	double c_theta = cos(theta);
	double c_psi = cos(psi);
	double s_phi = sin(phi);
	double s_theta = sin(theta);
	double s_psi = sin(psi);
	double t_theta = tan(theta);


	//Get mass and moment of inertia from aircraft model
	double m = vehicle_model["m_kg"];
	double Jxz_b = vehicle_model["Jxz_b"];
	double Jxx_b = vehicle_model["Jxx_b"];
	double Jyy_b = vehicle_model["Jyy_b"];
	double Jzz_b = vehicle_model["Jzz_b"];


	//Get reference dimensions
	double Aref = vehicle_model["Area_ref"];
	double b_m = vehicle_model["b_m"];
	double c_m = vehicle_model["c_m"];


	//get aerodynamic coefficients
	double Clp = vehicle_model["Clp"];
	double Clr = vehicle_model["Clr"];
	double Cmq = vehicle_model["Cmq"];
	double Cnp = vehicle_model["Cnp"];
	double Cnr = vehicle_model["Cnr"];


	//get current altitude
	double altitude = -p30_n;


	/////////////// Data from US Standard Atmosphere //////////////

	//values for rho
	//double rho_interp = 1.20;
	double rho_interp_2 = linearInterpolation(atmosphere_mod["alt_m"], atmosphere_mod["rho_kgpm3"], altitude);
	//std::cout << "rho interp: " << rho_interp_2 << std::endl;
	//double rho_interp = 1.225;   //later we will interpolate this

	//values for c_mps
	double c_interp = linearInterpolation(atmosphere_mod["alt_m"], atmosphere_mod["c_mps"], altitude);
	//double c_interp = 340.2939; //later we will interpolate this


	//Air data calculation(Mach, AoA, AoS)
	double translation_velocity = sqrt(u * u + v * v + w * w);
	double qbar = 0.5 * rho_interp_2 * translation_velocity * translation_velocity;
	double mach_number = translation_velocity / c_interp;
	double alpha_rad = atan2(w, u);
	double beta_rad = asin(v / translation_velocity);




	double w_over = 0.0000000001;
	if (u == 0 && w == 0) {
		w_over = 0.0;
	}

	else {
		w_over = w / u;
	}


	double v_over = 0.000000001;
	if (translation_velocity == 0 && v == 0) {
		v_over = 0;
	}
	else {
		v_over = v / translation_velocity;
	}


	double alpha = atan(w_over);
	double beta = asin(v_over);
	double s_alpha = sin(alpha);
	double c_alpha = cos(alpha);
	double s_beta = sin(beta);
	double c_beta = cos(beta);


	//gravity will act normal to earht tangent CS
	//double gs_n_mps2 = 9.81;
	//double gs_n_mps2 = 9.80665;
	double gs_interp = linearInterpolation(atmosphere_mod["alt_m"], atmosphere_mod["g_mps2"], altitude);

	//we need to transfrom gravity to body coordinate system
	double gx_b = -sin(theta) * gs_interp;
	double gy_b = sin(phi) * cos(theta) * gs_interp;
	double gz_b = cos(phi) * cos(theta) * gs_interp;


	//Aerodynamic forces
	//double drag = vehicle_model["CD_approx"] * qbar * vehicle_model["Aref"];
	double drag = 0.0;
	double side = 0.0;
	double lift = 0.0;

	double s_alpha_rad = sin(alpha_rad);
	double c_alpha_rad = cos(alpha_rad);
	double s_beta_rad = sin(beta_rad);
	double c_beta_rad = cos(beta_rad);


	//Components of DCM for body to stability axis
	double w2b_11 = c_alpha * c_beta;
	double w2b_12 = -c_alpha * s_beta;
	double w2b_13 = -s_alpha;
	double w2b_21 = s_beta;
	double w2b_22 = c_beta;
	double w2b_23 = 0.0;
	double w2b_31 = s_alpha * c_beta;
	double w2b_32 = -s_alpha * s_beta;
	double w2b_33 = c_alpha;

	//External forces
	double Fx_b = -(c_alpha * c_beta * drag - c_alpha * s_beta * side - s_alpha * lift);
	double Fy_b = -(s_beta * drag + c_beta * side + 0.0*lift);
	double Fz_b = -(s_alpha * c_beta * drag - s_alpha * s_beta * side + c_alpha * lift);

	//External moments
	double l_b = Clp * p * b_m / (2.0 * translation_velocity) + Clr * r * b_m / (2.0 * translation_velocity); //l_b is zero if b_m = 0 or Clp or Clr are 0.
	double m_b = Cmq * q * c_m / (2.0 * translation_velocity); //m_b is 0 if Cmq or q_b or c_m = 0.
	double n_b = Cnp * p * b_m / (2.0 * translation_velocity) + Cnr * r * b_m / (2.0 * translation_velocity); //n_b is 0 if Cnp or p_b or b_m = 0.


	//Denominator for roll and yaw rate equations
	double denominator = Jxx_b * Jzz_b - Jxz_b * Jxz_b;


	
	dx[0] = 1 / m * Fx_b + gx_b - w * q + v * r; //x-axis velocity eq
	dx[1] = 1 / m * Fy_b + gy_b - u * r + w * p; //y-axis velocity
	dx[2] = 1 / m * Fz_b + gz_b - v * p + u * q; //z-axis velocity


	//Roll eq
	dx[3] = (Jxz_b * (Jxx_b - Jyy_b + Jzz_b) * p * q - (Jzz_b * (Jzz_b - Jyy_b) + Jxz_b * Jxz_b) * q * r + Jzz_b * l_b + Jxz_b * n_b) / denominator;

	//Pitch eq
	dx[4] = ((Jzz_b - Jxx_b) * p * r - Jxz_b * (p * p - r * r) + m_b) / Jyy_b;

	//Yaw eq
	dx[5] = ((Jxx_b * (Jxx_b - Jyy_b) + Jxz_b * Jxz_b) * p * q - Jxz_b * (Jxx_b - Jyy_b + Jzz_b) * q * r + Jxz_b * l_b + Jxx_b * n_b) / denominator;

	//Kinematic equations
	Eigen::VectorXd EulerKinEq = EulersKinematicalEq(p, q, r, phi, theta, psi);
	dx[6] = EulerKinEq[0];
	dx[7] = EulerKinEq[1];
	dx[8] = EulerKinEq[2];


	//Eigen::VectorXd QuatKinEqs = QuaternionKinematicalEqs(phi, theta, psi, p, q, r);
	//dx[6] = QuatKinEqs[0];
	//dx[7] = QuatKinEqs[1];
	//dx[8] = QuatKinEqs[2];


	//Position (Navigation) equations
	dx[9] = c_theta * c_phi * u +(-c_phi * s_psi + s_phi * s_theta * c_psi) * v +  (s_phi * s_psi + c_phi * s_theta * c_psi) * w;
	dx[10] = c_phi * s_psi * u + (c_phi * c_psi + s_phi * s_theta * s_psi) * v +  (-s_phi * c_psi + c_phi * s_theta * s_psi) * w;
	dx[11] =      -s_theta * u +                           s_phi * c_theta * v +                             c_phi * c_theta * w;

	return dx;
}

Eigen::MatrixXd FourOrderRungeKutta(std::vector<double>& time_vector, Eigen::MatrixXd solution_matrix, std::map<std::string, double>& vehicle_model, double step_size, std::map<std::string, std::vector<double> >& atmosphere_mod)
{
	for (int i = 1; i < time_vector.size(); i++)
	{

		Eigen::VectorXd k1 = step_size * TwelveDofSimulation(time_vector[i - 1], solution_matrix.col(i - 1), vehicle_model, atmosphere_mod);
		Eigen::VectorXd k2 = step_size * TwelveDofSimulation(time_vector[i - 1] + 0.5 * step_size, solution_matrix.col(i - 1) + 0.5 * k1, vehicle_model, atmosphere_mod);
		Eigen::VectorXd k3 = step_size * TwelveDofSimulation(time_vector[i - 1] + 0.5 * step_size, solution_matrix.col(i - 1) + 0.5 * k2, vehicle_model, atmosphere_mod);
		Eigen::VectorXd k4 = step_size * TwelveDofSimulation(time_vector[i - 1] + step_size, solution_matrix.col(i - 1) + k3, vehicle_model, atmosphere_mod);

		Eigen::VectorXd increment = (1.0 / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

		solution_matrix.col(i) = solution_matrix.col(i - 1) + increment;

	}
	return solution_matrix;
}

int main() {

	long double pi = 3.141592653589793L;

	///////////////////////////////////////// Initialize variables for simulation ///////////////////////////////////////////////////

	//Simple conversions
	double radian2degrees = 180 / pi;
	double degrees2radians = 1 / radian2degrees;

	//// Initial conditions - expressed in the body frame - Euler angles in radians ////
	double u = 0.0000000001;
	double v = 0.0;
	double w = 0.0;
	double p = 10*degrees2radians;
	double q = 20 * degrees2radians;
	double r = 30 * degrees2radians;
	double phi = 0 * pi / 180.0;
	double theta = 0*degrees2radians;
	double psi = 0 * degrees2radians;
	double p10_n_m = 0;
	double p20_n_m = 0;
	double p30_n_m = -30000/3.28;
	Eigen::VectorXd initial_vector{ {u, v,w, p,q,r, phi,theta,psi,p10_n_m,p20_n_m,p30_n_m} };


	//// Time conditions ////
	double t_0 = 0.0;
	double tf_s = 30;
	double step = 0.01;

	std::vector<double> time_step_vector = generate_arange(t_0, tf_s + step, step); //Vector of points at which solution is approximated
	Eigen::MatrixXd solution_matrix(12, time_step_vector.size()); //Initilize matrix that will hold solutions
	solution_matrix.setZero().col(0) = initial_vector;

	///////////////////////////////////////////////// End of initializing variables ///////////////////////////////////////////////////////////////




	/////////////////////////////////// USSA 1976 - Read in the gravity and atmosphere data and we'll store it in a dictionary ///////////////////////////////

	std::string filePath_altitude = "gravity_atm_data_alt_m.csv"; //m
	std::string filePath_c_mps = "gravity_atm_data_c_mps.csv"; //mps
	std::string filePath_gravity = "gravity_atm_data_g_mps2.csv"; //mps2
	std::string filePath_air_density = "gravity_atm_data_rho_kgpm3.csv"; //kgpm3
	std::ifstream alt_m;

	alt_m.open(filePath_altitude);

	if (alt_m.fail()) {
		std::cout << "Failed to open file" << std::endl;
	}


	std::vector<double> altitude_vector = readNumericCSVColumn(filePath_altitude);
	std::vector<double> c_mps_vec = readNumericCSVColumn(filePath_c_mps);
	std::vector<double> gravity_vector = readNumericCSVColumn(filePath_gravity);
	std::vector<double> air_density_vector = readNumericCSVColumn(filePath_air_density);
	//////////// end of reading data for atmosphere ////////////////


	/// Store atmosphere data in dictionary
	std::map<std::string, std::vector<double>> atmosphere_map;
	atmosphere_map["alt_m"] = altitude_vector;
	atmosphere_map["c_mps"] = c_mps_vec;
	atmosphere_map["g_mps2"] = gravity_vector;
	atmosphere_map["rho_kgpm3"] = air_density_vector;
	
	//////////////////////////////////// End of gravity model extraction ////////////////////////////////////////////////////////




	////////////////////////////////////// Vehicles/Objects to simulate ////////////////////////////////////////////////////////

	//std::map<std::string, double> vehicle_map = NASA_Atmos01_Sphere();
	std::map<std::string, double> vehicle_map = NASA_Atmos02_Brick();
	std::cout << vehicle_map["Vterm"] << std::endl;

	////////////////////////////////////// End of vehicle selection ///////////////////////////////////////////////////////////




	////////////////////////////////////// Run simulation ////////////////////////////////////////////////////////////////////

	//Eigen::MatrixXd resultMatrix = ForwardEuler(time_step_vector, solution_matrix, vehicle_map, step, atmosphere_map);
	Eigen::MatrixXd resultMatrix = FourOrderRungeKutta(time_step_vector, solution_matrix, vehicle_map, step, atmosphere_map);
	//Eigen::MatrixXd resultMatrix = AdamsBashforth2(time_step_vector, solution_matrix, vehicle_map, step, atmosphere_map);

	///////////////////////////////////// End of simulation /////////////////////////////////////////////////////////////////



	///////////////////////////////////// Post processing of results ////////////////////////////////////////////////////////
	
	//Vectors to store air data specs
	Eigen::VectorXd Altitude(time_step_vector.size());
	Eigen::VectorXd InterpSpeedOfSound(time_step_vector.size());
	Eigen::VectorXd AirDensity(time_step_vector.size());
	Eigen::VectorXd TranslationalVelocity(time_step_vector.size());
	Eigen::VectorXd Mach(time_step_vector.size());


	Eigen::VectorXd AngleOfAttack(time_step_vector.size()); //radians
	Eigen::VectorXd AngleOfSideslip(time_step_vector.size()); //radians


	//Cosine and sine values for Euler angles from main simulation
	Eigen::VectorXd c_phi(time_step_vector.size());
	Eigen::VectorXd c_theta(time_step_vector.size());
	Eigen::VectorXd c_psi(time_step_vector.size());
	Eigen::VectorXd s_phi(time_step_vector.size());
	Eigen::VectorXd s_theta(time_step_vector.size());
	Eigen::VectorXd s_psi(time_step_vector.size());


	//Coordinate transformation vectors for body to NED frame - easier than having update a rotation matrix
	Eigen::VectorXd body2NED_11(time_step_vector.size());
	Eigen::VectorXd body2NED_12(time_step_vector.size());
	Eigen::VectorXd body2NED_13(time_step_vector.size());
	Eigen::VectorXd body2NED_21(time_step_vector.size());
	Eigen::VectorXd body2NED_22(time_step_vector.size());
	Eigen::VectorXd body2NED_23(time_step_vector.size());
	Eigen::VectorXd body2NED_31(time_step_vector.size());
	Eigen::VectorXd body2NED_32(time_step_vector.size());
	Eigen::VectorXd body2NED_33(time_step_vector.size());


	//u,v,w transformed to NED frame
	Eigen::VectorXd u_NED(time_step_vector.size());
	Eigen::VectorXd v_NED(time_step_vector.size());
	Eigen::VectorXd w_NED(time_step_vector.size());


	//phi, theta, psi in radians
	Eigen::VectorXd phi_rad(time_step_vector.size());
	Eigen::VectorXd theta_rad(time_step_vector.size());
	Eigen::VectorXd psi_rad(time_step_vector.size());


	//phi, theta, psi in degrees
	Eigen::VectorXd phi_deg(time_step_vector.size());
	Eigen::VectorXd theta_deg(time_step_vector.size());
	Eigen::VectorXd psi_deg(time_step_vector.size());


	Eigen::VectorXd roll_rate(time_step_vector.size());
	Eigen::VectorXd pitch_rate(time_step_vector.size());
	Eigen::VectorXd yaw_rate(time_step_vector.size());


	

	for (int i = 0; i < time_step_vector.size(); i++)
	{
		Altitude[i] = -resultMatrix(11, i);
		InterpSpeedOfSound[i] = linearInterpolation(altitude_vector, c_mps_vec, Altitude[i]);
		AirDensity[i] = linearInterpolation(altitude_vector, air_density_vector, Altitude[i]);
		TranslationalVelocity[i] = sqrt(resultMatrix(0, i) * resultMatrix(0, i) + resultMatrix(1, i) * resultMatrix(1, i) + resultMatrix(2, i) * resultMatrix(2, i));
		Mach[i] = TranslationalVelocity[i] / InterpSpeedOfSound[i];

		c_phi[i] = cos(resultMatrix(6, i));
		c_theta[i] = cos(resultMatrix(7, i));
		c_psi[i] = cos(resultMatrix(8, i));
		s_phi[i] = sin(resultMatrix(6, i));
		s_theta[i] = sin(resultMatrix(7, i));
		s_psi[i] = sin(resultMatrix(8, i));

		body2NED_11[i] = c_theta[i] * c_psi[i];
		body2NED_12[i] = -c_phi[i] * s_psi[i] + s_phi[i] * s_theta[i] * c_psi[i];
		body2NED_13[i] = s_phi[i] * s_psi[i] + c_phi[i] * s_theta[i] * c_psi[i];
		body2NED_21[i] = c_theta[i] * s_psi[i];
		body2NED_22[i] = c_phi[i] * c_psi[i] + s_phi[i] * s_theta[i]*s_psi[i];
		body2NED_23[i] = -s_phi[i] * c_psi[i] + c_phi[i] * s_theta[i] * s_psi[i];
		body2NED_31[i] = -s_theta[i];
		body2NED_32[i] = s_phi[i] * c_theta[i];
		body2NED_33[i] = c_phi[i] * c_theta[i];

		u_NED[i] = body2NED_11[i] * resultMatrix(0, i) + body2NED_12[i] * resultMatrix(1, i) + body2NED_13[i] * resultMatrix(2, i);
		v_NED[i] = body2NED_21[i] * resultMatrix(0, i) + body2NED_22[i] * resultMatrix(1, i) + body2NED_23[i] * resultMatrix(2, i);
		w_NED[i] = body2NED_31[i] * resultMatrix(0, i) + body2NED_32[i] * resultMatrix(1, i) + body2NED_33[i] * resultMatrix(2, i);

		phi_rad[i] = atan2(body2NED_32[i], body2NED_33[i]);
		theta_rad[i] = -asin(body2NED_31[i]);
		psi_rad[i] = atan2(body2NED_21[i],body2NED_11[i]);

		phi_deg[i] = (180 / 3.14) * phi_rad[i];
		theta_deg[i] = (180 / 3.14) * theta_rad[i];
		psi_deg[i] = (180 / 3.14) * psi_rad[i];

		roll_rate[i] = resultMatrix(3, i);
		pitch_rate[i] = resultMatrix(4, i);
		yaw_rate[i] = resultMatrix(5, i);

	}


	//Angle of attack
	for (int i = 0; i < time_step_vector.size(); i++)
	{
		double temp; //avoid division by zero

		if (resultMatrix(0,i) == 0) {
			temp = 0.0;
		}
		else {
			temp = resultMatrix(2,i)/resultMatrix(0,i);
		}
		AngleOfAttack[i] = atan(temp);
	}


	//Angle of sideslip
	for (int i = 0; i < time_step_vector.size(); i++)
	{
		double temp; //avoid division by zero

		if (TranslationalVelocity[i] == 0) {
			temp = 0.0;
		}
		else {
			temp = resultMatrix(1,i) / TranslationalVelocity[i];
		}

		AngleOfSideslip[i] = asin(temp);
	}




	//Vector for plotting purposes
	Eigen::VectorXd time_vector(time_step_vector.size());
	for (int i = 0; i < time_step_vector.size(); i++) {
		time_vector[i] = time_step_vector[i];
	}



	Eigen::MatrixXd postProcessMatrix(29, time_step_vector.size());
	postProcessMatrix.row(0) = time_vector;
	postProcessMatrix.row(1) = Altitude.transpose(); //m
	postProcessMatrix.row(2) = 3.28*Altitude.transpose(); // ft - mean sea level

	postProcessMatrix.row(3) = 3.28 * InterpSpeedOfSound.transpose(); // ft/sec
	postProcessMatrix.row(4) = 0.001941811*AirDensity.transpose(); // slug/ft3
	postProcessMatrix.row(5) = 1.94384*TranslationalVelocity.transpose(); // nmi_h
	postProcessMatrix.row(6) = Mach.transpose();

	postProcessMatrix.row(7) = AngleOfAttack.transpose(); // radians
	postProcessMatrix.row(8) = (180/3.14)*AngleOfAttack.transpose(); // degrees

	postProcessMatrix.row(9) = AngleOfSideslip.transpose(); // radians
	postProcessMatrix.row(10) = (180/3.14)*AngleOfSideslip.transpose(); // degrees

	postProcessMatrix.row(11) = u_NED.transpose(); // m/sec
	postProcessMatrix.row(12) = v_NED.transpose(); // m/sec
	postProcessMatrix.row(13) = w_NED.transpose(); // m/sec
	postProcessMatrix.row(14) = 3.28*u_NED.transpose(); // ft/sec
	postProcessMatrix.row(15) = 3.28*v_NED.transpose(); // ft/sec
	postProcessMatrix.row(16) = 3.28*w_NED.transpose(); // ft/sec

	postProcessMatrix.row(17) = phi_rad.transpose(); // rad
	postProcessMatrix.row(18) = theta_rad.transpose();
	postProcessMatrix.row(19) = psi_rad.transpose();
	postProcessMatrix.row(20) = phi_deg.transpose(); // deg
	postProcessMatrix.row(21) = theta_deg.transpose();
	postProcessMatrix.row(22) = psi_deg.transpose();
	
	postProcessMatrix.row(23) = roll_rate.transpose(); // rad/sec
	postProcessMatrix.row(24) = pitch_rate.transpose();
	postProcessMatrix.row(25) = yaw_rate.transpose();

	postProcessMatrix.row(26) = (180/3.14)*roll_rate.transpose(); // deg/sec
	postProcessMatrix.row(27) = (180/3.14)*pitch_rate.transpose();
	postProcessMatrix.row(28) = (180/3.14)*yaw_rate.transpose(); 


	//////////////////////////////////////////// end post processing ////////////////////////////////////////////




	 ////////////////////// Write solution matrix to CSV file //////////////////////
	std::ofstream outputFile("solution.csv");
	if (outputFile.is_open()) {
		outputFile << std::fixed << std::setprecision(10); // Set precision for output

		// Write header row (optional, but helpful)
		outputFile << "u_mps,v_mps,w_mps,p_rps,q_rps,r_rps,phi_rad,theta_rad,psi_rad,p1_n_m,p2_n_m,p3_n_m\n";

		// Write each row of the solution matrix to the CSV file
		for (int i = 0; i < resultMatrix.cols(); ++i) {
			for (int j = 0; j < resultMatrix.rows(); ++j) {
				outputFile << resultMatrix(j, i);
				if (j < resultMatrix.rows() - 1) {
					outputFile << ","; // Add comma as separator
				}
			}
			outputFile << "\n"; // Add newline after each row
		}
		outputFile.close();
		std::cout << "Solution matrix written to solution.csv" << std::endl;
	}
	else {
		std::cerr << "Unable to open file for writing: solution.csv" << std::endl;
		return 1;
	}



	std::ofstream outputPostProcessFile("postProcessSolution.csv");
	if (outputPostProcessFile.is_open()) {
		outputPostProcessFile << std::fixed << std::setprecision(10); // Set precision for output

		// Write header row (optional, but helpful)
		outputPostProcessFile << "Time, Altitude - m, Altitude - ft, SpeedofSound - ft/sec, AirDensity - slug/ft3,TransVel - ft/sec,Mach, AoA - rad, AoA - deg, AoS - rad, AoS - deg, u_NED, v_NED, w_NED, u_NED - ft/sec, v_NED - ft/sec, w_NED - ft/sec, phi_rad, theta_rad, psi_rad, phi_deg, theta_deg, psi_deg, roll_rate, pitch_rate, yaw_rate, roll_rate - deg, pitch_rate - deg, yaw_rate - deg \n";

		// Write each row of the solution matrix to the CSV file
		for (int i = 0; i < postProcessMatrix.cols(); ++i) {
			for (int j = 0; j < postProcessMatrix.rows(); ++j) {
				outputPostProcessFile << postProcessMatrix(j, i);
				if (j < postProcessMatrix.rows() - 1) {
					outputPostProcessFile << ","; // Add comma as separator
				}
			}
			outputPostProcessFile << "\n"; // Add newline after each row
		}
		outputPostProcessFile.close();
		std::cout << "Post process matrix written to postProcessSolution.csv" << std::endl;
	}
	else {
		std::cerr << "Unable to open file for writing: postProcessSolution.csv" << std::endl;
		return 1;
	}




}