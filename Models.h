#pragma once

#include <Eigen/Dense>
#include <map>
#include <cmath>
#include <numbers>
#include <math.h>

std::map<std::string, double> NASA_Atmos01_Sphere() {

	double pi = 3.14159;

	double in2m = 0.0254;

	double rho_NASA_Sphere = 7868.36;

	double r_sphere_in = 3.0;
	double r_sphere_m = r_sphere_in * in2m;

	double b_m = r_sphere_m; //ref wing span

	double c_m = r_sphere_m; //ref wing chord

	double CD_approx = 0.5; //Approximate coef of drag for subsonic regimes

	double Clp = 0.0;

	double Clr = 0.0;

	double Cmq = 0.0;

	double Cnp = 0.0;

	double Cnr = 0.0;

	long double vol_sphere = (4.0 / 3.0) * pi * pow(r_sphere_m, 3);

	double mass = rho_NASA_Sphere * vol_sphere;

	double J_sphere = 0.4 * mass * pow(r_sphere_m, 2);

	double A_ref = pi * pow(rho_NASA_Sphere, 2);

	double Vterm = sqrt((2 * mass * 9.81) / (1.2 * CD_approx * A_ref));

	std::map<std::string, double> _vehicle_map{};

	_vehicle_map["m_kg"] = mass;
	_vehicle_map["Jxz_b"] = 0;
	_vehicle_map["Jxx_b"] = J_sphere;
	_vehicle_map["Jyy_b"] = J_sphere;
	_vehicle_map["Jzz_b"] = J_sphere;
	_vehicle_map["r"] = r_sphere_m;
	_vehicle_map["mass"] = mass;
	_vehicle_map["CD_approx"] = CD_approx;
	_vehicle_map["Clp"] = Clp;
	_vehicle_map["Clr"] = Clr;
	_vehicle_map["Cmq"] = Cmq;
	_vehicle_map["Cnp"] = Cnp;
	_vehicle_map["Cnr"] = Cnr;
	_vehicle_map["Aref"] = A_ref;
	_vehicle_map["Vterm"] = Vterm;
	_vehicle_map["b_m"] = b_m;
	_vehicle_map["c_m"] = c_m;

	return _vehicle_map;
}
std::map<std::string, double> NASA_Atmos02_Brick() {

	double pi = 3.14159;

	double in2m = 0.0254;

	double slug2kg = 14.5939;
	double kg2slug = 1 / slug2kg;

	double mass_slug = 0.1554048;
	double mass_kg = kg2slug * mass_slug;


	double ft2m = 0.304878;
	double Jxx_slug_ft2 = 0.00189422;
	double Jxx_kgm2 = slug2kg * (ft2m * ft2m) * Jxx_slug_ft2;

	double Jyy_slug_ft2 = 0.00621102;
	double Jyy_kgm2 = slug2kg * (ft2m * ft2m) * Jyy_slug_ft2;

	double Jzz_slug_ft2 = 0.00719467;
	double Jzz_kgm2 = slug2kg * (ft2m * ft2m) * Jzz_slug_ft2;

	double Jzx_slug_ft2 = 0.0;
	double Jzx_kgm2 = slug2kg * (ft2m * ft2m) * Jzx_slug_ft2;


	double body_xpos_CM = 0.0;
	double body_ypos_CM = 0.0;
	double body_zpos_CM = 0.0;

	double length = 8*in2m; //m
	double width = 4 * in2m; //m
	double A_ref = length * width;

	double ref_wing_span = 0.33333 * ft2m;
	double ref_wing_chord = 0.66667 * ft2m;

	double Clp = 0.0; //roll rate roll damping
	double Clr = 0.0; //yaw rate roll damping
	double Cmq = 0.0; // pitch rate pitch damping
	double Cnp = 0.0; //roll rate yaw damping
	double Cnr = 0.0; //yaw rate yaw damping



	std::map<std::string, double> _vehicle_map{};

	_vehicle_map["m_kg"] = mass_kg;
	_vehicle_map["Jxz_b"] = Jzx_kgm2;
	_vehicle_map["Jxx_b"] = Jxx_kgm2;
	_vehicle_map["Jyy_b"] = Jyy_kgm2;
	_vehicle_map["Jzz_b"] = Jzz_kgm2;
	_vehicle_map["Area_ref"] = A_ref;
	_vehicle_map["Clp"] = Clp;
	_vehicle_map["Clr"] = Clr;
	_vehicle_map["Cmq"] = Cmq;
	_vehicle_map["Cnp"] = Cnp;
	_vehicle_map["Cnr"] = Cnr;
	_vehicle_map["b_m"] = ref_wing_span;
	_vehicle_map["c_m"] = ref_wing_chord;

	return _vehicle_map;
}
std::map<std::string, double> NASA_Atmos03_Brick() {

	double pi = 3.14159;

	double in2m = 0.0254;

	double slug2kg = 14.5939;
	double kg2slug = 1 / slug2kg;

	double mass_slug = 0.1554048;
	double mass_kg = kg2slug * mass_slug;


	double ft2m = 0.304878;
	double Jxx_slug_ft2 = 0.00189422;
	double Jxx_kgm2 = slug2kg * (ft2m * ft2m) * Jxx_slug_ft2;

	double Jyy_slug_ft2 = 0.00621102;
	double Jyy_kgm2 = slug2kg * (ft2m * ft2m) * Jyy_slug_ft2;

	double Jzz_slug_ft2 = 0.00719467;
	double Jzz_kgm2 = slug2kg * (ft2m * ft2m) * Jzz_slug_ft2;

	double Jzx_slug_ft2 = 0.0;
	double Jzx_kgm2 = slug2kg * (ft2m * ft2m) * Jzx_slug_ft2;


	double body_xpos_CM = 0.0;
	double body_ypos_CM = 0.0;
	double body_zpos_CM = 0.0;

	double length = 8 * in2m; //m
	double width = 4 * in2m; //m
	double A_ref = length * width;

	double ref_wing_span = 0.33333 * ft2m;
	double ref_wing_chord = 0.66667 * ft2m;

	double Clp = -1.0; //roll rate roll damping
	double Clr = 0.0; //yaw rate roll damping
	double Cmq = -1.0; // pitch rate pitch damping
	double Cnp = 0.0; //roll rate yaw damping
	double Cnr = -1.0; //yaw rate yaw damping



	std::map<std::string, double> _vehicle_map{};

	_vehicle_map["m_kg"] = mass_kg;
	_vehicle_map["Jxz_b"] = Jzx_kgm2;
	_vehicle_map["Jxx_b"] = Jxx_kgm2;
	_vehicle_map["Jyy_b"] = Jyy_kgm2;
	_vehicle_map["Jzz_b"] = Jzz_kgm2;
	_vehicle_map["Area_ref"] = A_ref;
	_vehicle_map["Clp"] = Clp;
	_vehicle_map["Clr"] = Clr;
	_vehicle_map["Cmq"] = Cmq;
	_vehicle_map["Cnp"] = Cnp;
	_vehicle_map["Cnr"] = Cnr;
	_vehicle_map["b_m"] = ref_wing_span;
	_vehicle_map["c_m"] = ref_wing_chord;

	return _vehicle_map;
}


//To be implemented
std::map<std::string, double> F16() {


	std::map<std::string, double> _vehicle_map{};

	/*_vehicle_map["m_kg"] = mass_kg;
	_vehicle_map["Jxz_b"] = Jzx_kgm2;
	_vehicle_map["Jxx_b"] = Jxx_kgm2;
	_vehicle_map["Jyy_b"] = Jyy_kgm2;
	_vehicle_map["Jzz_b"] = Jzz_kgm2;
	_vehicle_map["Area_ref"] = A_ref;
	_vehicle_map["Clp"] = Clp;
	_vehicle_map["Clr"] = Clr;
	_vehicle_map["Cmq"] = Cmq;
	_vehicle_map["Cnp"] = Cnp;
	_vehicle_map["Cnr"] = Cnr;
	_vehicle_map["b_m"] = ref_wing_span;
	_vehicle_map["c_m"] = ref_wing_chord;*/

	return _vehicle_map;
}


class X15
{

public:

	//Aerodynamic coefficients
	//The drag, lift and sideforce coefficent build ups are from "Adaptive Control and the NASA X-15-3 Flight revisited".
	//When these functions are called, they will be dimensionalized
	double drag_coefficient(double CD_wing_body, double CD_delta_elevator, double CD_delta_speed_brake, double elevator_angle, double speed_brake_angle) {

		double CD = CD_wing_body + CD_delta_elevator * elevator_angle + CD_delta_speed_brake * speed_brake_angle;
		return CD;
	}


	//From Walker60, CL_dele is a function of AoA, dele and Mach
	double lift_coefficient(double CL_wing_body, double CL_alpha, double CL_dele, double alpha, double elevator_angle) {

		double CL = CL_wing_body + CL_alpha * alpha + CL_dele * (alpha + elevator_angle);
		return CL;
	}
	double sideforce_coefficient(double CY_beta_sideslip, double CY_p_roll, double CY_r_yaw, double CY_beta_dot_accel, double CY_delta_aileron, double CY_delta_rudder, double beta, double p, double r, double beta_dot, double aileron_angle, double rudder_angle, double true_velocity, double b_m)
	{
		double CY = CY_beta_sideslip * beta + CY_p_roll * p * b_m / (2 * true_velocity) + CY_r_yaw * r * b_m / (2 * true_velocity) +
			CY_beta_dot_accel * beta_dot * b_m / (2 * true_velocity) + CY_delta_aileron * aileron_angle + CY_delta_rudder * rudder_angle;
		return CY;
	}




	//Moment coefficients - roll, pitch and yaw. Also from "Adaptive Control and the NASA X-15-3 Flight revisited"
	double roll_coefficient() {

	}

	double pitch_coefficient() {

	}

	double yaw_coefficient() {

	}



private:

};


class Vehicle
{
private:

	const double in2m = 0.0254; //inches to meter conversion
	std::string _name{};
	double _rho{}; //density
	double _radius_in{};
	double _radius_m{};
	double _coef_drag{};
	std::map<std::string, double> _vehicle_map{};

	long double pi = 3.14159;

	long double volume(long double radius) {

		return (4.0 / 3.0) * pi * pow(radius, 3);

	}

	double mass(double radius, double rho) {
		return rho * volume(radius);
	}

	double J(double radius, double rho) {
		return 0.4 * mass(radius, rho) * mass(radius, rho);
	}

	double Aref(double radius, double rho) {
		return pi * radius * radius;
	}


public:

	Vehicle(double rho, double radius, double coef_drag, std::string name) : _name(name), _rho(rho), _radius_in(radius), _coef_drag(coef_drag) {

		_radius_m = _radius_in * in2m;

		double Vterm = sqrt((2 * mass(_radius_m, _rho) * 9.81) / (1.2 * _coef_drag * Aref(_radius_m, _rho)));

		//_vehicle_map["name"] = _name;
		_vehicle_map["m_kg"] = mass(_radius_m, _rho);
		_vehicle_map["Jxz_b"] = 0;
		_vehicle_map["Jxx_b"] = J(_radius_m, _rho);
		_vehicle_map["Jyy_b"] = J(_radius_m, _rho);
		_vehicle_map["Jzz_b"] = J(_radius_m, _rho);
		_vehicle_map["r"] = _radius_m;
		_vehicle_map["mass"] = mass(_radius_m, rho);
		_vehicle_map["CD_approx"] = _coef_drag;
		_vehicle_map["Aref"] = Aref(_radius_m, rho);
		_vehicle_map["Vterm"] = Vterm;

		_vehicle_map["Clp"] = 0;
		_vehicle_map["Clr"] = 0;
		_vehicle_map["Cmq"] = 0;
		_vehicle_map["Cnp"] = 0;
		_vehicle_map["Cnr"] = 0;
		_vehicle_map["b_m"] = 0;
		_vehicle_map["c_m"] = 0;



		std::cout << "radius_m: " << _radius_m << std::endl;
		std::cout << "Volume: " << volume(_radius_m) << std::endl;
		std::cout << "mass: " << mass(_radius_m, _rho) << std::endl;
		std::cout << "Aref: " << Aref(_radius_m, rho) << std::endl;

	}

	std::map<std::string, double> getVehicleMap() {
		return _vehicle_map;
	}

};