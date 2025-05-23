#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

void printVector(std::vector<double>& vec, std::string& name)
{

	std::cout << "------------ " << name << " -------------------" << std::endl;
	for (int i = 0; i < vec.size(); i++) {
		std::cout << "Time i: " << i << " Value: " << vec[i] << std::endl;
	}

}