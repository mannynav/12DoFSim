#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

std::vector<double> readNumericCSVColumn(const std::string& filename, int columnIndex = 0) {
	std::vector<double> data;
	std::ifstream file(filename);

	if (!file.is_open()) {
		throw std::runtime_error("Could not open file: " + filename);
	}

	std::string line;
	//Skip the header
	if (std::getline(file, line)) {}

	while (std::getline(file, line)) {
		std::stringstream lineStream(line);
		std::string cell;
		std::vector<std::string> rowData;

		while (std::getline(lineStream, cell, ',')) {
			rowData.push_back(cell);
		}

		if (columnIndex >= rowData.size()) {
			throw std::runtime_error("Column index out of range in row: " + line);
		}

		try {
			double value = std::stod(rowData[columnIndex]);
			data.push_back(value);
		}
		catch (const std::invalid_argument& e) {
			throw std::runtime_error("Invalid numeric value in row: " + line + ", column: " + std::to_string(columnIndex));
		}
		catch (const std::out_of_range& e) {
			throw std::runtime_error("Numeric value out of range in row: " + line + ", column: " + std::to_string(columnIndex));
		}
	}

	return data;
}