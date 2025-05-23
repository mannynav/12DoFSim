#pragma once

std::vector<double> generate_arange(double start, double stop, double step) {
	std::vector<double> result;
	for (double i = start; i < stop; i += step) {
		result.push_back(i);
	}
	return result;
}

double linearInterpolation(const std::vector<double>& x_values, const std::vector<double>& y_values, double x) {
	// Error handling: Check if input vectors are valid
	if (x_values.empty() || y_values.empty()) {
		throw std::invalid_argument("Input vectors cannot be empty.");
	}
	if (x_values.size() != y_values.size()) {
		throw std::invalid_argument("Input vectors must have the same size.");
	}

	// Find the interval containing x
	auto it = std::upper_bound(x_values.begin(), x_values.end(), x);

	if (it == x_values.begin()) {
		// x is smaller than the first x_value
		return y_values.front();
	}
	if (it == x_values.end()) {
		// x is larger than the last x_value
		return y_values.back();
	}

	// Get the indices of the two points surrounding x
	size_t index_upper = std::distance(x_values.begin(), it);
	size_t index_lower = index_upper - 1;

	// Get the x and y values of the surrounding points
	double x1 = x_values[index_lower];
	double y1 = y_values[index_lower];
	double x2 = x_values[index_upper];
	double y2 = y_values[index_upper];

	// Perform linear interpolation
	if (x2 == x1) {
		// Avoid division by zero if x1 and x2 are the same
		return y1;
	}
	return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
}