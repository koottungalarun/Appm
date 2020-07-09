#pragma once

#include <vector>
#include <mkl.h>
#include <Eigen/Dense>

class InterpolationTable
{
	InterpolationTable();
	InterpolationTable(const std::vector<double> x, const std::vector<double> y);
	~InterpolationTable();

	const Eigen::VectorXd interpolate(const Eigen::VectorXd & sites);

private:
	std::vector<double> x;
	std::vector<double> y;
};

