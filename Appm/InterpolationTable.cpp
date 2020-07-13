#include "InterpolationTable.h"

InterpolationTable::InterpolationTable()
{
}

InterpolationTable::InterpolationTable(const std::vector<double> x, const std::vector<double> y)
{
	assert(x.size() > 0);
	assert(y.size() > 0);
	this->x = x;
	this->y = y;
}

InterpolationTable::~InterpolationTable()
{
}

const Eigen::VectorXd InterpolationTable::interpolate(const Eigen::VectorXd & sites)
{
	const int N = sites.size();
	Eigen::VectorXd result;
	result = Eigen::VectorXd(N);
	result.setZero();
	return result;
}
