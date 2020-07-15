#include "InterpolationTable.h"

InterpolationTable::InterpolationTable()
{
}

InterpolationTable::InterpolationTable(const std::vector<double> x, const std::vector<double> y)
{
	assert(x.size() > 0);
	assert(y.size() > 0);
	assert(x.size() == y.size());
	const int N = x.size();
	this->x = Eigen::VectorXd(N);
	this->y = Eigen::VectorXd(N);
	for (int i = 0; i < N; i++) {
		this->x(i) = x[i];
		this->y(i) = y[i];
	}
	init();
}

InterpolationTable::InterpolationTable(Eigen::VectorXd x, Eigen::VectorXd y)
{
	assert(x.size() > 0);
	assert(y.size() > 0);
	assert(x.size() == y.size());
	this->x = x;
	this->y = y;
	init();
}

InterpolationTable::~InterpolationTable()
{
	if (task != nullptr) {
		int status = -1;
		status = dfDeleteTask(&task);
		assert(errorCheck(status));
		task = nullptr;
	}
}

const Eigen::VectorXd InterpolationTable::interpolate(const Eigen::VectorXd & sites)
{
	const int N = sites.size();
	Eigen::VectorXd result;
	result = Eigen::VectorXd(N);
	result.setZero();

	const double sitesMax = sites.array().maxCoeff();
	const double sitesMin = sites.array().minCoeff();
	const double xMax = x.array().maxCoeff();
	const double xMin = x.array().minCoeff();
	if (sitesMax > xMax || sitesMin < xMin) {
		std::cout << "Warning: interpolation table range exceeded." << std::endl;
		std::cout << "sites min,max: " << sitesMin << ",\t" << sitesMax << std::endl;
		std::cout << "xdata min,max: " << xMin << ",\t" << xMax << std::endl;
	}

	int status;
	status = dfdInterpolate1D(task, DF_INTERP, DF_METHOD_PP, N, sites.data(), sitehint, ndorder, &dorder, datahint, result.data(), rhint, cell);
	assert(errorCheck(status));
	return result;
}

const Eigen::VectorXd InterpolationTable::getXdata() const
{
	return this->x;
}

const Eigen::VectorXd InterpolationTable::getYdata() const
{
	return this->y;
}

void InterpolationTable::init()
{
	int status;
	MKL_INT nx = x.size(); // size of partition
	MKL_INT xhint = DF_NON_UNIFORM_PARTITION; // Additional info on x-data 
	MKL_INT ny = 1; // function dimension (ny = 1 is a scalar-valued function)
	MKL_INT yhint = DF_NO_HINT; // additional info on y-values 

	status = dfdNewTask1D(&task, nx, x.data(), xhint, ny, y.data(), yhint);
	assert(errorCheck(status));

	scoeff = std::vector<double>((nx - 1) * s_order); // spline coefficients
	assert(scoeff.size() > 0);
	status = dfdEditPPSpline1D(task, s_order, s_type, bc_type, bc, ic_type, ic, scoeff.data(), scoeffhint);
	errorCheck(status);

	status = dfdConstruct1D(task, DF_PP_SPLINE, DF_METHOD_STD);
	assert(errorCheck(status));
}

bool InterpolationTable::errorCheck(const int status)
{
	bool isValid = status == DF_STATUS_OK;
	if (!isValid) {
		std::cerr << "status = " << status << std::endl;
	}
	return isValid;
}
