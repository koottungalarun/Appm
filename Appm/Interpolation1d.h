#pragma once

#include "interpolation.h"
#include "DataTransform.h"
#include <Eigen/Dense>
#include <fstream>

class Interpolation1d
{
public:
	Interpolation1d();
	Interpolation1d(const std::string & filename, const DataTransform xTrans, const DataTransform fTrans);

	~Interpolation1d();

	const Eigen::VectorXd cubicInterp(const Eigen::VectorXd & xSites) const;



private:
	alglib::real_1d_array x;
	alglib::ae_int_t n = 0;
	DataTransform xTrans;

	alglib::real_1d_array f;
	DataTransform fTrans;

	alglib::spline1dinterpolant s;

	void readCsvFile(const std::string & filename, const DataTransform & xTrans, const DataTransform & fTrans);


	template<class T>
	T applyTransform(const T array, const DataTransform & trans);
};

