#pragma once

#include "interpolation.h"
#include "DataTransform.h"
#include <Eigen/Dense>

class Interpolation1d
{
public:
	Interpolation1d();
	Interpolation1d(const std::string & filename, const DataTransform xTrans, const DataTransform fTrans);

	~Interpolation1d();

	Eigen::VectorXd bicubicInterp(const Eigen::VectorXd & xSites);



private:
	void readCsvFile(const std::string & filename, const DataTransform & xTrans, const DataTransform & fTrans);


	template<class T>
	T applyTransform(const T array, const DataTransform & trans);
};

