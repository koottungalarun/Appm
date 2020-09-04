#include "Interpolation1d.h"

Interpolation1d::Interpolation1d()
{
}

Interpolation1d::Interpolation1d(const std::string & filename, const DataTransform xTrans, const DataTransform fTrans)
{
	readCsvFile(filename, xTrans, fTrans);
}

Interpolation1d::~Interpolation1d()
{
}

Eigen::VectorXd Interpolation1d::bicubicInterp(const Eigen::VectorXd & xSites)
{
	return Eigen::VectorXd(xSites.size());
}

void Interpolation1d::readCsvFile(const std::string & filename, const DataTransform & xTrans, const DataTransform & fTrans)
{
}


template <class T>
T Interpolation1d::applyTransform(const T array, const DataTransform & trans)
{
	T transformed = array;
	if (trans == DataTransform::INVERSE) {
		transformed = array.array().inverse();
	}
	if (trans == DataTransform::LOG) {
		transformed = array.array().log();
	}
	return transformed;
}

