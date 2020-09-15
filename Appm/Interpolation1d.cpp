#include "Interpolation1d.h"

Interpolation1d::Interpolation1d()
{
}

Interpolation1d::Interpolation1d(const std::string & filename, const DataTransform xTrans_, const DataTransform fTrans_, const double xScale, const double fScale)
	: xTrans(xTrans_), fTrans(fTrans_)
{
	readCsvFile(filename, xTrans, fTrans, xScale, fScale);

	assert(n > 0);
	assert(x.length() > 0);
	assert(f.length() > 0);
	assert(x.length() == f.length());

	alglib::spline1dbuildcubic(x, f, s);
}

Interpolation1d::~Interpolation1d()
{
}

const Eigen::VectorXd Interpolation1d::cubicInterp(const Eigen::VectorXd & xSites) const
{
	const int n = xSites.size();
	Eigen::VectorXd result = Eigen::VectorXd(n);
	Eigen::VectorXd xSitesTrans = xSites;
	Eigen::VectorXd resultTrans = result;

	// Apply transformation to xSites
	switch (xTrans) {
	case DataTransform::NONE:
		// do nothing
		break;
		
	case DataTransform::INVERSE:
		xSitesTrans = xSites.array().inverse();
		break;
		
	case DataTransform::LOG:
		xSitesTrans = xSites.array().exp();
		break;

	default:
		assert(false);
	}

	// Interpolate
	for (int i = 0; i < n; i++) {
		result(i) = alglib::spline1dcalc(s, xSitesTrans(i));
	}

	// Apply transformation on function data
	switch (fTrans) {
	case DataTransform::NONE:
		// do nothing
		break;

	case DataTransform::LOG:
		resultTrans = result.array().exp();
		break;

	case DataTransform::INVERSE:
		resultTrans = result.array().exp();
		break;
		
	default:
		assert(false);
	}
	return resultTrans;
}

void Interpolation1d::readCsvFile(const std::string & filename, const DataTransform & xTrans, const DataTransform & fTrans, const double xScale, const double fScale)
{
	assert(xTrans == DataTransform::NONE);

	std::ifstream file(filename);
	assert(file.is_open());

	// Read each line of data file
	std::string line;
	int rows = 0;
	int cols = 0;
	std::vector<double> data;

	while (std::getline(file, line))
	{
		if (line.empty() || line.front() == '#') {
			continue; // skip empty lines and comment lines
		}

		// Read a single line 
		std::string chunk;
		std::stringstream ss(line);
		while (std::getline(ss, chunk, ',')) {
			double value = std::stod(chunk);
			data.push_back(value);
		}

		// get number of columns
		if (rows == 0) {
			cols = data.size();
		}
		rows++;
	}
	assert(data.size() == (rows * cols));
	assert(cols == 2);

	// Note that we read data in row-major format. 
	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXdRowMajor;
	Eigen::Map<MatrixXdRowMajor> map(data.data(), rows, cols);

	Eigen::VectorXd firstRow = map.row(0);
	Eigen::VectorXd firstCol = map.col(0);

	Eigen::VectorXd TeVec = firstCol.segment(1, firstCol.size() - 1);
	TeVec *= xScale;

	MatrixXdRowMajor dataMatrix = map.bottomRightCorner(rows - 1, cols - 1);
	dataMatrix *= fScale;
	if (fTrans != DataTransform::NONE) {
		dataMatrix = applyTransform(dataMatrix, fTrans);
	}
	if (xTrans != DataTransform::NONE) {
		TeVec = applyTransform(TeVec, xTrans);
	}

	this->f.setcontent(dataMatrix.size(), dataMatrix.data());
	this->x.setcontent(TeVec.size(), TeVec.data());

	//Eigen::VectorXd TeVec_inv = TeVec.array().inverse();
	//MatrixXdRowMajor dataMatrix_log = dataMatrix.array().log();
	//this->f.setcontent(dataMatrix_log.size(), dataMatrix_log.data());
	//this->y.setcontent(TeVec_inv.size(), TeVec_inv.data());


	this->n = this->x.length();
	assert(this->n > 0);
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

