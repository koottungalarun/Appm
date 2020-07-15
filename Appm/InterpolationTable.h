#pragma once

#include <vector>
#include <mkl.h>
#include <Eigen/Dense>
#include <iostream>

class InterpolationTable
{
public:
	InterpolationTable();
	InterpolationTable(const std::vector<double> x, const std::vector<double> y);
	InterpolationTable(Eigen::VectorXd x, Eigen::VectorXd y);
	~InterpolationTable();

	const Eigen::VectorXd interpolate(const Eigen::VectorXd & sites);

	const Eigen::VectorXd getXdata() const;
	const Eigen::VectorXd getYdata() const;


private:
	DFTaskPtr task = nullptr;
	MKL_INT s_order = DF_PP_CUBIC;   // spline order
	MKL_INT s_type = DF_PP_NATURAL;  // spline type 
	MKL_INT ic_type = DF_NO_IC; // internal conditions
	MKL_INT bc_type = DF_BC_NOT_A_KNOT; // boundry conditions

	double * bc = nullptr;
	double * ic = nullptr;

	std::vector<double> scoeff;
	MKL_INT scoeffhint = DF_NO_HINT;

	MKL_INT sitehint = DF_NON_UNIFORM_PARTITION;

	MKL_INT ndorder = 1; // maximal derivative order computed, increased by 1
	MKL_INT dorder = 1;  // order of derivatives to be computed

	double * datahint = DF_NO_APRIORI_INFO; // additional info on partition and interpolation sites
	MKL_INT rhint = DF_MATRIX_STORAGE_ROWS;  // additional info on results
	MKL_INT * cell = nullptr; // array of cell indices


	Eigen::VectorXd x;
	Eigen::VectorXd y;

	void init();
	bool errorCheck(const int status);
};

