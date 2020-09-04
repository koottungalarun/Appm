#pragma once

#include "interpolation.h"
#include "DataTransform.h"
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <cassert>

/**
* Interpolation of scalar-valued function in two dimension, i.e., 
* v = f(x,y) where x, y, and v are scalars
*
* This class bluids on the open-source implementation of ALGLIB project.
*
* The data is provided ona 
* 
*/

class Interpolation2d
{

public:

	//enum class DataTransform {
	//	NONE, INVERSE, LOG
	//};


	Interpolation2d();
	

	/**
	* Read 2D data from csv file.
	*
	* The data is assumed to be given on a cartesian grid. 
	* - The first row has the x-coordinates of the grid
	* - The first column has the y-coordiantes of the grid
	* - The data is given at f_ij = f(y(i), x(j)), where f_ij is the data point at i-th row and j-th column
	* - The top left corner is a dummy value
	*
	* @param filename    name of data file
	* @param xTrans      transformation to data on x-axis
	* @param yTrans      transformation to data on y-axis
	* @param fTrans      transformation to function data
	*/
	Interpolation2d(const std::string & filename, DataTransform xTrans, DataTransform yTrans, DataTransform fTrans);


	~Interpolation2d();

	/** 
	* Interpolate data using bicubic splines z = S(x,y). 
	*
	* @param xSites vector of x(i)
	* @param ySites vector of y(i)
	* @return vector z(i) = S(x(i), y(i))
	*/
	Eigen::VectorXd bicubicInterp(const Eigen::VectorXd & xSites, const Eigen::VectorXd & ySites);
	


private:

	// Cartesian grid axis
	alglib::real_1d_array x;
	// length of data
	alglib::ae_int_t n = 0;
	// transformation to data applied when reading from file
	DataTransform xTrans;

	// Cartesian grid axis
	alglib::real_1d_array y;
	// length of data
	alglib::ae_int_t m = 0;
	// transformation to data applied when reading from file
	DataTransform yTrans;


	// Function values
	alglib::real_1d_array f;
	// transformation to data applied when reading from file
	DataTransform fTrans;

	// Function dimension (d = 1 means a scalar-valued function)
	alglib::ae_int_t d = 1;

	// Data structure that stores information of spline interpolation
	alglib::spline2dinterpolant s;


	/**
	* Read data from csv file into internal data structure.
	*/
	void readCsvFile(const std::string & filename, const DataTransform & xTrans, const DataTransform & yTrans, const DataTransform & fTrans);

	void writeCsvFile(const std::string & filename);

	template<class T>
	T applyTransform(const T array, const DataTransform & trans);
};

