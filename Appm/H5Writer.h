#pragma once

#include "H5Cpp.h"
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>
using namespace H5;


class H5Writer
{
public:
	H5Writer();
	H5Writer(const std::string & filename);
	~H5Writer();
	
	void writeData(const Eigen::MatrixXd & matrix, const std::string & dataname);

private:
	std::string filename;
	H5File file;
};

