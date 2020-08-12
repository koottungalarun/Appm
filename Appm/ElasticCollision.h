#pragma once

#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include "InterpolationTable.h"

class ElasticCollision
{

public:
	ElasticCollision();
	ElasticCollision(const std::string & filename, const int idxA, const int idxB);
	ElasticCollision(const std::string & filename, const int idxA, const int idxB, const double Tscale, const double crossSecScale);
	~ElasticCollision();

	const int getFluidIdxA() const;
	void setFluidIdxA(const int idxA);
	const int getFluidIdxB() const;
	void setFluidIdxB(const int idxB);

	const Eigen::VectorXd getAvgMomCrossSection(const Eigen::VectorXd & T) const;

	const Eigen::MatrixXd getData() const;

private:
	int fluidxA = -1;
	int fluidxB = -1;
	InterpolationTable * table = nullptr;

};

