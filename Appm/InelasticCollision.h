#pragma once

#include "InterpolationTable.h"
#include <fstream>
#include <Eigen/Dense>

class InelasticCollision
{
public:
	InelasticCollision();
	InelasticCollision(const std::string & filename);

	~InelasticCollision();

	Eigen::VectorXd getIonizationRate(const Eigen::VectorXd & T);

	Eigen::MatrixXd getData() const;

	int getElectronFluidx() const;
	int getAtomFluidx() const;
	int getIonFluidx() const;

private:
	InterpolationTable * table = nullptr;
	int fluidxElectrons = -1;
	int fluidxAtoms = -1;
	int fluidxIons = -1;

};

