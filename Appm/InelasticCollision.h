#pragma once

#include "InterpolationTable.h"
#include <fstream>
#include <Eigen/Dense>
#include "PhysicsConstants.h"
#include "ScalingParameters.h"

class InelasticCollision
{
public:
	InelasticCollision();
	InelasticCollision(const std::string & filename);
	InelasticCollision(const std::string & filename, const double kScale, const double Tscale);

	~InelasticCollision();

	Eigen::VectorXd getIonizationRate(const Eigen::VectorXd & T);
	Eigen::VectorXd getRecombinationRate_Saha(const Eigen::VectorXd & ki, const Eigen::VectorXd & Te);

	Eigen::MatrixXd getData() const;

	int getElectronFluidx() const;
	void setElectronFluidx(const int idx);
	int getAtomFluidx() const;
	void setAtomFluidx(const int idx);
	int getIonFluidx() const;
	void setIonFluidx(const int idx);
	void setNumberDensityScale(const double nbar);
	void setElectronMassRatio(const double electronMassRatio);
	void setMassScale(const double mbar);
	void setTemperatureScale(const double Tbar);
	void setScalingParameters(const ScalingParameters & params);

private:
	InterpolationTable * table = nullptr;
	int fluidxElectrons = -1;
	int fluidxAtoms = -1;
	int fluidxIons = -1;
	//double nbar = 0;
	//double electronMassRatio = 0;
	//double mbar = 0;
	//double Tbar = 0;
	//double deBroglieScaled_pow3 = 0;
	double recomb_Saha_coeff = 0;

	const double getRecombSahaCoeff() const;
	void setRecombSahaCoeff(const double a);

	//const double getNumberDensityScale() const;
	//const double getElectronMassRatio() const;
	//const double getDeBroglieScaled_pow3() const;
};

