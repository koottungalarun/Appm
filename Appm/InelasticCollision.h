#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include "Interpolation1d.h"
#include "Interpolation2d.h"
#include "DataTransform.h"
#include "Physics.h"
#include "PhysicsConstants.h"
#include "ScalingParameters.h"
#include <fstream>
#include <Eigen/Dense>

class InelasticCollision
{
public:
	InelasticCollision();

	/* Read data for inelastic collisions from folder */
	InelasticCollision(const std::string & folderPath, 
		const int idxE_, const int idxA_, const int idxI_, 
		const ScalingParameters & scales);

	~InelasticCollision();

	const Eigen::VectorXd getGion(const Eigen::VectorXd & nE, const Eigen::VectorXd & nA, const Eigen::VectorXd & vthE, const Eigen::VectorXd & TeVec, const Eigen::VectorXd & lambdaIon) const;
	const Eigen::VectorXd getGrec(const Eigen::VectorXd & nI, const Eigen::VectorXd & nE, const Eigen::VectorXd & vthE, const Eigen::VectorXd & xStar, const double mE, const Eigen::VectorXd & Te, const Eigen::VectorXd & lambdaVec) const;
	const Eigen::VectorXd getR0ion(const Eigen::VectorXd & nE, const Eigen::VectorXd & nA, const Eigen::VectorXd & vthE, const Eigen::VectorXd & TeVec, const Eigen::VectorXd & lambdaIon) const;
	const Eigen::VectorXd getJ00ion(const Eigen::VectorXd & nE, const Eigen::VectorXd & nA, const Eigen::VectorXd & vthE, const Eigen::VectorXd & TeVec, const Eigen::VectorXd & lambdaIon) const;
	const Eigen::VectorXd getR1rec(const Eigen::VectorXd & nE, const Eigen::VectorXd & nI, const Eigen::VectorXd & vthE, const Eigen::VectorXd & xStar, const Eigen::VectorXd & TeVec, const Eigen::VectorXd & lambdaRec) const;
	const Eigen::VectorXd getR2rec(const Eigen::VectorXd & nE, const Eigen::VectorXd & nI, const Eigen::VectorXd & vthE, const Eigen::VectorXd & xStar, const Eigen::VectorXd & TeVec, const Eigen::VectorXd & lambdaRec) const;
	const Eigen::VectorXd getJ11rec(const Eigen::VectorXd & nE, const Eigen::VectorXd & nI, const Eigen::VectorXd & vthE, const Eigen::VectorXd & xStar, const Eigen::VectorXd & TeVec, const Eigen::VectorXd & lambdaRec) const;
	const Eigen::VectorXd getJ22rec(const Eigen::VectorXd & nE, const Eigen::VectorXd & nI, const Eigen::VectorXd & vthE, const Eigen::VectorXd & xStar, const Eigen::VectorXd & TeVec, const Eigen::VectorXd & lambdaRec) const;
	const Eigen::VectorXd getJ12rec(const Eigen::VectorXd & nE, const Eigen::VectorXd & nI, const Eigen::VectorXd & vthE, const Eigen::VectorXd & xStar, const Eigen::VectorXd & TeVec, const Eigen::VectorXd & lambdaRec) const;

	const int getElectronFluidx() const;
	const int getIonFluidx() const;
	const int getAtomFluidx() const;

	const double getIonizationEnergyScaled() const;

	const Eigen::VectorXd getGionInterpolated(const Eigen::VectorXd & Te) const;
	const Eigen::VectorXd getR0ionInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda) const;
	const Eigen::VectorXd getJ00ionInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda) const;
	const Eigen::VectorXd getGrecInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda) const;
	const Eigen::VectorXd getR1recInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda) const;
	const Eigen::VectorXd getR2recInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda) const;
	const Eigen::VectorXd getJ11recInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda) const;
	const Eigen::VectorXd getJ22recInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda) const;
	const Eigen::VectorXd getJ12recInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda) const;

private:
	const int idxE = -1;
	const int idxA = -1;
	const int idxI = -1;

	double ionizationEnergyScaled = 0;

	int g0 = 0;
	int g1 = 0;
	double nScale = 0;
	double thermalDeBroglieScale = 0;

	/* Scale for cross sections */
	double sigmaScale = 0;

	Interpolation1d data_Gion;
	Interpolation2d data_R0ion;
	Interpolation2d data_J00ion;

	Interpolation2d data_Grec;
	Interpolation2d data_R1rec;
	Interpolation2d data_R2rec;
	Interpolation2d data_J11rec;
	Interpolation2d data_J22rec;
	Interpolation2d data_J12rec;

};


