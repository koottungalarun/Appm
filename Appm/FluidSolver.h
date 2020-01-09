#pragma once

#include <Eigen/Dense>

#include "DualMesh.h"
#include "FluidState.h"
#include "Numerics.h"
#include "H5Writer.h"

class FluidSolver
{
public:
	FluidSolver();
	FluidSolver(const DualMesh * mesh);
	~FluidSolver();

	virtual const double updateFluidState() = 0;

	void writeStates(H5Writer & writer) const;


protected:
	const DualMesh * mesh = nullptr;

	Eigen::MatrixXd fluidStates;
	Eigen::MatrixXd fluidFluxes;


private:
	void init();
};

