#pragma once
#include "FluidSolver.h"
class SingleFluidSolver :
	public FluidSolver
{
public:
	SingleFluidSolver();
	SingleFluidSolver(const DualMesh * dualMesh);
	~SingleFluidSolver();

	const double updateFluidState() override;
};

