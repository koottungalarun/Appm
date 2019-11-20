#pragma once

#include "PrimalMesh.h"
#include "DualMesh.h"


class AppmSolver
{
public:
	AppmSolver();
	~AppmSolver();

	void run();

private:
	PrimalMesh primalMesh;
	DualMesh dualMesh;


	void writeMesh();
};

