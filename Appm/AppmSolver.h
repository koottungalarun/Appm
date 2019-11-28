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

	Eigen::VectorXd bvec, dvec, evec, hvec, jvec;

	std::vector<double> timeStamps;

	void writeXdmf();
	void writeXdmfPrimalMesh();
	void writeXdmfDualMesh();

	void writeOutput(const int iteration, const double time);
};

