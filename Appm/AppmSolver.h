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

	XdmfGrid getOutputPrimalEdgeGrid(const int iteration, const double time, const std::string & dataFilename);
	XdmfGrid getOutputPrimalSurfaceGrid(const int iteration, const double time, const std::string & dataFilename);
	
	XdmfGrid getOutputDualEdgeGrid(const int iteration, const double time, const std::string & dataFilename);
	XdmfGrid getOutputDualSurfaceGrid(const int iteration, const double time, const std::string & dataFilename);
};

