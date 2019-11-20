#include "AppmSolver.h"



AppmSolver::AppmSolver()
{
}


AppmSolver::~AppmSolver()
{
}

void AppmSolver::run()
{
	primalMesh = PrimalMesh();
	dualMesh = DualMesh(primalMesh);

	writeMesh();
}

void AppmSolver::writeMesh()
{
	primalMesh.writeToFile();
	primalMesh.writeXdmf();
	dualMesh.writeToFile();
}
