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
	primalMesh.init();
	//dualMesh = DualMesh(primalMesh);

	writeMesh();
}

void AppmSolver::writeMesh()
{
	primalMesh.writeToFile();
	primalMesh.writeXdmf();
	//dualMesh.writeToFile();
}
