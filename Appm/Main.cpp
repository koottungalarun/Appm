#include "Main.h"

int main() {
	std::cout << "***********************" << std::endl;
	std::cout << "*    APPM             *" << std::endl;
	std::cout << "***********************" << std::endl;
	Main main;
	main.run();
	std::cout << "TERMINATED" << std::endl;
	return EXIT_SUCCESS;
}


Main::Main()
{
}


Main::~Main()
{
}

void Main::run()
{
	FluidPrimitiveState primitive_left;
	primitive_left.p = 1.0;
	primitive_left.rho = 1.0;
	primitive_left.u = Eigen::Vector3d(0, 0, 0);

	FluidPrimitiveState primitive_right;
	primitive_right.p = 0.1;
	primitive_right.rho = 0.125;
	primitive_right.u = Eigen::Vector3d(0, 0, 0);

	FluidState qL(primitive_left);
	FluidState qR(primitive_right);

	Eigen::VectorXd qLvec = qL.getVecState();
	Eigen::VectorXd qRvec = qR.getVecState();
	
	Eigen::Vector3d fn(0, 0, 1);
	const double dx = 1.0;
	double dt_loc = 0;
	Eigen::VectorXd flux = Numerics::fluidFlux_rusanov(qLvec, qRvec, fn, dx, dt_loc);

	std::cout << "flux: " << flux.transpose() << std::endl;
	std::cout << "dt:   " << dt_loc << std::endl;

	AppmSolver * appm = new AppmSolverCrankNichol();
	appm->run();
	delete appm;
}
