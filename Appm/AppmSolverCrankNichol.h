#pragma once
#include "AppmSolver.h"
class AppmSolverCrankNichol :
	public AppmSolver
{
public:
	AppmSolverCrankNichol();
	~AppmSolverCrankNichol();

protected:
	void init_maxwell() override;
	void init_maxwell(const double dt) override;
	void update_maxwell(const double dt, const double time) override;

private:
	// State vector for Maxwell equations. 
	// For Crank-Nicholson scheme, we have state = [H_h, x_h] with E_h = Q*x_h.
	Eigen::VectorXd maxwellState; 


};

