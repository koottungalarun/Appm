#pragma once
#include "AppmSolver.h"
class AppmSolverLeapfrog :
	public AppmSolver
{
public:
	AppmSolverLeapfrog();
	~AppmSolverLeapfrog();

protected:
	void init_maxwell(const double dt) override;
	void update_maxwell(const double dt, const double time) override;


private:
	Eigen::SparseMatrix<double> A;
	Eigen::SparseMatrix<double> C;
	Eigen::SparseMatrix<double> M;
	Eigen::SparseMatrix<double> M_d, M_f;
	Eigen::SparseMatrix<double> Q;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> maxwellSolver;
	Eigen::VectorXd x_m, x_mm1;

};

