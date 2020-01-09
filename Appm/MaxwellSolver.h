#pragma once

#include "PrimalMesh.h"
#include "DualMesh.h"
#include "H5Writer.h"

#include <Eigen/Dense>

class MaxwellSolver
{
public:
	MaxwellSolver();
	MaxwellSolver(const PrimalMesh * primal, const DualMesh * dual);
	~MaxwellSolver();

	bool isMaxwellCurrentSource = false;


	virtual void updateMaxwellState(const double dt, const double time);

	void writeStates(H5Writer & writer) const;

	const Eigen::VectorXd & getBstate() const;

	void setUniformMagneticFluxField(const Eigen::Vector3d & fieldVector);
	void setAzimuthalMagneticFluxField();
	void setTorusCurrent(const double x1, const double x2, const double z1, const double z2);

	Eigen::VectorXd electricPotentialTerminals(const double time);


protected:
	const PrimalMesh * primal = nullptr;
	const DualMesh * dual = nullptr;

	Mesh::MeshInfo primalMeshInfo;
	Mesh::MeshInfo dualMeshInfo;

	Eigen::VectorXd maxwellState;

	Eigen::VectorXd B_h, E_h, H_h, J_h;

	Eigen::SparseMatrix<double> get_Mnu();
	Eigen::SparseMatrix<double> get_Meps();
	Eigen::SparseMatrix<double> get_Msigma();
	Eigen::SparseMatrix<double> get_bndryInclOp();

	Eigen::SparseMatrix<int> setupOperatorQ();
	Eigen::SparseMatrix<double> setupOperatorMeps();
	Eigen::SparseMatrix<double> setupOperatorMnu();


private:
	void init();

	Eigen::SparseMatrix<double> Mnu;
	Eigen::SparseMatrix<double> Meps;
	Eigen::SparseMatrix<double> Msigma;
	Eigen::SparseMatrix<double> bndryInclOp;


};

