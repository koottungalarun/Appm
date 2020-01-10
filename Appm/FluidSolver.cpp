#include "FluidSolver.h"



FluidSolver::FluidSolver()
{
}


FluidSolver::FluidSolver(const DualMesh * mesh)
{
	this->mesh = mesh;
}

FluidSolver::~FluidSolver()
{
}

const double FluidSolver::updateFluidState()
{
	const int nDualFaces = mesh->getNumberOfFaces();
	Eigen::VectorXd dt_local(nDualFaces);
	Eigen::MatrixXd dq = Eigen::MatrixXd::Zero(fluidStates.rows(), fluidStates.cols());

	// Compute fluid fluxes at faces
	for (int i = 0; i < nDualFaces; i++) {
		const Face * face = mesh->getFace(i);
		const Eigen::Vector3d faceNormal = face->getNormal();
		const double faceArea = face->getArea();

		if (face->isBoundary()) {
			const std::vector<Cell*> faceCells = face->getCellList();
			assert(faceCells.size() == 1);
			Cell * cell = faceCells[0];
			const int idxC = cell->getIndex();
			const int orientation = (face->getCenter() - cell->getCenter()).dot(faceNormal) > 0 ? 1 : -1;

			Eigen::VectorXd qL, qR;
			// TODO open boundary conditions
			if (orientation > 0) {
				qL = fluidStates.col(idxC);
				qR = qL;
				qR.segment(1, 3) *= -1;
			}
			else {
				qR = fluidStates.col(idxC);
				qL = qR;
				qL.segment(1, 3) *= -1;
			}
			const Eigen::Vector3d cc = cell->getCenter();
			const Eigen::Vector3d fc = face->getCenter();
			const double dx = std::abs((fc - cc).dot(faceNormal));
			double dt_loc = 0;
			//const Eigen::VectorXd flux = Numerics::fluidFlux_rusanov(qL, qR, faceNormal, dx, dt_loc);
			const Eigen::VectorXd flux = getRusanovFlux(qL, qR, faceNormal, dx, dt_loc);
			assert(dt_loc > 0);
			dt_local(i) = dt_loc;
			dq.col(idxC) += orientation * faceArea * flux;
		}
		else {
			// Get left and right cell of this face
			std::vector<Cell*> faceCells = face->getCellList();
			assert(faceCells.size() == 2);
			const Cell * cell = faceCells[0];
			const int orientation = (face->getCenter() - cell->getCenter()).dot(faceNormal) > 0 ? 1 : -1;
			Cell * leftCell = nullptr;
			Cell * rightCell = nullptr;
			if (orientation > 0) {
				leftCell = faceCells[0];
				rightCell = faceCells[1];
			}
			else {
				leftCell = faceCells[1];
				rightCell = faceCells[0];
			}
			const int idxL = leftCell->getIndex();
			const int idxR = rightCell->getIndex();

			// distance between cell center and face center
			const Eigen::Vector3d fc = face->getCenter();
			const Eigen::Vector3d ccL = leftCell->getCenter();
			const Eigen::Vector3d ccR = rightCell->getCenter();
			const double dxL = std::abs((fc - ccL).dot(faceNormal));
			const double dxR = std::abs((fc - ccR).dot(faceNormal));
			const double dx = std::min(dxL, dxR); // minimum distance between cell center and face center
			assert(dx > 0);
			const Eigen::VectorXd qL = fluidStates.col(idxL);
			const Eigen::VectorXd qR = fluidStates.col(idxR);

			double dt_loc = 0;
			//const Eigen::VectorXd faceFlux = Numerics::fluidFlux_rusanov(qL, qR, faceNormal, dx, dt_loc);
			const Eigen::VectorXd faceFlux = getRusanovFlux(qL, qR, faceNormal, dx, dt_loc);
			assert(dt_loc > 0);
			dt_local(i) = dt_loc;

			dq.col(idxL) += faceFlux * faceArea;
			dq.col(idxR) -= faceFlux * faceArea;
		}
	}
	bool allPositive = (dt_local.array() > 0.).all();
	assert(allPositive);

	// get global timestep size
	const double dt = dt_local.minCoeff();
	assert(dt > 0);

	for (int i = 0; i < mesh->getNumberOfCells(); i++) {
		dq.col(i) *= 1. / (mesh->getCell(i)->getVolume());
	}

	fluidStates -= dt * dq;

	return dt;
}