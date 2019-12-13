#include "AppmSolverCrankNichol.h"



AppmSolverCrankNichol::AppmSolverCrankNichol()
{
}


AppmSolverCrankNichol::~AppmSolverCrankNichol()
{
}

void AppmSolverCrankNichol::init_maxwell()
{
	// number of degrees of freedom
	const int ndof =
		primalMeshInfo.nFacesInner
		+ primalMeshInfo.nEdgesInner
		+ primalMeshInfo.nVerticesBoundary;
	maxwellState = Eigen::VectorXd::Zero(ndof);

	// Call function of base class
	AppmSolver::init_maxwell();

	// Extend initialization for Maxwell equations



}


void AppmSolverCrankNichol::init_maxwell(const double dt)
{
	// number of degrees of freedom
	const int ndof =
		primalMeshInfo.nFacesInner
		+ primalMeshInfo.nEdgesInner
		+ primalMeshInfo.nVerticesBoundary;
	maxwellState = Eigen::VectorXd::Zero(ndof);

	// Call function of base class
	AppmSolver::init_maxwell(dt);

	// Extend initialization for Maxwell equations

}

void AppmSolverCrankNichol::update_maxwell(const double dt, const double time)
{
	const Eigen::VectorXd prevState = maxwellState;

	// TODO move variables to initialization that are repeatedly used
	const int n = maxwellState.size();
	const int nFacesInterior = primalMeshInfo.nFacesInner;
	const int nEdgesInterior = primalMeshInfo.nEdgesInner;
	const int nEdges = primalMeshInfo.nEdges;
	const int nVerticesTerminal = primalMeshInfo.nVerticesTerminal;
	const int nVerticesBoundary = primalMeshInfo.nVerticesBoundary;
	const int nVertices = primalMeshInfo.nVertices;

	// Discrete Operators
	const Eigen::SparseMatrix<double> M_mu    = hodgeOperatorDualEdgeToPrimalFace();
	const Eigen::SparseMatrix<double> M_eps   = hodgeOperatorPrimalEdgeToDualFace();
	const Eigen::SparseMatrix<double> M_sigma = hodgeOperatorElectricalConductivity();

	const Eigen::SparseMatrix<double> G = primalMesh.get_e2vMap().cast<double>();
	const Eigen::SparseMatrix<double> C = primalMesh.get_f2eMap().cast<double>();

	const Eigen::SparseMatrix<double> X = inclusionOperatorBoundaryVerticesToAllVertices();
	
	// Boundary gradient operator
	const Eigen::SparseMatrix<double> GX = G * X;

	// Operator from inner edges and boundary vertices to all edges
	const int nRowsQ = nEdgesInterior + nVerticesBoundary;
	Eigen::SparseMatrix<double> Q = speye(nEdges, nRowsQ);
	Q.rightCols(GX.cols()) = -GX;

	// Curl operator in interior (inner faces, inner edges)
	const Eigen::SparseMatrix<double> C_inner = C.topLeftCorner(nFacesInterior, nEdgesInterior);
	assert(C_inner.rows() == nFacesInterior);
	assert(C_inner.cols() == nEdgesInterior);
	Eigen::SparseMatrix<double> P(nFacesInterior, nEdgesInterior + nVerticesBoundary);
	P.leftCols(C_inner.cols()) = C_inner;
	
	// Solve system of equations: A*u(k+1) = B*u(k) for u(k+1)
	Eigen::SparseMatrix<double> A(n,n);
	Eigen::SparseMatrix<double> B(n,n);
	Eigen::SparseMatrix<double> A_11, A_12, A_21, A_22;
	A_11 = M_mu;
	A_12 = +0.5 * dt * P;
	A_21 = -0.5 * dt * P.transpose();
	A_22 = Q.transpose() * (M_eps + 0.5 * dt * M_sigma) * Q;

	const Eigen::SparseMatrix<double> Aleft = Eigen::vertcat(A_11, A_21);
	const Eigen::SparseMatrix<double> Aright = Eigen::vertcat(A_12, A_22);
	A.leftCols(Aleft.cols()) = Aleft;
	A.rightCols(Aright.cols()) = Aright;

	Eigen::SparseMatrix<double> B_11, B_12, B_21, B_22;
	B_11 = M_mu;
	B_12 = -0.5 * dt * P;
	B_21 = +0.5 * dt * P.transpose();
	B_22 = Q.transpose() * (M_eps - 0.5 * dt * M_sigma) * Q;

	const Eigen::SparseMatrix<double> Bleft  = Eigen::vertcat(B_11, B_21);
	const Eigen::SparseMatrix<double> Bright = Eigen::vertcat(B_12, B_22);
	B.leftCols(Bleft.cols()) = Bleft;
	B.rightCols(Bright.cols()) = Bright;

	const int nDirichlet = nVerticesTerminal; // number of fixed values by boundary condition
	const int nFree = n - nDirichlet; // number of free unknowns 

	const Eigen::SparseMatrix<double> M_f = A.topLeftCorner(nFree, nFree);
	const Eigen::SparseMatrix<double> M_d = A.rightCols(nDirichlet);

	// Solve system of equations for free unknowns
	Eigen::SparseLU<Eigen::SparseMatrix<double>> maxwellSolver(M_f);
	const Eigen::VectorXd x_d   = electricPotentialTerminals(time); // voltage boundary condition on terminals at new timestep
	const Eigen::VectorXd rhs   = B * prevState - M_d * x_d;        // get right hand side of equation
	const Eigen::VectorXd rhs_f = rhs.topRows(M_f.rows());          // discard equations to terminal vertices
	Eigen::VectorXd x_f(nFree);
	x_f = maxwellSolver.solve(rhs_f);       // solve system of equations

	// assign solution to state vector
	maxwellState.topRows(nFree) = x_f; 
	maxwellState.bottomRows(nDirichlet) = x_d; 

	// map data to vectors
	const int nH = primalMeshInfo.nFacesInner;
	H_h.topRows(nH) = maxwellState.topRows(nH);
	const int nx = Q.cols();
	E_h = Q * maxwellState.bottomRows(nx);
}

