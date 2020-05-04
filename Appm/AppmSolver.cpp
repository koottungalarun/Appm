#include "AppmSolver.h"



AppmSolver::AppmSolver() 
	: AppmSolver(PrimalMesh::PrimalMeshParams())
{
}

AppmSolver::AppmSolver(const PrimalMesh::PrimalMeshParams & primalMeshParams)
{
	isStateWrittenToOutput = true;
	readParameters("AppmSolverParams.txt");
	init_meshes(primalMeshParams);  // Initialize primal and dual meshes
	if (primalMesh.getNumberOfCells() == 0) {
		return;
	}
	std::cout << "Dual mesh has " << dualMesh.getNumberOfVertices() << " vertices" << std::endl;


	init_multiFluid("particleParameters.txt");
	
	switch (initType) {
	case 1:
	{
		const double zRef = 0.;
		init_SodShockTube(zRef);
	}
	break;

	case 2:
	{
		double p = 1.0;
		double n = 1.0;
		double u = 0.0;
		init_Uniformly(n, p, u);
	}
		break;

	case 3:
	{
		const Eigen::Vector3d refPos = Eigen::Vector3d(0, 0, 0);
		const double radius = 0.2;
		init_Explosion(refPos, radius);
	}
	break;

	default:
		exit(-1);
	}

	const int nFluids = this->getNFluids();
	const int nFaces = dualMesh.getNumberOfFaces();
	faceFluxes = Eigen::MatrixXd::Zero(5 * nFluids, nFaces);
	faceFluxesImExRusanov = Eigen::MatrixXd::Zero(nFaces, nFluids);


	init_maxwellStates();
	B_vertex = Eigen::Matrix3Xd::Zero(3, primalMesh.getNumberOfVertices());
	init_RaviartThomasInterpolation();


	if (isMaxwellEnabled) {
		isWriteBfield = true;
		isWriteEfield = true;
	}

	// Check implementation of Ohm's law in implicit-consistent formulation
	// TODO: to be removed after testing
	Eigen::SparseMatrix<double> Msigma;
	double dt = timestepSize;

	set_Bfield_azimuthal();
	interpolateMagneticFluxToPrimalVertices();

	get_Msigma_consistent(dt, Msigma, J_h_aux);

	set_Efield_uniform(Eigen::Vector3d::UnitZ());
	J_h = Msigma * E_h + J_h_aux;


}


AppmSolver::~AppmSolver()
{
}

void AppmSolver::run()
{
	if (primalMesh.getNumberOfCells() == 0) {
		return;
	}

	double dt = timestepSize;
	double dt_previous = dt;
	double time = 0;
	int iteration = 0;


	// For testing of RT interpolation:
	//setAzimuthalMagneticFluxField();
	//setUniformMagneticFluxField(Eigen::Vector3d(1,1,1));
	//interpolateMagneticFluxToPrimalVertices();

	// initialize current flow
	//this->isMaxwellCurrentSource = false;
	//maxwellSolver->isMaxwellCurrentSource = isMaxwellCurrentSource;
	//if (isMaxwellCurrentSource) {
	//	const double x1 = -0.5;
	//	const double x2 = 0.5;
	//	const double z1 = 0.24;
	//	const double z2 = 0.76;
	//	maxwellSolver->setTorusCurrent(x1, x2, z1, z2);
	//}



	// Number of primal vertices on terminals
	const int nPrimalTerminalVertices = primalMesh.getMeshInfo().nVerticesTerminal;
	assert(nPrimalTerminalVertices > 0);

	// Number of primal edges
	const int nEdges = primalMesh.getNumberOfEdges();


	const int nFluids = this->getNFluids();
	const int nFaces = dualMesh.getNumberOfFaces();
	//faceFluxes = Eigen::MatrixXd::Zero(5*nFluids, nFaces);
	//faceFluxesImExRusanov = Eigen::MatrixXd::Zero(nFaces, nFluids);
	
	// write initial data to file (iteration = 0, time = 0)
	writeOutput(iteration, time);

	const std::string stopFilename = "stop.txt";
	std::ofstream(stopFilename) << 0 << std::endl;

	auto timer_startAppmSolver = std::chrono::high_resolution_clock::now();


	std::vector<double> currentBCValue;
	/*
	* Time integration loop 
	*/
	while (iteration < maxIterations && time < maxTime) {
		std::cout << "Iteration " << iteration << ",\t time = " << time << std::endl;

		dt_previous = dt;

		fluidFluxes.setZero();
		fluidSources.setZero();
		faceFluxes.setZero();
		faceFluxesImExRusanov.setZero();

		// Determine timestep
		if (isFluidEnabled) {
			dt = getNextFluidTimestepSize();
		}


		// Maxwell equations
		if (isMaxwellEnabled) {

			// New state vector, solution of the implicit system of equation that we are setting up in the next lines
			Eigen::VectorXd x(maxwellState.size());

			Eigen::SparseMatrix<double> Msigma;
			get_Msigma_consistent(dt , Msigma, J_h_aux);

			Eigen::SparseMatrix<double> Msigma_inner;
			Msigma_inner = Msigma.topLeftCorner(nEdges, nEdges);
			
			// System matrix on left hand side, i.e., M*x = rhs, to be solved for x
			Eigen::SparseMatrix<double> M;
			assert(M1.size() > 0);
			assert(M1.nonZeros() > 0);
			assert(M2.size() > 0);
			assert(M2.nonZeros() > 0);
			M = M1 + pow(dt, 2) * M2 + dt * Q.transpose() * Msigma_inner * Q;
			M.makeCompressed();
			//Eigen::sparseMatrixToFile(M, "M.dat");

			// The system of equations has fixed and free values; 
			// - fixed values: electric potential at terminals (Dirichlet boundary condition)
			// - free  values: electric potential at non-terminal vertices, and electric voltages at non-boundary edges
			int nDirichlet = nPrimalTerminalVertices;
			//if (maxwellSolverBCType == MaxwellSolverBCType::CURRENT_BC) {
			//	nDirichlet = nPrimalTerminalVertices / 2;
			//}
			const int nFree = maxwellState.size() - nDirichlet;

			// The vector of degrees of freedom (DoF) is sorted such that free values are in front of fixed values:
			// x = [freeValues, fixedValues]
			// Therefore, the system of free DoF is given as: (without considering proper array sizes)
			// M_free * x_free = -M_fixed * x_fixed + rhs

			// Dirichlet conditions
			Eigen::SparseMatrix<double> Md = M.rightCols(nDirichlet);
			Eigen::VectorXd xd = Eigen::VectorXd::Zero(nDirichlet);
			
			//switch (maxwellSolverBCType) {
			//	case MaxwellSolverBCType::VOLTAGE_BC:
			//	{
					assert(nDirichlet % 2 == 0);
					const double t0 = 1;
					const double t1 = 11;
					const double tscale = 0.2;
					xd.topRows(nDirichlet / 2) = terminalVoltageBC_sideA(time, t0, t1, tscale) * Eigen::VectorXd::Ones(nDirichlet / 2);
					xd.bottomRows(nDirichlet / 2) = terminalVoltageBC_sideB(time) * Eigen::VectorXd::Ones(nDirichlet / 2);
					//break;
			//	}

			//	case MaxwellSolverBCType::CURRENT_BC:
			//	{
			//		xd.bottomRows(nDirichlet) = terminalVoltageBC_sideB(time) * Eigen::VectorXd::Ones(nDirichlet);
			//		break;
			//	}

			//	default:
			//	{
			//		std::cout << "Maxwell Solver Boundary Type not implemented" << std::endl;
			//		assert(false);
			//	}
			//}


			Eigen::VectorXd src = Eigen::VectorXd::Zero(maxwellState.size());

			// Vector on right hand side
			Eigen::VectorXd rhs(x.size());
			rhs.setZero();

			// Setup of rhs vector
			double dt_ratio = dt / dt_previous;
			rhs += M1 * (1 + dt_ratio) * maxwellState 
				+ dt * Q.transpose() * Msigma_inner * Q * maxwellState 
				- dt_ratio * M1 * maxwellStatePrevious;

			// TODO: current source in Ampere equation
			//if (maxwellSolverBCType == MaxwellSolverBCType::CURRENT_BC) {
			//	J_h_previous = J_h;
			//	
			//	const double currentValue = currentDensityBC(time);
			//	currentBCValue.push_back(currentValue);
			//	const Eigen::VectorXd zUnitVec = Eigen::Vector3d::UnitZ();
			//	Eigen::VectorXd dualFacesInZDirection(dualMesh.getNumberOfFaces());
			//	dualFacesInZDirection.setZero();
			//	for (int i = 0; i < dualFacesInZDirection.size(); i++) {
			//		const Face * face = dualMesh.getFace(i);
			//		const Eigen::Vector3d fc = face->getCenter();
			//		if (true || fc.segment(0, 2).norm() < 0.35) {
			//			const Eigen::Vector3d fn = face->getNormal();
			//			const double fA = face->getArea();
			//			dualFacesInZDirection(i) = fn.dot(zUnitVec) * fA;
			//		}
			//	}
			//	J_h = currentValue * dualFacesInZDirection;
			//	const Eigen::VectorXd deltaJ = dt * Q.transpose() * (J_h - J_h_previous).segment(0, nEdges);
			//	rhs -= deltaJ;
			//}

			// TODO: this is the implicit-consistent face current due to fluid momentum flux
			// update data vector for J_h_aux
			//for (int i = 0; i < nEdges; i++) {
			//	assert(false);
			//	assert(getNFluids() == 1);
			//	const Face * face = dualMesh.getFace(i);
			//	//assert(face->hasFluidCells());
			//	const std::vector<Cell*> adjacientCells = face->getCellList();
			//	assert(adjacientCells.size() == 2);
			//	int nAdjacientFluidCells = 0;
			//	for (auto cell : adjacientCells) {
			//		//assert(cell->getFluidType() == Cell::FluidType::FLUID);
			//		if (cell->getFluidType() == Cell::FluidType::FLUID) {
			//			nAdjacientFluidCells++;
			//		}
			//	}
			//	if (nAdjacientFluidCells < 2) { // skip faces that are adjacient to a solid cell
			//		continue;
			//	}

			//	assert(massFluxScheme == MassFluxScheme::IMPLICIT_EXPLICIT);
			//	const double fA = face->getArea();
			//	const int faceIdx = face->getIndex();
			//	double value = 0;

			//	// Loop over all species
			//	for (int fluidIdx = 0; fluidIdx < getNFluids(); fluidIdx++) {
			//		const int q = particleParams[fluidIdx].electricCharge; 
			//		const Eigen::Vector3d implicitFlux = getRusanovFluxImEx(faceIdx, fluidIdx, dt);
			//		const double implicitMassFlux = implicitFlux(0); // first component is the mass flux
			//		value += q * implicitMassFlux;
			//	}
			//	value *= fA;
			//	J_h_aux(i) = value;
			//}

			

			rhs -= dt * Q.transpose() * (J_h_aux - J_h_aux_mm1).segment(0, nEdges);

			rhs -= Md * xd;

			// Vector of free coefficients
			Eigen::VectorXd xf(nFree);

			// Load vector for free coefficients
			Eigen::VectorXd rhsFree = rhs.topRows(nFree);

			// Matrix of free coefficients
			Eigen::SparseMatrix<double> Mf = M.topLeftCorner(nFree, nFree);
			Mf.makeCompressed();
			

			// Solve system
			std::cout << "Setup Maxwell solver" << std::endl;
			Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> maxwellSolver;

			auto timer_start = std::chrono::high_resolution_clock::now();
			maxwellSolver.compute(Mf);			
			auto timer_afterCompute = std::chrono::high_resolution_clock::now();
			auto delta_afterCompute = std::chrono::duration<double>(timer_afterCompute - timer_start);
			std::cout << "Solver time for compute step: " << delta_afterCompute.count() << std::endl;
			if (maxwellSolver.info() != Eigen::Success) {
				std::cout << "Maxwell solver setup failed" << std::endl;
			}
			auto timer_startSolve = std::chrono::high_resolution_clock::now();
			xf = maxwellSolver.solve(rhsFree);
			auto timer_endSolve = std::chrono::high_resolution_clock::now();
			std::cout << "Maxwell solver finished" << std::endl;
			auto delta_solve = std::chrono::duration<double>(timer_endSolve - timer_startSolve);
			std::cout << "Solver time for solve step: " << delta_solve.count() << std::endl;
			if (maxwellSolver.info() != Eigen::Success) {
				std::cout << "Maxwell solver solving failed" << std::endl;
			}
			assert(maxwellSolver.info() == Eigen::Success);

			std::cout << "#iterations:     " << maxwellSolver.iterations() << std::endl;
			std::cout << "estimated error: " << maxwellSolver.error() << std::endl;

			// assemble new state vector
			x.topRows(nFree) = xf;
			x.bottomRows(nDirichlet) = xd;
			//std::ofstream("x.dat") << x << std::endl;

			// Update state vectors for discrete states
			E_h = Q * x;
			B_h -= dt * C * E_h;

			// Get electric current due to Ohms law (see implicit and consistent formulation of electric current)
			//assert(Msigma.rows() == nEdges);
			//assert(Msigma.cols() == E_h.size());
			//J_h.segment(0, nEdges) = Msigma * E_h;
			J_h = Msigma * E_h + J_h_aux;

			//std::ofstream("E_h.dat") << E_h << std::endl;

			// Set new state vector
			maxwellStatePrevious = maxwellState;
			maxwellState = x;

			J_h_aux_mm1 = J_h_aux;


			// Get electric field at cell center (interpolated from primal edge values)
			E_cc = getEfieldAtCellCenter();

			// Interpolation of B-field to dual cell centers
			interpolateMagneticFluxToPrimalVertices();
		}


		// Fluid equations
		if (isFluidEnabled) {
			const int nFluids = this->getNFluids();
			const int nCells = dualMesh.getNumberOfCells();
			const int nFaces = dualMesh.getNumberOfFaces();
			
			// Calculate flux at each dual face
			for (int fidx = 0; fidx < nFaces; fidx++) {
				const Face * face = dualMesh.getFace(fidx);
				const Eigen::Vector3d faceNormal = face->getNormal();
				assert(std::fabs(faceNormal.norm() -  1) < 16*std::numeric_limits<double>::epsilon());
				const double faceArea = face->getArea();
				const std::vector<Cell*> faceCells = face->getCellList();
				assert(faceCells.size() >= 1);

				// skip faces that have no adjacient fluid cell
				if (!face->hasFluidCells()) {
					continue;
				}

				for (int fluidIdx = 0; fluidIdx < nFluids; fluidIdx++) {
					Eigen::Vector3d qL, qR;  // left and right cell state at this face
					double s = 0;            // wavespeed at face
					const int orientation = getOrientation(faceCells[0], face);

					if (fidx == faceIdxRef) {
						std::cout << "face idx ref: " << faceIdxRef << std::endl;
						std::cout << "Face center: " << face->getCenter().transpose() << std::endl;
						std::cout << "Face normal: " << faceNormal.transpose() << std::endl;
						std::cout << "orientation: " << orientation << std::endl;
					}

					// Explicit Rusanov scheme
					Eigen::Vector3d flux = Eigen::Vector3d::Zero();
					bool isCollinearZ = faceNormal.cross(Eigen::Vector3d::UnitZ()).norm() < (std::numeric_limits<double>::epsilon() * 128);
					const Face::FluidType faceFluidType = face->getFluidType();
					//if (true || isCollinearZ) {
						//flux = getRusanovFluxExplicit(fidx, fluidIdx);
					flux = getRusanovFluxImEx(fidx, fluidIdx, dt);
					//}

					// 3D face flux data vector 
					Eigen::VectorXd faceFlux3d(5);
					faceFlux3d(0) = flux(0);
					faceFlux3d.segment(1, 3) = flux(1) * faceNormal;
					faceFlux3d(4) = flux(2);
					
					// store data
					faceFluxes.col(fidx) = faceFlux3d; 
					if (fidx == faceIdxRef) {
						std::cout << "face flux 3d: " << faceFlux3d.transpose() << std::endl;
					}

					// multiply flux by face area
					faceFlux3d *= faceArea;

					if (face->getFluidType() == Face::FluidType::INTERIOR) {
						assert(faceCells.size() == 2);
						const int idx0 = faceCells[0]->getIndex();
						const int idx1 = faceCells[1]->getIndex();
						fluidFluxes.col(idx0).segment(5 * fluidIdx, 5) += orientation * faceFlux3d;
						fluidFluxes.col(idx1).segment(5 * fluidIdx, 5) -= orientation * faceFlux3d;
					}
					if (face->getFluidType() == Face::FluidType::OPENING) {
						assert(faceCells.size() == 1);
						const Cell * cell = faceCells[0];
						assert(cell->getFluidType() == Cell::FluidType::FLUID);
						fluidFluxes.col(cell->getIndex()).segment(5 * fluidIdx, 5) += orientation * faceFlux3d;
					}
					if (face->getFluidType() == Face::FluidType::WALL) {
						assert(faceCells.size() == 1 || faceCells.size() == 2);
						const Cell * cell = faceCells[0];
						if (cell->getFluidType() == Cell::FluidType::FLUID) {
							const int idx = cell->getIndex();
							fluidFluxes.col(idx).segment(5 * fluidIdx, 5) += orientation * faceFlux3d;
						}
						else {
							assert(false);
						}
					}
				}
			}

			// Set source term
			for (int cIdx = 0; cIdx < nCells; cIdx++) {
				const Cell * cell = dualMesh.getCell(cIdx);

				// Skip cell that ore not of type Fluid
				if (cell->getFluidType() != Cell::FluidType::FLUID) {
					continue;
				}
				// For all fluids ... 
				for (int fluidIdx = 0; fluidIdx < nFluids; fluidIdx++) {
					Eigen::VectorXd srcLocal(5);
					srcLocal.setZero();

					//// Define geometric position of source region: ball of radius r around reference position
					//const Eigen::Vector3d cc = cell->getCenter();
					//const Eigen::Vector3d srcRefPos(0, 0, 0);
					//const double srcRadius = 0.2;
					//if ((cc - srcRefPos).norm() < srcRadius) {
					//	srcLocal(0) = 0; // mass source
					//	srcLocal.segment(1, 3) = Eigen::Vector3d::Zero(); // momentum source
					//	srcLocal(4) = 1; // energy source
					//}

					// TODO: Lorentz force (electrostatic)
					// const Eigen::Vector3d Efield_atCellCenter = getEfieldAtDualCellCenter(cIdx);


					fluidSources.col(cIdx).segment(5 * fluidIdx, 5) = srcLocal;
				}
			}
			
			// Update to next timestep: U(m+1) = U(m) - dt / volume * sum(fluxes)
			for (int i = 0; i < nCells; i++) {
				const Cell * cell = dualMesh.getCell(i);
				double cellVolume = cell->getVolume();
				fluidStates_new.col(i) = fluidStates.col(i) - dt / cellVolume * fluidFluxes.col(i) + dt * fluidSources.col(i);
			}
			fluidStates = fluidStates_new; // update data storage to new timestep
			fluidStates_new.setZero();     // clear auxiliary data storage
		}
		std::cout << "dt = " << dt << std::endl;


		iteration++;
		time += dt;
		writeOutput(iteration, time);

		int stopValue = 0;
		std::ifstream(stopFilename) >> stopValue;
		if (stopValue > 0) {
			std::cout << "Stop because of value set in stop-file (" << stopFilename << ")" << std::endl;
			std::cout << "Stop value: " << stopValue << std::endl;
			break;
		}
		else {
			std::cout << "Continue; value in stop-file is non-positive" << std::endl;
		}
	}
	std::cout << "Final time:      " << time << std::endl;
	std::cout << "Final iteration: " << iteration << std::endl;

	auto timer_endAppmSolver = std::chrono::high_resolution_clock::now();
	auto delta_appmSolver = std::chrono::duration<double>(timer_endAppmSolver - timer_startAppmSolver);
	std::cout << "Elapsed time for APPM solver: " << delta_appmSolver.count() << std::endl;

	if (currentBCValue.size() > 0) {
		std::string filename = "currentBC.dat";
		std::cout << "Write data to file: " << filename << std::endl;
		std::ofstream currentBCoutput(filename);
		for (auto value : currentBCValue) {
			currentBCoutput << value << std::endl;
		}
	}

	// test_raviartThomas();

	// Use Paraview (version 5.6.0) to visualize.
	// Light data is given in XDMF files (version 3). 
	// Heavy data is stored in HDF5 files.

	// Note that polyhedrons have been added to XDMF in version 3; version 2 does not support polyhedrons.
	// Moreover, Paraview has 3 readers: XDMF Reader, Xdmf3ReaderS, Xdmf3ReaderT.
	// Xdmf3Reader3: can     read polyhedrons, but not grid-of-grids (i.e., grids with GridType=Tree).
	// XDMF Reader:  can not read polyhedrons, but     grid-of-grids.
	// Therefore, the mesh data for vertices, edges, and faces, are separated from the volume data.
	writeXdmf();
	writeXdmfDualVolume();
	//writeXdmfDualFaceFluxes();
}

const int AppmSolver::getNFluids() const
{
	return particleParams.size();
}

void AppmSolver::init_maxwellStates()
{
	std::cout << "Initialize Maxwell states" << std::endl;


	E_h = Eigen::VectorXd::Zero(primalMesh.getNumberOfEdges());
	B_h = Eigen::VectorXd::Zero(primalMesh.getNumberOfFaces());
	J_h = Eigen::VectorXd::Zero(dualMesh.getNumberOfFaces());
	J_h_previous = J_h;

	const Eigen::VectorXi vertexTypes = primalMesh.getVertexTypes();
	const Eigen::VectorXi edgeTypes = primalMesh.getEdgeTypes();
	const Eigen::VectorXi faceTypes = primalMesh.getFaceTypes();

	// Number of primal vertices on surface boundary (terminals & free floating)
	const int nPvb = 
		(vertexTypes.array() == static_cast<int>(Vertex::Type::Boundary)).count() +
		(vertexTypes.array() == static_cast<int>(Vertex::Type::Terminal)).count();

	// Number of primal interior edges (not on boundary)
	const int nPei =
		(edgeTypes.array() == static_cast<int>(Edge::Type::Interior)).count() +
		(edgeTypes.array() == static_cast<int>(Edge::Type::InteriorToBoundary)).count();

	// Number of primal interior faces 
	const int nPfi = (faceTypes.array() == 0).count();

	// Number of degrees of freedom for Maxwell's equation
	const int nDof = nPei + nPvb;
	this->maxwellState = Eigen::VectorXd::Zero(nDof);
	this->maxwellStatePrevious = Eigen::VectorXd::Zero(nDof);

	this->Q = getBoundaryGradientInnerInclusionOperator();
	assert(Q.nonZeros() > 0);
	assert(Q.rows() == primalMesh.getNumberOfEdges());
	assert(Q.cols() == nDof);

	// Material laws
	this->Meps = getElectricPermittivityOperator();
	assert(Meps.nonZeros() > 0);

	this->Mnu = getMagneticPermeabilityOperator();
	assert(Mnu.nonZeros() > 0);

	// Curl operator on all primal faces and all primal edges
	this->C = primalMesh.get_f2eMap().cast<double>();
	assert(C.nonZeros() > 0);

	// Curl operator on inner primal faces and inner primal edges
	const Eigen::SparseMatrix<double> C_inner = C.topLeftCorner(nPfi, nPei);

	Eigen::SparseMatrix<double> P(nPfi, nDof);
	assert(P.cols() == nPei + nPvb);
	P.leftCols(nPei) = C_inner;

	M1 = lambdaSquare * Q.transpose() * Meps * Q;
	M1.makeCompressed();

	const Eigen::SparseMatrix<double> Mnu_inner = Mnu.topLeftCorner(nPfi, nPfi);
	M2 = P.transpose() * Mnu_inner * P;
	M2.makeCompressed();


	// Initialize matrix for interpolating electric field from primal edges to dual cell centers
	initPerotInterpolationMatrix();
	
	// Electric field at cell centers
	E_cc = Eigen::Matrix3Xd::Zero(3, dualMesh.getNumberOfCells());

	// Initialize matrix that maps electric field on primal edges to electric current on dual faces; see Lorentz force in Fluid equations
	//initMsigma();


	const int nFaces = dualMesh.getNumberOfFaces();

	J_h_aux = Eigen::VectorXd(nFaces);
	J_h_aux.setZero();

	J_h_aux_mm1 = Eigen::VectorXd(nFaces);
	J_h_aux_mm1.setZero();

}

Eigen::SparseMatrix<double> AppmSolver::getBoundaryGradientInnerInclusionOperator()
{
	const Eigen::VectorXi vertexTypes = primalMesh.getVertexTypes();
	const int nPvi = (vertexTypes.array() == static_cast<int>(Vertex::Type::Inner)).count();
	const int nPvb =
		(vertexTypes.array() == static_cast<int>(Vertex::Type::Boundary)).count() +
		(vertexTypes.array() == static_cast<int>(Vertex::Type::Terminal)).count();

	Eigen::SparseMatrix<int> G = primalMesh.get_e2vMap();

	// Boundary vertex inclusion operator 
	Eigen::SparseMatrix<int> X(nPvi + nPvb, nPvb);
	typedef Eigen::Triplet<int> T;
	std::vector<T> triplets;
	for (int i = 0; i < nPvb; i++) {
		triplets.push_back(T(nPvi + i, i, 1));
	}
	X.setFromTriplets(triplets.begin(), triplets.end());
	X.makeCompressed();


	// Inner edges inclusion operator
	const Eigen::VectorXi edgeTypes = primalMesh.getEdgeTypes();
	const int nPe = primalMesh.getNumberOfEdges();
	const int nPei =
		(edgeTypes.array() == static_cast<int>(Edge::Type::Interior)).count() +
		(edgeTypes.array() == static_cast<int>(Edge::Type::InteriorToBoundary)).count();

	Eigen::SparseMatrix<int> innerEdgeInclusionOperator(nPe, nPei);
	triplets = std::vector<T>();
	for (int i = 0; i < nPei; i++) {
		triplets.push_back(T(i, i, 1));
	}
	innerEdgeInclusionOperator.setFromTriplets(triplets.begin(), triplets.end());
	innerEdgeInclusionOperator.makeCompressed();

	// Negative gradient on boundary vertices
	Eigen::SparseMatrix<int> negGX = -G * X;

	// Q = [id, -grad] = inner edges inclusion & negative gradient on boundary vertices operator
	Eigen::SparseMatrix<double> Q(nPe, nPei + nPvb);
	assert(Q.rows() == innerEdgeInclusionOperator.rows());
	assert(Q.rows() == negGX.rows());
	assert(Q.cols() == innerEdgeInclusionOperator.cols() + negGX.cols());
	Q.leftCols(innerEdgeInclusionOperator.cols()) = innerEdgeInclusionOperator.cast<double>();
	Q.rightCols(negGX.cols()) = negGX.cast<double>();
	return Q;
}

Eigen::SparseMatrix<double> AppmSolver::getElectricPermittivityOperator()
{
	const int nPe = primalMesh.getNumberOfEdges();
	Eigen::SparseMatrix<double> Meps(nPe, nPe);
	Meps.setIdentity();
	for (int i = 0; i < nPe; i++) {
		const double dualFaceArea = dualMesh.getFace(i)->getArea();
		const double primalEdgeLength = primalMesh.getEdge(i)->getLength();
		Meps.coeffRef(i, i) = dualFaceArea / primalEdgeLength;
	}
	Meps.makeCompressed();
	return Meps;
}

Eigen::SparseMatrix<double> AppmSolver::getMagneticPermeabilityOperator()
{
	const int nFaces = primalMesh.getNumberOfFaces();
	Eigen::SparseMatrix<double> Mnu(nFaces, nFaces);
	Mnu.setIdentity();
	for (int i = 0; i < nFaces; i++) {
		const double dualEdgeLength = dualMesh.getEdge(i)->getLength();
		const double primalFaceArea = primalMesh.getFace(i)->getArea();
		Mnu.coeffRef(i, i) = dualEdgeLength / primalFaceArea;
	}
	return Mnu;
}

void AppmSolver::init_multiFluid(const std::string & filename)
{
	// Read parameter file
	assert(filename.size() > 4);
	std::ifstream file(filename);
	if (!file.is_open()) {
		std::cout << "File not opened: " << filename << std::endl;
		exit(-1);
	}
	std::string line;
	while (std::getline(file, line)) {
		if (line.empty() || line.substr(0, 1) == "#") {
			continue;
		}
		std::string fluidName;
		double mass;
		int charge;
		char delimiter;
		std::stringstream ss(line);
		ss >> fluidName;
		ss >> mass >> delimiter;
		ss >> charge;

		ParticleParameters params;
		params.name = fluidName.substr(0, fluidName.size() - 1);
		params.mass = mass;
		params.electricCharge = charge;

		particleParams.push_back(params);
	}

	const int nFluids = particleParams.size();

	// Print parameters to output

	std::cout << std::endl;
	std::cout << "Fluid parameters: " << std::endl;
	std::cout << "Fluid #: Name,\tMass,\tCharge" << std::endl;
	std::cout << "=====================" << std::endl;
	for (int i = 0; i < nFluids; i++) {
		std::cout << "Fluid " << i << ": "
			<< particleParams[i].name << "\t"
			<< particleParams[i].mass << "\t"
			<< particleParams[i].electricCharge << std::endl;
	}
	std::cout << std::endl;

	// Initialize data

	const int nCells = dualMesh.getNumberOfCells();
	const int fluidStateLength = getFluidStateLength();
	fluidStates = Eigen::MatrixXd::Zero(fluidStateLength, nCells);
	fluidStates_new = Eigen::MatrixXd::Zero(fluidStateLength, nCells);
	fluidSources = Eigen::MatrixXd::Zero(fluidStateLength, nCells);
	fluidFluxes = Eigen::MatrixXd::Zero(fluidStateLength, nCells);
}

const int AppmSolver::getFluidStateLength() const {
	return 5 * this->getNFluids();
}

/* 
* Initialize fluid states (see Sod's shock tube problem).
*/
void AppmSolver::init_SodShockTube(const double zRef) {
	const int nCells = dualMesh.getNumberOfCells();
	const int nFluids = getNFluids();
	const int fluidStateLength = getFluidStateLength();
	Eigen::VectorXd leftState(fluidStateLength);
	Eigen::VectorXd rightState(fluidStateLength);

	for (int k = 0; k < nFluids; k++) {
		// TODO: edit such that number density is equal among the fluids (= quasi-neutral plasma state), not the mass density.
		const double epsilon2 = particleParams[k].mass;

		const double pL = 1;
		const double nL = 1;
		const double rhoL = epsilon2 * nL;
		const Eigen::Vector3d uL = Eigen::Vector3d::Zero();

		const double pR = 0.1;
		const double nR = 0.125;
		const double rhoR = epsilon2 * nR;
		const Eigen::Vector3d uR = Eigen::Vector3d::Zero();

		Eigen::VectorXd singleFluidStateLeft(5);
		singleFluidStateLeft.setZero();
		singleFluidStateLeft(0) = nL;
		singleFluidStateLeft(4) = 1. / epsilon2 * (pL / (Physics::gamma - 1) + 0.5 * rhoL * uL.squaredNorm());

		Eigen::VectorXd singleFluidStateRight(5);
		singleFluidStateRight.setZero();
		singleFluidStateRight(0) = nR;
		singleFluidStateRight(4) = 1. / epsilon2 * (pR / (Physics::gamma - 1) + 0.5 * rhoR * uR.squaredNorm());

		leftState.segment(5 * k, 5) = singleFluidStateLeft;
		rightState.segment(5 * k, 5) = singleFluidStateRight;
	}

	for (int i = 0; i < nCells; i++) {
		const Cell * cell = dualMesh.getCell(i);
		Eigen::VectorXd cellState(fluidStateLength);
		cellState.setConstant(std::nan("")); // NaN

		if (cell->getFluidType() == Cell::FluidType::FLUID) {
			const Eigen::Vector3d cellCenterPos = cell->getCenter();
			cellState = (cellCenterPos(2) < zRef) ? leftState : rightState;
		}
		else {
			//const double a = std::nan(""); // value of Not-A-Number
			//cellState.setConstant(a);
		}
		fluidStates.col(i) = cellState;
	}
}

void AppmSolver::init_Uniformly(const double n, const double p, const double u)
{
	const int nCells = dualMesh.getNumberOfCells();
	const int nFluids = getNFluids();
	const int fluidStateLength = getFluidStateLength();
	Eigen::VectorXd state(fluidStateLength);
	for (int k = 0; k < nFluids; k++) {
		const double epsilon2 = particleParams[k].mass;
		const double e = p / ((Physics::gamma - 1) * epsilon2 * n);
		const double etot = e + 0.5 * pow(u, 2);
		Eigen::VectorXd singleFluidState(5);
		singleFluidState(0) = n;
		singleFluidState.segment(1, 3) = n * u * Eigen::Vector3d::UnitZ();
		singleFluidState(4) = n * etot;
		state.segment(5 * k, 5) = singleFluidState;
	}
	fluidStates.setConstant(std::nan(""));
	for (int i = 0; i < nCells; i++) {
		const Cell * cell = dualMesh.getCell(i);
		if (cell->getFluidType() == Cell::FluidType::FLUID) {
			fluidStates.col(i) = state;
		}
	}
}

void AppmSolver::init_Explosion(const Eigen::Vector3d refPos, const double radius)
{
	Eigen::VectorXd state_lo = Physics::primitive2state(1, 1, Eigen::Vector3d::Zero());
	Eigen::VectorXd state_hi = Physics::primitive2state(1, 5, Eigen::Vector3d::Zero());

	const int nCells = dualMesh.getNumberOfCells();
	const int nFluids = getNFluids();
	for (int cIdx = 0; cIdx < nCells; cIdx++) {
		const Cell * cell = dualMesh.getCell(cIdx);
		if (cell->getFluidType() != Cell::FluidType::FLUID) { continue; } // Skip cells that are not Fluid type
		const Eigen::Vector3d cc = cell->getCenter();
		for (int fluidIdx = 0; fluidIdx < nFluids; fluidIdx++) {
			const bool isInside = (cc - refPos).norm() < radius;
			const Eigen::VectorXd state = isInside ? state_hi : state_lo;
			fluidStates.col(cIdx).segment(5 * fluidIdx, 5) = state;
		}
	}
}

/**
* For Testing purposes. 
* Set E_h uniform in a given direction. 
*/

void AppmSolver::set_Efield_uniform(const Eigen::Vector3d direction)
{
	E_h.setZero();
	assert(E_h.size() == primalMesh.getNumberOfEdges());
	for (int i = 0; i < E_h.size(); i++) {
		const Edge * edge = primalMesh.getEdge(i);
		const Eigen::Vector3d edgeDir = edge->getDirection().normalized();
		const double edgeLength = edge->getLength();
		E_h(i) = edgeDir.dot(direction) * edgeLength;
	}
}

/** 
* For testing purposes. 
* Set B-field vectors azimuthally with respect to z-unit vector. 
*/
void AppmSolver::set_Bfield_azimuthal()
{
	const int nPrimalFaces = primalMesh.getNumberOfFaces();
	const Eigen::Vector3d zvec = Eigen::Vector3d::UnitZ();	// z-Unit vector

	for (int idx = 0; idx < nPrimalFaces; idx++) {
		const Face * face = primalMesh.getFace(idx);
		const Eigen::Vector3d fc = face->getCenter();
		const Eigen::Vector3d fn = face->getNormal();
		const double fA = face->getArea();

		Eigen::Vector3d rvec = fc;
		rvec(2) = 0; 		// Radial position

		Eigen::Vector3d thetaVec;
		thetaVec = zvec.cross(rvec);       // Azimuthal vector
		B_h(idx) = thetaVec.dot(fn) * fA;  // Projection of azimuthal vector onto face normal
	}
}

/** 
* Initialize matrix for interpolating electric field from primal edges to dual cell centers.*
*/
void AppmSolver::initPerotInterpolationMatrix()
{
	const std::string text = "Initialize Perot Interpolation Matrix";
	std::cout << text << std::endl;

	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;

	const int nCells = dualMesh.getNumberOfCells();
	const int nEdges = primalMesh.getNumberOfEdges();
	assert(nCells > 0);
	assert(nEdges > 0);
	assert(nEdges == E_h.rows());

	for (int cIdx = 0; cIdx < nCells; cIdx++) {
		const Cell * cell = dualMesh.getCell(cIdx);
		const std::vector<Face*> cellFaces = cell->getFaceList();
		for (auto face : cellFaces) {
			const int faceIdx = face->getIndex();
			const int edgeIdx = dualMesh.getAssociatedPrimalEdgeIndex(faceIdx);
			const Edge * edge = primalMesh.getEdge(edgeIdx);

			const Eigen::Vector3d fc = face->getCenter();
			const Eigen::Vector3d cc = cell->getCenter();
			const double dV = cell->getVolume();
			const double dA = face->getArea();
			const double dL = edge->getLength();

			const Eigen::Vector3d r = cc - fc;
			const Eigen::Vector3d L = edge->getDirection().normalized();
			const Eigen::Vector3d n = face->getNormal();
			const int incidence = r.dot(n) > 0 ? 1 : -1;

			const double value = 1. / dL * L.dot(n) * incidence * dA / dV;

			const int col = edgeIdx;
			for (int i = 0; i < 3; i++) {
				int row = 3 * cIdx + i;
				triplets.push_back(T(row, col, value * r(i)));
			}
		}
	}
	M_perot = Eigen::SparseMatrix<double>(3*nCells, nEdges);
	M_perot.setFromTriplets(triplets.begin(), triplets.end());
	M_perot.makeCompressed();
	
	std::cout << text << ": DONE" << std::endl;
	// Eigen::sparseMatrixToFile(M_perot, "Mperot.dat");
}

const Eigen::Matrix3Xd AppmSolver::getEfieldAtCellCenter()
{
	const int nCells = dualMesh.getNumberOfCells();
	assert(M_perot.size() > 0);
	assert(M_perot.nonZeros() > 0);

	Eigen::VectorXd E_cc = M_perot * E_h;

	//// Eigen library uses Maps to perform reshaping operations; more precisely, it is another 'view' on the underlying data. 
	//Eigen::Map<Eigen::Matrix3Xd> E_reshaped(E_cc.data(), 3, nCells);
	//Eigen::MatrixXd E_cellCenter(3, nCells);
	//E_cellCenter.setZero();
	//assert(E_cc.rows() == 3 * nCells);
	//for (int i = 0; i < nCells; i++) {
	//	E_cellCenter.col(i) = E_cc.segment(3*i, 3);
	//	assert((E_cellCenter.col(i).array() == E_reshaped.col(i).array()).all());
	//}
	//E_cellCenter = E_reshaped;

	// The class Map allows to 'view' the underlying data in different manner (e.g., doing a reshape operation)
	Eigen::Map<Eigen::Matrix3Xd> result(E_cc.data(), 3, nCells);
	return result;
}

//void AppmSolver::initMsigma()
//{
//	init_Msigma_solid();
//}

/**
* Ohm's law for a solid body (j = Msigma * e).
*/
Eigen::SparseMatrix<double> AppmSolver::init_Msigma_solid()
{
	assert(false);
	// outdated; do not use this function anymore!

	const int nEdges = primalMesh.getNumberOfEdges();
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;

	for (int i = 0; i < nEdges; i++) {
		const Face * face = dualMesh.getFace(i);
		const Eigen::Vector3d fc = face->getCenter();
		const double faceArea = face->getArea();
		const double edgeLength = primalMesh.getEdge(i)->getLength();

		//const bool isInside = fc.segment(0, 2).norm() < 0.35;
		//const double sigma = isInside ? 10 : 1;
		const double sigma = 0;

		const double value = sigma * faceArea / edgeLength;
		triplets.push_back(T(i, i, value));
	}
	Eigen::SparseMatrix<double> Msigma = Eigen::SparseMatrix<double>(nEdges, nEdges);
	Msigma.setFromTriplets(triplets.begin(), triplets.end());
	Msigma.makeCompressed();
	return Msigma;
}

Eigen::SparseMatrix<double> AppmSolver::init_Msigma_fluid()
{
	assert(false);
	// outdated; Do not use this function anymore!

	// Setup of equation J = M * E + Jb
	assert(J_h.size() > E_h.size());

	// TODO extend to multi-fluid model
	std::cout << "Number of fluids: " << getNFluids() << std::endl;
	assert(getNFluids() == 1);
	const int fluidIdx = 0;
	const double q = 1;

	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;

	// For each dual face
	for (int dualFaceIdx = 0; dualFaceIdx < J_h.size(); dualFaceIdx++) {
		const Face * dualFace = dualMesh.getFace(dualFaceIdx);
		const Eigen::Vector3d dualFaceNormal = dualFace->getNormal();

		// The source term for Lorentz force is defined only in fluid cells
		if (!dualFace->hasFluidCells()) {
			continue;
		}

		// get list of adjacient cells
		const std::vector<Cell*> adjacientCells = dualFace->getCellList();
		assert(adjacientCells.size() >= 1);
		assert(adjacientCells.size() <= 2);
		assert((adjacientCells.size() == 1 && dualFace->isBoundary())
			|| (adjacientCells.size() == 2 && !dualFace->isBoundary()));

		// For each adjacient cell
		for (auto cell : adjacientCells) {
			// Lorentz force is defined only for fluid cells
			if (cell->getFluidType() != Cell::FluidType::FLUID) {
				continue;
			}
			const double cellVolume = cell->getVolume();
			assert(cellVolume > 0);

			const std::vector<Face*> cellFaces = cell->getFaceList();
			const Eigen::Vector3d cc = cell->getCenter();
			const double n = fluidStates(5 * fluidIdx, cell->getIndex()); // number density in cell 

			// For each face of that cell
			for (auto face : cellFaces) {
				const int faceIdx = face->getIndex();

				// Get corresponding primal edge
				const int edgeIdx = dualMesh.getAssociatedPrimalEdgeIndex(faceIdx);
				assert(edgeIdx >= 0);
				assert(edgeIdx < primalMesh.getNumberOfEdges());
				const Edge * edge = primalMesh.getEdge(edgeIdx);
				//const Eigen::Vector3d edgeDirection = edge->getDirection().normalized();

				const Eigen::Vector3d nj = face->getNormal();
				const Eigen::Vector3d L = edge->getDirection().normalized();
				//assert(abs(abs(L.dot(nj)) - 1) < 0.001); // check if edge direction and face normal are (almost) parallel

				const Eigen::Vector3d fc = face->getCenter();
				const Eigen::Vector3d r = fc - cc; // Position vector of face center relative to cell center
				assert(r.norm() > 0);

				if (isElectricLorentzForceActive) {

					// Use Perot's interpolation method to obtain a cell-centered value 
					// of the electric field that is defined at dual faces (= primal edges)
					const double r_dot_dualFaceNormal = r.dot(dualFaceNormal);
					const int incidence = r.dot(nj) > 0 ? 1 : -1;
					const double dL = edge->getLength();
					const double dA = face->getArea();
					const double dV = cell->getVolume();
					double value = n * r_dot_dualFaceNormal * 1. / dL * (L.dot(nj)) * incidence * dA / dV;

					// factor due to definition of fluid flux: 
					// - at interior faces: Rusanov scheme f = 0.5 * ( (nu)_k + (nu)_k+1 ) - 0.5 * (n_k+1 - n_k)
					// - at fluid boundary faces: 
					//      f = (nu)_k (opening condition)				
					//      f = ...    (wall condition) (TODO)				
					value *= (dualFace->getFluidType() == Face::FluidType::OPENING) ? 1 : 0.5;

					// factor due to consistency of electric current and fluid flux: 
					// J_h = j * A = q * f * A, where: 
					//   J_h: face-integrated current
					//   j:   current density at face
					//   q:   ionization degree of fluid
					//   f:   fluid flux at face
					//   A:   face area
					value *= q * dualFace->getArea();

					// Store coefficient in sparse matrix
					const int row = dualFaceIdx;
					const int col = edgeIdx;
					triplets.push_back(T(row, col, value));
				} // if isElectricLorentzForceActive
			} // for cellFaces
		} // for adjacientCells
	} // for dualFaceIdx
	Eigen::SparseMatrix<double> Msigma;
	Msigma = Eigen::SparseMatrix<double>(J_h.size(), E_h.size());
	Msigma.setFromTriplets(triplets.begin(), triplets.end());
	Msigma.makeCompressed();
	return Msigma;
}

/** 
* Setup of equation J_h = Msigma * E_h + J_h_aux, which results of the consistent-implicit formulation.
*/
void AppmSolver::get_Msigma_consistent(const double dt, Eigen::SparseMatrix<double>& Msigma, Eigen::VectorXd & jaux)
{
	const int nPrimalEdges = primalMesh.getNumberOfEdges();
	const int nDualFaces = dualMesh.getNumberOfFaces();
	const int nFluids = getNFluids();

	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;

	std::cout << "isElectricLorentzForceActive: " << isElectricLorentzForceActive << std::endl;
	std::cout << "isMagneticLorentzForceActive: " << isMagneticLorentzForceActive << std::endl;
	std::cout << "isMomentumFluxActive: " << isMomentumFluxActive << std::endl;
	assert(isMomentumFluxActive);

	std::cout << "Warning: check definition of implicit-consistent current: factors due to Rusanov scheme in magnetic Lorentz force" << std::endl;
	std::cout << "Warning: check definition of implicit-consistent current: factors due to Rusanov scheme at domain boundary" << std::endl;

	// Set parameters to zero
	Msigma = Eigen::SparseMatrix<double>(nDualFaces, nPrimalEdges);
	jaux = Eigen::VectorXd::Zero(nDualFaces);

	// For each fluid ...
	for (int fluidIdx = 0; fluidIdx < nFluids; fluidIdx++) {
		// species ionization degree 
		const int q = particleParams[fluidIdx].electricCharge;
		if (q == 0) { // skip species that are neutral
			continue; 
		}

		// For each face in dual mesh at which we evaulate implicit-consistent current value
		for (int jdx = 0; jdx < nDualFaces; jdx++) {
			const Face * dualFace = dualMesh.getFace(jdx);
			if (!dualFace->hasFluidCells()) { // skip faces that are not adjacient to a fluid cell
				continue;
			}

			const Eigen::Vector3d fn = dualFace->getNormal();
			const double fA = dualFace->getArea();
			const std::vector<Cell*> adjacientCells = dualFace->getCellList();
			assert(adjacientCells.size() >= 1);
			assert(adjacientCells.size() <= 2);
			assert((adjacientCells.size() == 1 &&  dualFace->isBoundary())
				|| (adjacientCells.size() == 2 && !dualFace->isBoundary()));

			// For each adjacient cell
			for (auto cell : adjacientCells) {
				// model is defined only for fluid cells; skip other cells
				if (cell->getFluidType() != Cell::FluidType::FLUID) {
					continue;
				}

				// Number density in dual cell
				const double numberDensity = fluidStates(5 * fluidIdx, cell->getIndex());

				// List of faces that make this cell
				std::vector<Face*> cellFaces = cell->getFaceList();
				for (auto face : cellFaces) {
					// Index of primal edge that is associated with dual face
					const int edx = dualMesh.getAssociatedPrimalEdgeIndex(face->getIndex());
					assert(edx >= 0 && edx < primalMesh.getNumberOfEdges());
					const Edge * edge = primalMesh.getEdge(edx);
					
					const Eigen::Vector3d fc = face->getCenter();
					const Eigen::Vector3d cc = cell->getCenter();
					
					// Use Perot's interpolation method to obtain cell-centered value
					// of the electric field that is defined at dual faces (= primal edges)

					const Eigen::Vector3d rvec = fc - cc; // Position vector of face center relative to cell center
					const Eigen::Vector3d nj = face->getNormal();
					const int incidence = rvec.dot(nj) > 0 ? 1 : -1;
					const Eigen::Vector3d Lhat = edge->getDirection().normalized();
					const double rvec_proj_fn = rvec.dot(fn);
					const double dL = edge->getLength();
					const double dA = face->getArea();
					const double dV = cell->getVolume();

					
					// Consistent-implicit electric current due to electric Lorentz force (F = q*n*E)
					if (isElectricLorentzForceActive) {

						// Electric field projected in direction of dual face and interpolated from dual face centers
						double value = q * numberDensity * rvec_proj_fn * Lhat.dot(nj) * incidence * dA / dV * (1. / dL);
						value *= dt * q * fA;

						// Factor due to definition of reconstructed fluid flux:
						// - At inner faces (Rusanov scheme): f = 0.5 * ( (nu)_k + (nu)_{k+1} ) - 0.5 * s * (n_{k+1} - n_k)
						// - At opening faces: f = (nu)_k
						// - At boundary faces: ???                  <<<---- TODO
						const bool isFaceTypeOpening = dualFace->getFluidType() == Face::FluidType::OPENING;
						value *= isFaceTypeOpening ? 1 : 0.5;

						triplets.push_back(T(jdx, edx, value));
					} // end if isElectricLorentzForceActive

						
					// Consistent-implicit electric current due to momentum flux 
					if (isMomentumFluxActive) {
						const double momentumFlux = getRusanovFluxExplicit(face->getIndex(), fluidIdx)(1);
						jaux(jdx) += -dt / dV * dA * nj.dot(fn) * momentumFlux * q * fA;
					}


				} // end for cellFaces			

				// Consistent-implicit electric current due to magnetic Lorentz force (F = q * (nu) x B)
				if (isMagneticLorentzForceActive) {
					const Eigen::Vector3d B = B_vertex.col(cell->getIndex());
					const Eigen::Vector3d nu = fluidStates.col(cell->getIndex()).segment(5 * fluidIdx + 1, 3);
					double value = q * nu.cross(B).dot(fn);
					value *= dt * q * fA;
					value *= 0.5; // this factor is due to defintion of reconstructed flux (Rusanov scheme)
					jaux(jdx) += value;
				}

			} // end for adjacient cells
		} // end for dual faces
	} // end for fluidIdx

	Msigma.setFromTriplets(triplets.begin(), triplets.end());
	Msigma.makeCompressed();

	assert(Msigma.nonZeros() > 0);
}

const double AppmSolver::terminalVoltageBC_sideA(const double time, const double t0, const double t1, const double tscale) const
{
	double f1 = 0.5 * (1 + tanh((time - t0) / tscale)) ; 
	double f2 = 0.5 * (1 + tanh(- (time - t1) / tscale));
	return f1 * f2;
}

const double AppmSolver::terminalVoltageBC_sideB(const double time) const
{
	return 0.0;
}

const double AppmSolver::currentDensityBC(const double time) const
{
	const double t0 = 2;
	const double t1 = 5;
	const double tscale = 0.5;
	const double currentValue = std::exp(-pow((time - t0) / tscale, 2)) - std::exp(-pow((time - t1) / tscale, 2));
	return currentValue;
}

/** 
* Get timestep size for next iteration.
*/
const double AppmSolver::getNextFluidTimestepSize() const
{
	const int nFaces = dualMesh.getNumberOfFaces();
	const int nFluids = this->getNFluids();

	Eigen::VectorXd dt_faces(nFaces);
	dt_faces.setConstant(std::numeric_limits<double>::max());

	/*
	* For each face that has an adjacient fluid cell,
	* determine timestep size such that a travelling wave does not
	* exceed the distance from cell center to face.
	*/
	for (int fidx = 0; fidx < nFaces; fidx++) {
		const Face * face = dualMesh.getFace(fidx); 
		const Eigen::Vector3d fc = face->getCenter();
		const Eigen::Vector3d faceNormal = face->getNormal();
		if (!face->hasFluidCells()) {
			continue;
		}
		// Get list of cells (one or two) that are adjacient to this face.
		const std::vector<Cell*> faceCells = face->getCellList();

		Eigen::VectorXd dt_local(faceCells.size());
		dt_local.setConstant(std::numeric_limits<double>::max());
		bool is_timestepDefined = false;

		for (int i = 0; i < faceCells.size(); i++) {
			const Cell * cell = faceCells[i];
			if (cell->getFluidType() != Cell::FluidType::FLUID) {
				continue;
			}
			const int cellIdx = cell->getIndex();
			const Eigen::Vector3d cc = cell->getCenter();
			const double dx = std::abs((fc - cc).dot(faceNormal));
			assert(dx > 1e-6);
			const Eigen::VectorXd cellState = fluidStates.col(cell->getIndex());
			Eigen::VectorXd wavespeeds(nFluids);
			for (int fluidIdx = 0; fluidIdx < nFluids; fluidIdx++) {
				Eigen::Vector3d q = getFluidState(cellIdx, fluidIdx, faceNormal);
				const double s = Physics::getMaxWavespeed(q);
				wavespeeds(fluidIdx) = s;
			}
			assert((wavespeeds.array() > 0).all());
			const double smax = wavespeeds.maxCoeff();
			assert(smax > 0);
			assert(smax > 1e-6);
			assert(smax < 1e6);
			dt_local(i) = dx / smax;
		}
		assert(dt_local.allFinite());
		dt_faces(fidx) = dt_local.minCoeff();
	}
	const double dt = dt_faces.minCoeff();
	assert(dt > 1e-12);
	return dt;
}

//const Eigen::VectorXd AppmSolver::getFluidState(const int cellIdx, const int fluidIdx) const 
//{
//	assert(cellIdx >= 0);
//	assert(cellIdx < fluidStates.cols());
//	assert(fluidIdx >= 0);
//	assert(fluidIdx < getNFluids());
//	return fluidStates.col(cellIdx).segment(5*fluidIdx, 5);
//}

const Eigen::Vector3d AppmSolver::getFluidState(const int cellIdx, const int fluidIdx, const Eigen::Vector3d & faceNormal) const
{
	assert(cellIdx >= 0);
	assert(cellIdx < fluidStates.cols());
	assert(fluidIdx >= 0);
	assert(fluidIdx < this->getNFluids());
	const Eigen::VectorXd state = fluidStates.col(cellIdx).segment(5 * fluidIdx, 5);
	return Eigen::Vector3d(state(0), state.segment(1,3).dot(faceNormal), state(4));
}

const int AppmSolver::getOrientation(const Cell * cell, const Face * face) const
{
	assert(cell != nullptr);
	assert(face != nullptr);
	assert(face->isAdjacient(cell));
	const Eigen::Vector3d fc = face->getCenter();
	const Eigen::Vector3d fn = face->getNormal();
	const Eigen::Vector3d cc = cell->getCenter();
	int orientation = (fc - cc).dot(fn) > 0 ? 1 : -1;
	return orientation;
}

//const Eigen::Vector3d AppmSolver::getFluidFluxFromState(const Eigen::Vector3d & q) const
//{
//	assert(q.norm() > 0);
//	Eigen::Vector3d flux;
//	flux(0) = q(1);
//	flux(1) = 0.5 * (3 - Physics::gamma) * pow(q(1), 2) / q(0) + (Physics::gamma - 1) * q(2);
//	flux(2) = Physics::gamma * q(1) * q(2) / q(0) - 0.5 * (Physics::gamma - 1) * pow(q(1) / q(0), 2) * q(1);
//	return flux;
//}

/**
* Get the face flux with explicit Rusanov scheme. 
*/
const Eigen::Vector3d AppmSolver::getRusanovFluxExplicit(const int faceIdx, const int fluidIdx) const
{
	Eigen::Vector3d qL, qR;
	getAdjacientCellStates(faceIdx, fluidIdx, qL, qR);
	const double showOutput = faceIdx == faceIdxRef;
	return Physics::getRusanovFlux(qL, qR, showOutput); 
}

const Eigen::Vector3d AppmSolver::getRusanovFluxImEx(const int faceIdx, const int fluidIdx, const double dt) 
{
	const Face * face = dualMesh.getFace(faceIdx);
	const Eigen::Vector3d faceNormal = face->getNormal();
	const Face::FluidType faceFluidType = face->getFluidType();

	Eigen::Vector3d qL, qR;
	const std::pair<int, int> adjacientCellIdx = getAdjacientCellStates(faceIdx, fluidIdx, qL, qR);

	double implicitExtraTerms = 0;
	int cellIndex;

	if (faceFluidType != Face::FluidType::INTERIOR) {
		assert(adjacientCellIdx.first == -1 || adjacientCellIdx.second == -1);
	}
	const bool isMassFluxImexScheme = this->massFluxScheme == MassFluxScheme::IMPLICIT_EXPLICIT;
	const bool isInteriorFace = faceFluidType == Face::FluidType::INTERIOR;
	const bool isBoundaryFace = !isInteriorFace;
	if (isMassFluxImexScheme) {
		cellIndex = adjacientCellIdx.first;
		double extra_L = 0;
		if (cellIndex >= 0) {
			extra_L = getImplicitExtraTermMomentumFlux(cellIndex, faceNormal, fluidIdx);
		}
		cellIndex = adjacientCellIdx.second;
		double extra_R = 0;
		if (cellIndex >= 0) {
			extra_R = getImplicitExtraTermMomentumFlux(cellIndex, faceNormal, fluidIdx);
		}
		implicitExtraTerms -= dt * extra_L;
		implicitExtraTerms -= dt * extra_R;
		faceFluxesImExRusanov.coeffRef(faceIdx, fluidIdx) = implicitExtraTerms;
	}
	

	Eigen::Vector3d faceFlux = getRusanovFluxExplicit(faceIdx, fluidIdx);
	faceFlux(0) += implicitExtraTerms;
	return faceFlux;
}

const double AppmSolver::getImplicitExtraTermMomentumFlux(const int cellIdx, const Eigen::Vector3d & faceNormal, const int fluidIdx) const
{
	assert(cellIdx >= 0); 
	assert(std::fabs(faceNormal.norm() - 1) <= 4*std::numeric_limits<double>::epsilon());

	const Cell * cell = dualMesh.getCell(cellIdx);
	const double cellVolume = cell->getVolume();
	const std::vector<Face*> cellFaces = cell->getFaceList();
	double sumFaceFluxes = 0;
	for (auto face : cellFaces) {
		const int faceIdx = face->getIndex();
		const Eigen::Vector3d faceFlux = getRusanovFluxExplicit(faceIdx, fluidIdx);
		const double momentumFlux = faceFlux(1);
		const double faceArea = face->getArea();
		const Eigen::Vector3d fn = face->getNormal();
		const double orientation = fn.dot(faceNormal) * getOrientation(cell, face);
		sumFaceFluxes += momentumFlux * faceArea * orientation;
	}
	const double result = 0.5 * 1./cellVolume * sumFaceFluxes;
	return result;
}

const std::pair<int,int> AppmSolver::getAdjacientCellStates(const int faceIdx, const int fluidIdx, Eigen::Vector3d & qL, Eigen::Vector3d & qR) const
{
	const Face * face = dualMesh.getFace(faceIdx);
	const std::vector<Cell*> faceCells = face->getCellList();
	assert(faceCells.size() >= 1);

	const Eigen::Vector3d faceNormal = face->getNormal();
	const int orientation = getOrientation(faceCells[0], face);
	qL.setZero();
	qR.setZero();
	int idxL = -1;
	int idxR = -1;

	const Face::FluidType faceFluidType = face->getFluidType();
	switch (faceFluidType) {
	case Face::FluidType::INTERIOR:
		assert(faceCells.size() == 2);
		assert(faceCells[0]->getFluidType() == Cell::FluidType::FLUID);
		assert(faceCells[1]->getFluidType() == Cell::FluidType::FLUID);
		if (orientation > 0) {
			idxL = faceCells[0]->getIndex();
			idxR = faceCells[1]->getIndex();
		}
		else {
			idxL = faceCells[1]->getIndex();
			idxR = faceCells[0]->getIndex();
		}
		qL = getFluidState(idxL, fluidIdx, faceNormal);
		qR = getFluidState(idxR, fluidIdx, faceNormal);
		break;

	case Face::FluidType::OPENING:
		assert(faceCells.size() == 1);
		if (orientation > 0) {
			idxL = faceCells[0]->getIndex();
			qL = getFluidState(idxL, fluidIdx, faceNormal);
			qR = qL;
		}
		else {
			idxR = faceCells[0]->getIndex();
			qR = getFluidState(idxR, fluidIdx, faceNormal);
			qL = qR;
		}
		break;

	case Face::FluidType::WALL:
		assert(faceCells.size() >= 1 && faceCells.size() <= 2);
		if (faceCells[0]->getFluidType() == Cell::FluidType::FLUID) {
			const int idx = faceCells[0]->getIndex();

			if (orientation > 0) {
				idxL = faceCells[0]->getIndex();
				qL = getFluidState(idxL, fluidIdx, faceNormal);
				qR = qL.cwiseProduct(Eigen::Vector3d(1, -1, 1));
			}
			else {
				idxR = faceCells[0]->getIndex();
				qR = getFluidState(idxR, fluidIdx, faceNormal);
				qL = qR.cwiseProduct(Eigen::Vector3d(1, -1, 1));
			}
		}
		else {
			assert(faceCells.size() >= 2);
			assert(faceCells[1]->getFluidType() == Cell::FluidType::FLUID);
			assert(false); // Not yet implemented
		}
		break;

	default:
		std::cout << "Face Fluid Type not implemented" << std::endl;
		exit(-1);
	}

	return std::pair<int, int>(idxL, idxR);
}

const double AppmSolver::getMomentumUpdate(const int k, const Eigen::Vector3d & nvec, const int fluidIdx) const
{
	double result = 0;
	const Cell * cell = dualMesh.getCell(k);
	const std::vector<Face*> cellFaces = cell->getFaceList();
	for (auto face : cellFaces) {
		const double faceArea = face->getArea();
		const Eigen::Vector3d faceNormal = face->getNormal();
		const Eigen::Vector3d faceFlux = getRusanovFluxExplicit(face->getIndex(), fluidIdx);
		const double momentumFlux = faceFlux(1);
		double localResult = momentumFlux * faceArea * faceNormal.dot(nvec);
		result += localResult;
	}
	result /= cell->getVolume();
	return result;
}


void AppmSolver::interpolateMagneticFluxToPrimalVertices()
{
	B_vertex.setZero();
	// For each cell ...
	//const int nCells = primalMesh.getNumberOfCells();
	// For each cell that has a Piola map
	const int nCells = rt_piolaMatrix.size();
	Eigen::VectorXi countVertexVisits = Eigen::VectorXi::Zero(primalMesh.getNumberOfVertices());

	//const Eigen::VectorXd B_h = maxwellSolver->getBstate();
	const int nPrimalFaces = primalMesh.getNumberOfFaces();
	//Eigen::VectorXd B_h(nPrimalFaces); 
	//B_h.setZero();

	for (int cidx = 0; cidx < nCells; cidx++) {
		const Cell * cell = primalMesh.getCell(cidx);
		const std::vector<Face*> cellFaces = cell->getFaceList();
		const std::vector<Vertex*> bottomVertices = cellFaces[3]->getVertexList();
		const std::vector<Vertex*> topVertices = cellFaces[4]->getVertexList();

		// Piola map
		const Eigen::Matrix3d & BK = rt_piolaMatrix[cidx];
		const Eigen::Vector3d & bK = rt_piolaVector[cidx];
		const double detBK = BK.determinant();

		// Get prism vertices
		std::vector<Vertex*> cellVertices(6);
		for (int i = 0; i < 3; i++) {
			cellVertices[i] = bottomVertices[i];
			cellVertices[i + 3] = topVertices[i];
		}

		//Eigen::Matrix3cd vertexCoords(3, 6);
		Eigen::MatrixXd refCoords(3, 6);
		Eigen::VectorXi vertexIdx(6);
		for (int i = 0; i < 6; i++) {
			const Vertex * v = cellVertices[i];
			const Eigen::Vector3d pos = v->getPosition();
			//vertexCoords.col(i) = pos;
			vertexIdx(i) = v->getIndex();
			refCoords.col(i) = BK.inverse() * (pos - bK);
		}
		const double tolerance = 16 * std::numeric_limits<double>::epsilon();
		// Check that reference coordinates are in [0,1]
		assert((refCoords.array() >= -tolerance).all());
		assert(((refCoords.array() - 1) <= tolerance).all());

		// Get coefficients for RT interpolation
		Eigen::VectorXd coeff(5);
		for (int idx = 0; idx < 5; idx++) {
			const Face * face = cellFaces[idx];
			const int faceIncidence = cell->getOrientation(face);
			const int fidx = face->getIndex();
			const int factor = 1;// face->isBoundary() ? 2 : 1;
			coeff(idx) = faceIncidence * B_h(fidx) * factor;
		}

		// Interpolate B-field with RT-basis functions in actual space
		Eigen::Matrix3Xd B_vertex_local = Eigen::Matrix3Xd::Zero(3, 6);
		for (int idx = 0; idx < 5; idx++) {
			for (int i = 0; i < 6; i++) {
				const int vIdx = vertexIdx(i);
				const Eigen::Vector3d refPos = refCoords.col(i);
				const Eigen::Vector3d rt = Numerics::raviartThomasBasis(idx, refPos);
				const int directionFactor = (idx >= 3) ? 2 : 1; // Because: prism in z-direction (i.e., face idx < 3) are quads with normal vector perpendicular to z-axis, and face idx >= 3 are triangles. Hence, det(B) acts differently.
				B_vertex_local.col(i) += directionFactor * coeff(idx) / detBK * BK * rt;
			}
		}
		for (int i = 0; i < 6; i++) {
			const int vIdx = vertexIdx(i);
			countVertexVisits(vIdx) += 1;
			B_vertex.col(vIdx) += B_vertex_local.col(i);
		}
	}
	for (int i = 0; i < primalMesh.getNumberOfVertices(); i++) {
		if (countVertexVisits(i) > 0) {
			B_vertex.col(i) /= countVertexVisits(i);
		}
	}
}


const Eigen::Matrix3Xd AppmSolver::getPrismReferenceCoords(const int nSamples)
{
	Eigen::Matrix3Xd refCoords(3, (nSamples + 2) * (nSamples + 1) / 2);
	int sIdx = 0;
	for (int i = 0; i <= nSamples; i++) {
		if (i == 0) {
			refCoords.col(sIdx++) = Eigen::Vector3d(0, 0, 0);
			continue;
		}
		Eigen::Vector3d a, b;
		a = 1.0*i / nSamples * Eigen::Vector3d(1, 0, 0);
		b = 1.0*i / nSamples * Eigen::Vector3d(0, 1, 0);
		Eigen::Vector3d d = b - a;
		const Eigen::VectorXd samples = Eigen::VectorXd::LinSpaced(i + 1, 0., 1.);
		for (int k = 0; k < samples.size(); k++) {
			const Eigen::Vector3d v = a + samples(k) * (b - a);
			refCoords.col(sIdx++) = v;
		}
	}
	assert(sIdx == refCoords.cols());
	const int n = refCoords.cols();
	Eigen::Matrix3Xd refCoords3d = refCoords.replicate(1, nSamples + 1);
	for (int i = 0; i <= nSamples; i++) {
		const int startCol = i * n;
		refCoords3d.block(2, startCol, 1, n).array() = 1.0 * i / nSamples;
	}
	return refCoords3d;
}

void AppmSolver::init_meshes(const PrimalMesh::PrimalMeshParams & primalParams)
{
	std::cout << "Init primal mesh" << std::endl;

	primalMesh = PrimalMesh(primalParams);
	primalMesh.init();
	primalMesh.writeToFile();
	primalMesh.writeXdmf();
	primalMesh.check();

	if (primalMesh.getNumberOfCells() == 0) {
		std::cout << "Primal mesh has no cells" << std::endl;
		return;
	}

	std::cout << "Init dual mesh" << std::endl;
	dualMesh = DualMesh();
	dualMesh.init_dualMesh(primalMesh, primalParams.getElectrodeRadius());
	dualMesh.writeToFile();
	dualMesh.writeXdmf();

	std::cout << "Dual mesh has " << dualMesh.getNumberOfVertices() << " vertices" << std::endl;
}

void AppmSolver::writeXdmf()
{
	const int nTimesteps = timeStamps.size();
	
	const std::string filename = "appm.xdmf";
	std::string gridPrimalEdges;
	std::string gridPrimalFaces;
	std::string gridDualEdges;
	std::string gridDualFaces;

	std::ofstream file(filename);
	file << "<?xml version = \"1.0\" ?>" << std::endl;
	file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
	file << "<Xdmf Version=\"3.0\" xmlns:xi=\"[http://www.w3.org/2001/XInclude]\">" << std::endl;
	file << "<Domain>" << std::endl;
	file << "<Grid Name=\"Time Grid\" GridType=\"Collection\" CollectionType=\"Temporal\">" << std::endl;
	for (int i = 0; i < nTimesteps; i++) {
		const double time = this->timeStamps[i];

		file << "<Grid Name=\"Grid of Grids\" GridType=\"Tree\">" << std::endl;
		file << "<Time Value=\"" << time << "\" />" << std::endl;
		file << xdmf_GridPrimalEdges(i) << std::endl;
		file << xdmf_GridPrimalFaces(i) << std::endl;
		file << xdmf_GridDualEdges(i) << std::endl;
		file << xdmf_GridDualFaces(i) << std::endl;
		file << "</Grid>" << std::endl;
	}
	file << "</Grid>" << std::endl;
	file << "</Domain>" << std::endl;
	file << "</Xdmf>" << std::endl;
}

void AppmSolver::writeXdmfDualVolume()
{
	const int nTimesteps = timeStamps.size();
	const std::string filename = "appm-volume.xdmf";
	std::ofstream file(filename);
	file << "<?xml version = \"1.0\" ?>" << std::endl;
	file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
	file << "<Xdmf Version=\"3.0\" xmlns:xi=\"[http://www.w3.org/2001/XInclude]\">" << std::endl;
	file << "<Domain>" << std::endl;
	file << "<Grid Name=\"Time Grid\" GridType=\"Collection\" CollectionType=\"Temporal\">" << std::endl;
	for (int i = 0; i < nTimesteps; i++) {
		const double time = this->timeStamps[i];
		file << "<Time Value=\"" << time << "\" />" << std::endl;
		file << xdmf_GridDualCells(i) << std::endl;
	}
	file << "</Grid>" << std::endl;
	file << "</Domain>" << std::endl;
	file << "</Xdmf>" << std::endl;
}

//void AppmSolver::writeXdmfDualFaceFluxes()
//{
//	const int nTimesteps = timeStamps.size();
//	const std::string filename = "appm-faceFluxes.xdmf";
//	std::ofstream file(filename);
//	file << "<?xml version = \"1.0\" ?>" << std::endl;
//	file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
//	file << "<Xdmf Version=\"3.0\" xmlns:xi=\"[http://www.w3.org/2001/XInclude]\">" << std::endl;
//	file << "<Domain>" << std::endl;
//	file << "<Grid Name=\"Time Grid\" GridType=\"Collection\" CollectionType=\"Temporal\">" << std::endl;
//	for (int i = 0; i < nTimesteps; i++) {
//		const double time = this->timeStamps[i];
//		file << "<Time Value=\"" << time << "\" />" << std::endl;
//		file << xdmf_GridDualFaces(i) << std::endl;
//	}
//	file << "</Grid>" << std::endl;
//	file << "</Domain>" << std::endl;
//	file << "</Xdmf>" << std::endl;
//}


void AppmSolver::writeOutput(const int iteration, const double time)
{
	timeStamps.push_back(time);
	std::cout << "Write output at iteration " << iteration << ", time = " << time << std::endl;
	const std::string filename = (std::stringstream() << "appm-" << iteration << ".h5").str();

	const int nPrimalFaces = primalMesh.getNumberOfFaces();
	const int nDualEdges = dualMesh.getNumberOfEdges();
	const int nDualFaces = dualMesh.getNumberOfFaces();
	const int nDualCells = dualMesh.getNumberOfCells();

	H5Writer h5writer(filename);

	// Fluid states
	//fluidSolver->writeStates(h5writer);
	writeFluidStates(h5writer);

	// Maxwell states
	//maxwellSolver->writeStates(h5writer);
	writeMaxwellStates(h5writer);

	// Interpolated values of B-field to primal vertices
	h5writer.writeData(B_vertex, "/Bvertex");

	Eigen::VectorXd timeVec(1);
	timeVec(0) = time;
	h5writer.writeData(timeVec, "/time");
	Eigen::VectorXi iterVec(1);
	iterVec(0) = iteration;
	h5writer.writeData(iterVec, "/iteration");

	for (int nf = 0; nf < getNFluids(); nf++) {
		//std::cout << "Face fluxes size: ";
		//std::cout << faceFluxes.rows() << " x " << faceFluxes.cols() << std::endl;
		//assert(nDualFaces == faceFluxes.cols());
		assert(getNFluids() == 1);
		assert(faceFluxes.cols() > 0);
		assert(faceFluxes.rows() > 0);
		assert(nDualFaces == faceFluxes.cols());
		{
			Eigen::VectorXd faceFluxMass(nDualFaces);
			faceFluxMass = faceFluxes.row(0);
			h5writer.writeData(faceFluxMass, "/faceFluxMass" + nf);
		}
		{
			Eigen::MatrixXd faceFluxMomentum(3, nDualFaces);
			faceFluxMomentum = faceFluxes.block(1, 0, 3, nDualFaces);
			h5writer.writeData(faceFluxMomentum, "/faceFluxMomentum" + nf);
		} 
		{
			Eigen::VectorXd faceFluxEnergy(nDualFaces);
			faceFluxEnergy = faceFluxes.row(4);
			h5writer.writeData(faceFluxEnergy, "/faceFluxEnergy" + nf);
		}
	}

	h5writer.writeData(faceFluxesImExRusanov, "/faceFluxesImExRusanov");
}

void AppmSolver::writeFluidStates(H5Writer & writer)
{
	const int nFluids = this->getNFluids();
	const int nCells = dualMesh.getNumberOfCells();

	for (int k = 0; k < nFluids; k++) {
		const std::string fluidTag = (std::stringstream() << "/fluid" << k << "-").str();
		const std::string pressureTag = fluidTag + "pressure";
		const std::string velocityTag = fluidTag + "velocity";
		const std::string densityTag = fluidTag + "density";
		const std::string numberDensityTag = fluidTag + "numberDensity";
		Eigen::VectorXd numberDensity = fluidStates.row(5 * k);
		Eigen::VectorXd density(nCells);
		Eigen::MatrixXd velocity(3, nCells);
		Eigen::VectorXd pressure(nCells);

		const std::string stateN = fluidTag + "stateN";
		const std::string stateU = fluidTag + "stateU";
		const std::string stateE = fluidTag + "stateE";
		Eigen::VectorXd qN(nCells);
		Eigen::MatrixXd qU(3, nCells);
		Eigen::VectorXd qE(nCells);

		const double epsilon2 = particleParams[k].mass;
		for (int i = 0; i < nCells; i++) {
			const Eigen::VectorXd state = fluidStates.block(5 * k, i, 5, 1);

			if (isStateWrittenToOutput) {
				qN(i) = state(0);
				qU.col(i) = state.segment(1, 3);
				qE(i) = state(4);
			}

			const double n = state(0);
			const double rho = epsilon2 * n;
			const Eigen::Vector3d u = epsilon2 * state.segment(1, 3) / n;
			const double p = epsilon2 * (Physics::gamma - 1) * (state(4) - 0.5 * n * u.squaredNorm());

			density(i) = rho;
			velocity.col(i) = u;
			pressure(i) = p;
		}

		writer.writeData(density, densityTag);
		writer.writeData(numberDensity, numberDensityTag);
		writer.writeData(pressure, pressureTag);
		writer.writeData(velocity, velocityTag);

		if (isStateWrittenToOutput) {
			writer.writeData(qN, stateN);
			writer.writeData(qU, stateU);
			writer.writeData(qE, stateE);
		}
	}
}

void AppmSolver::writeMaxwellStates(H5Writer & writer)
{
	writer.writeData(maxwellState, "/x");

	assert(B_h.size() > 0);
	writer.writeData(B_h, "/bvec");

	const int nPrimalFaces = primalMesh.getNumberOfFaces();
	Eigen::Matrix3Xd B(3, nPrimalFaces);
	for (int i = 0; i < nPrimalFaces; i++) {
		const Face * face = primalMesh.getFace(i);
		const Eigen::Vector3d fn = face->getNormal();
		const double fA = face->getArea();
		B.col(i) = (B_h(i) / fA) * fn;
	}
	writer.writeData(B, "/B");

	assert(E_h.size() > 0);
	writer.writeData(E_h, "/evec");

	const int nPrimalEdges = primalMesh.getNumberOfEdges();
	Eigen::Matrix3Xd E(3, nPrimalEdges);
	for (int i = 0; i < nPrimalEdges; i++) {
		const Edge * edge = primalMesh.getEdge(i);
		E.col(i) = E_h(i) / edge->getLength() * edge->getDirection();
	}
	writer.writeData(E, "/E");

	writer.writeData(E_cc, "/Ecc");

	assert(J_h.size() > 0);
	assert(J_h.size() == dualMesh.getNumberOfFaces());
	//assert(J_h.size() == dualMesh.getNumberOfFaces() - 1);
	Eigen::Matrix3Xd currentDensity(3, dualMesh.getNumberOfFaces());
	currentDensity.setZero();
	for (int i = 0; i < J_h.size(); i++) {
		const Face * face = dualMesh.getFace(i);
		const double fA = face->getArea();
		const Eigen::Vector3d fn = face->getNormal();
		currentDensity.col(i) = J_h(i) / fA * fn;
	}
	writer.writeData(currentDensity, "/CurrentDensity");


	assert(dualMesh.getNumberOfFaces() == J_h_aux.size());
	Eigen::Matrix3Xd J_h_aux_vector(3, J_h_aux.size());
	for (int i = 0; i < J_h_aux.size(); i++) {
		const Face * face = dualMesh.getFace(i);
		const double fA = face->getArea();
		const Eigen::Vector3d fn = face->getNormal();
		J_h_aux_vector.col(i) = J_h_aux(i) / fA * fn;
	}
	writer.writeData(J_h_aux_vector, "/J_h_aux_vector");
	writer.writeData(J_h_aux, "/J_h_aux");

}

XdmfGrid AppmSolver::getOutputPrimalEdgeGrid(const int iteration, const double time, const std::string & dataFilename)
{
	XdmfGrid primalEdgeGrid(XdmfGrid::Tags("Primal Edges"));
	XdmfTopology primalEdgeTopology(
		XdmfTopology::Tags(XdmfTopology::TopologyType::Polyline, primalMesh.getNumberOfEdges(), 2)
	);
	primalEdgeTopology.addChild(
		XdmfDataItem(XdmfDataItem::Tags(
			{ 2 * primalMesh.getNumberOfEdges() },
			XdmfDataItem::NumberType::Int,
			XdmfDataItem::Format::HDF),
			(std::stringstream() << primalMesh.getPrefix() << "-mesh.h5:/edge2vertex").str()
		));
	primalEdgeGrid.addChild(primalEdgeTopology);

	XdmfGeometry primalEdgeGeometry;
	primalEdgeGeometry.addChild(
		XdmfDataItem(XdmfDataItem::Tags(
			{ primalMesh.getNumberOfVertices(), 3 },
			XdmfDataItem::NumberType::Float,
			XdmfDataItem::Format::HDF),
			(std::stringstream() << primalMesh.getPrefix() << "-mesh.h5:/vertexPos").str())
	);
	primalEdgeGrid.addChild(primalEdgeGeometry);

	// Attribute: Edge index
	XdmfAttribute edgeIndexAttribute(
		XdmfAttribute::Tags("Edge index", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell)
	);
	edgeIndexAttribute.addChild(
		XdmfDataItem(XdmfDataItem::Tags(
			{ primalMesh.getNumberOfEdges() },
			XdmfDataItem::NumberType::Int,
			XdmfDataItem::Format::HDF),
			(std::stringstream() << primalMesh.getPrefix() << "-mesh.h5:/edgeIdx").str()
		));
	primalEdgeGrid.addChild(edgeIndexAttribute);

	// Attribute: Electric field E
	XdmfAttribute efieldAttribute(
		XdmfAttribute::Tags("Electric field", XdmfAttribute::Type::Vector, XdmfAttribute::Center::Cell)
	);
	efieldAttribute.addChild(
		XdmfDataItem(XdmfDataItem::Tags(
			{ primalMesh.getNumberOfEdges(), 3 },
			XdmfDataItem::NumberType::Float,
			XdmfDataItem::Format::HDF),
			(std::stringstream() << dataFilename << ":/E").str()
		));
	primalEdgeGrid.addChild(efieldAttribute);
	return primalEdgeGrid;
}

XdmfGrid AppmSolver::getOutputPrimalSurfaceGrid(const int iteration, const double time, const std::string & dataFilename)
{
	H5Reader h5reader;
	h5reader = H5Reader("primal-mesh.h5");
	const int nElements = h5reader.readDataSize("/face2vertex");
	assert(nElements > 0);

	XdmfGrid primalSurfaceGrid(XdmfGrid::Tags("Primal Faces"));
	XdmfTopology topology(
		XdmfTopology::Tags(XdmfTopology::TopologyType::Mixed, primalMesh.getNumberOfFaces())
	);
	topology.addChild(
		XdmfDataItem(XdmfDataItem::Tags(
			{ nElements },
			XdmfDataItem::NumberType::Int,
			XdmfDataItem::Format::HDF),
			(std::stringstream() << primalMesh.getPrefix() << "-mesh.h5:/face2vertex").str()
		));
	primalSurfaceGrid.addChild(topology);

	XdmfGeometry primalFaceGeometry;
	primalFaceGeometry.addChild(
		XdmfDataItem(XdmfDataItem::Tags(
			{ primalMesh.getNumberOfVertices(), 3 },
			XdmfDataItem::NumberType::Float,
			XdmfDataItem::Format::HDF),
			(std::stringstream() << primalMesh.getPrefix() << "-mesh.h5:/vertexPos").str())
	);
	primalSurfaceGrid.addChild(primalFaceGeometry);

	// Attribute: Face index
	XdmfAttribute faceIndexAttribute(
		XdmfAttribute::Tags("Face index", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell)
	);
	faceIndexAttribute.addChild(
		XdmfDataItem(XdmfDataItem::Tags(
			{ primalMesh.getNumberOfFaces() },
			XdmfDataItem::NumberType::Int,
			XdmfDataItem::Format::HDF),
			(std::stringstream() << primalMesh.getPrefix() << "-mesh.h5:/faceIndex").str()
		));
	primalSurfaceGrid.addChild(faceIndexAttribute);

	// Attribute: Magnetic flux B
	XdmfAttribute BfieldAttribute(
		XdmfAttribute::Tags("Magnetic Flux", XdmfAttribute::Type::Vector, XdmfAttribute::Center::Cell)
	);
	BfieldAttribute.addChild(
		XdmfDataItem(XdmfDataItem::Tags(
			{ primalMesh.getNumberOfFaces(), 3 },
			XdmfDataItem::NumberType::Float,
			XdmfDataItem::Format::HDF),
			(std::stringstream() << dataFilename << ":/B").str()
		));
	primalSurfaceGrid.addChild(BfieldAttribute);

	{
		// Attribute: Magnetic flux B
		XdmfAttribute attribute(
			XdmfAttribute::Tags("Magnetic Flux Interpolated", XdmfAttribute::Type::Vector, XdmfAttribute::Center::Node)
		);
		attribute.addChild(
			XdmfDataItem(XdmfDataItem::Tags(
				{ primalMesh.getNumberOfVertices(), 3 },
				XdmfDataItem::NumberType::Float,
				XdmfDataItem::Format::HDF),
				(std::stringstream() << dataFilename << ":/Bvertex").str()
			));
		primalSurfaceGrid.addChild(attribute);

	}
	return primalSurfaceGrid;
}

XdmfGrid AppmSolver::getOutputDualEdgeGrid(const int iteration, const double time, const std::string & dataFilename)
{
	XdmfGrid grid(XdmfGrid::Tags("Dual Edges"));
	XdmfTopology topology(
		XdmfTopology::Tags(XdmfTopology::TopologyType::Polyline, dualMesh.getNumberOfEdges(), 2)
	);
	topology.addChild(
		XdmfDataItem(XdmfDataItem::Tags(
			{ 2 * dualMesh.getNumberOfEdges() },
			XdmfDataItem::NumberType::Int,
			XdmfDataItem::Format::HDF),
			(std::stringstream() << dualMesh.getPrefix() << "-mesh.h5:/edge2vertex").str()
		));
	grid.addChild(topology);

	XdmfGeometry geometry;
	geometry.addChild(
		XdmfDataItem(XdmfDataItem::Tags(
			{ dualMesh.getNumberOfVertices(), 3 },
			XdmfDataItem::NumberType::Float,
			XdmfDataItem::Format::HDF),
			(std::stringstream() << dualMesh.getPrefix() << "-mesh.h5:/vertexPos").str())
	);
	grid.addChild(geometry);

	// Attribute: Edge index
	{
		XdmfAttribute attribute(
			XdmfAttribute::Tags("Edge index", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell)
		);
		attribute.addChild(
			XdmfDataItem(XdmfDataItem::Tags(
				{ dualMesh.getNumberOfEdges() },
				XdmfDataItem::NumberType::Int,
				XdmfDataItem::Format::HDF),
				(std::stringstream() << dualMesh.getPrefix() << "-mesh.h5:/edgeIdx").str()
			));
		grid.addChild(attribute);
	}

	// Attribute: Magnetic Field H
	{
		XdmfAttribute attribute(
			XdmfAttribute::Tags("Magnetic field", XdmfAttribute::Type::Vector, XdmfAttribute::Center::Cell)
		);
		attribute.addChild(
			XdmfDataItem(XdmfDataItem::Tags(
				{ dualMesh.getNumberOfEdges(), 3 },
				XdmfDataItem::NumberType::Float,
				XdmfDataItem::Format::HDF),
				(std::stringstream() << dataFilename << ":/H").str()
			));
		grid.addChild(attribute);
	}
	return grid;
}

XdmfGrid AppmSolver::getOutputDualSurfaceGrid(const int iteration, const double time, const std::string & dataFilename)
{
	std::string datafilename = dualMesh.getMeshDataFilename();
	H5Reader h5reader;
	h5reader = H5Reader(datafilename);
	const int nElements = h5reader.readDataSize("/face2vertex");
	assert(nElements > 0);

	XdmfGrid grid(XdmfGrid::Tags("Dual Faces"));
	XdmfTopology topology(
		XdmfTopology::Tags(XdmfTopology::TopologyType::Mixed, dualMesh.getNumberOfFaces())
	);
	topology.addChild(
		XdmfDataItem(XdmfDataItem::Tags(
			{ nElements },
			XdmfDataItem::NumberType::Int,
			XdmfDataItem::Format::HDF),
			(std::stringstream() << dualMesh.getPrefix() << "-mesh.h5:/face2vertex").str()
		));
	grid.addChild(topology);

	XdmfGeometry geometry;
	geometry.addChild(
		XdmfDataItem(XdmfDataItem::Tags(
			{ dualMesh.getNumberOfVertices(), 3 },
			XdmfDataItem::NumberType::Float,
			XdmfDataItem::Format::HDF),
			(std::stringstream() << dualMesh.getPrefix() << "-mesh.h5:/vertexPos").str())
	);
	grid.addChild(geometry);

	// Attribute: Face index
	XdmfAttribute faceIndexAttribute(
		XdmfAttribute::Tags("Face index", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell)
	);
	faceIndexAttribute.addChild(
		XdmfDataItem(XdmfDataItem::Tags(
			{ dualMesh.getNumberOfFaces() },
			XdmfDataItem::NumberType::Int,
			XdmfDataItem::Format::HDF),
			(std::stringstream() << dualMesh.getPrefix() << "-mesh.h5:/faceIndex").str()
		));
	grid.addChild(faceIndexAttribute);

	//// Attribute: Displacement Field D
	//{
	//	XdmfAttribute attribute(
	//		XdmfAttribute::Tags("Displacement Field", XdmfAttribute::Type::Vector, XdmfAttribute::Center::Cell)
	//	);
	//	attribute.addChild(
	//		XdmfDataItem(XdmfDataItem::Tags(
	//			{ dualMesh.getNumberOfFaces(), 3 },
	//			XdmfDataItem::NumberType::Float,
	//			XdmfDataItem::Format::HDF),
	//			(std::stringstream() << dataFilename << ":/D").str()
	//		));
	//	grid.addChild(attribute);
	//}

	// Attribute: Electric current J
	{
		XdmfAttribute attribute(
			XdmfAttribute::Tags("Electric Current", XdmfAttribute::Type::Vector, XdmfAttribute::Center::Cell)
		);
		attribute.addChild(
			XdmfDataItem(XdmfDataItem::Tags(
				{ dualMesh.getNumberOfFaces(), 3 },
				XdmfDataItem::NumberType::Float,
				XdmfDataItem::Format::HDF),
				(std::stringstream() << dataFilename << ":/J").str()
			));
		grid.addChild(attribute);
	}

	return grid;
}

XdmfGrid AppmSolver::getOutputDualVolumeGrid(const int iteration, const double time, const std::string & dataFilename)
{
	std::string datafilename = dualMesh.getMeshDataFilename();
	H5Reader h5reader;
	h5reader = H5Reader(datafilename);

	const int nElements = h5reader.readDataSize("/cell2vertex");
	assert(nElements > 0);

	XdmfGrid grid(XdmfGrid::Tags("Dual Cells"));
	XdmfTopology topology(
		XdmfTopology::Tags(XdmfTopology::TopologyType::Mixed, dualMesh.getNumberOfCells())
	);
	topology.addChild(
		XdmfDataItem(XdmfDataItem::Tags(
			{ nElements },
			XdmfDataItem::NumberType::Int,
			XdmfDataItem::Format::HDF),
			(std::stringstream() << datafilename << ":/cell2vertex").str()
		));
	grid.addChild(topology);

	XdmfGeometry geometry;
	geometry.addChild(
		XdmfDataItem(XdmfDataItem::Tags(
			{ dualMesh.getNumberOfVertices(), 3 },
			XdmfDataItem::NumberType::Float,
			XdmfDataItem::Format::HDF),
			(std::stringstream() << datafilename << ":/vertexPos").str())
	);
	grid.addChild(geometry);

	// Attribute: Cell index
	XdmfAttribute cellIndexAttribute(
		XdmfAttribute::Tags("Cell index", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell)
	);
	cellIndexAttribute.addChild(
		XdmfDataItem(XdmfDataItem::Tags(
			{ dualMesh.getNumberOfCells() },
			XdmfDataItem::NumberType::Int,
			XdmfDataItem::Format::HDF),
			(std::stringstream() << datafilename << ":/cellIndex").str()
		));
	grid.addChild(cellIndexAttribute);

	// Attribute: B-field at primal vertices = dual cell centers
	{
		XdmfAttribute attribute(
			XdmfAttribute::Tags("Magnetic Flux", XdmfAttribute::Type::Vector, XdmfAttribute::Center::Node)
		);
		attribute.addChild(
			XdmfDataItem(XdmfDataItem::Tags(
				{ dualMesh.getNumberOfVertices(), 3 },
				XdmfDataItem::NumberType::Float,
				XdmfDataItem::Format::HDF),
				(std::stringstream() << dataFilename << ":/Bvertex").str()
			));
		grid.addChild(attribute);
	}

	// Attribute: Density
	XdmfAttribute densityAttribute(
		XdmfAttribute::Tags("Density", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell)
	);
	densityAttribute.addChild(
		XdmfDataItem(XdmfDataItem::Tags(
			{ dualMesh.getNumberOfCells() },
			XdmfDataItem::NumberType::Float,
			XdmfDataItem::Format::HDF),
			(std::stringstream() << dataFilename << ":/density").str()
		));
	grid.addChild(densityAttribute);

	// Attribute: Pressure
	XdmfAttribute pressureAttribute(
		XdmfAttribute::Tags("Pressure", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell)
	);
	pressureAttribute.addChild(
		XdmfDataItem(XdmfDataItem::Tags(
			{ dualMesh.getNumberOfCells() },
			XdmfDataItem::NumberType::Float,
			XdmfDataItem::Format::HDF),
			(std::stringstream() << dataFilename << ":/pressure").str()
		));
	grid.addChild(pressureAttribute);

	// Attribute: velocity
	XdmfAttribute velocityAttribute(
		XdmfAttribute::Tags("Velocity", XdmfAttribute::Type::Vector, XdmfAttribute::Center::Cell)
	);
	velocityAttribute.addChild(
		XdmfDataItem(XdmfDataItem::Tags(
			{ dualMesh.getNumberOfCells(), 3 },
			XdmfDataItem::NumberType::Float,
			XdmfDataItem::Format::HDF),
			(std::stringstream() << dataFilename << ":/velocity").str()
		));
	grid.addChild(velocityAttribute);




	return grid;
}




void AppmSolver::init_RaviartThomasInterpolation()
{
	const int nCells = primalMesh.getNumberOfCells();
	rt_piolaMatrix.reserve(nCells);
	rt_piolaVector.reserve(nCells);

	Eigen::VectorXi isPiolaMapDefined(nCells);
	isPiolaMapDefined.setZero();

	bool isInfoPrinted = true;

	for (int cidx = 0; cidx < nCells; cidx++) {
		const Cell * cell = primalMesh.getCell(cidx);
		const std::vector<Face*> cellFaces = cell->getFaceList();
		if (cellFaces.size() == 5) {
			isPiolaMapDefined(cidx) = 1;
		}
		else {
			isPiolaMapDefined(cidx) = 0;
			if (isInfoPrinted) {
				std::cout << "Raviart-Thomas interpolation is only implemented for triangular prisms!" << std::endl;
				isInfoPrinted = false; // show info only once
			}
			continue;
		}
		assert(cellFaces.size() == 5); // prism cells

		// get top and bottom face of prism
		const Face * bottomFace = cellFaces[3];
		const Face * topFace = cellFaces[4];
		assert(bottomFace->getVertexList().size() == 3); // check if faces are triangular
		assert(topFace->getVertexList().size() == 3);

		// get vertices of triangle faces
		const std::vector<Vertex*> bottomVertices = bottomFace->getVertexList();
		const std::vector<Vertex*> topVertices = topFace->getVertexList();

		// check if vertices of triangle faces have same ordering
		for (int i = 0; i < 3; i++) {
			const Eigen::Vector3d v = topVertices[i]->getPosition() - bottomVertices[i]->getPosition();
			assert(v.normalized().cross(Eigen::Vector3d::UnitZ()).norm() < 16 * std::numeric_limits<double>::epsilon());
		}

		// check if bottom triangle vertices form a right-handed system
		const Eigen::Vector3d fc = bottomFace->getCenter();
		for (int i = 0; i < 3; i++) {
			const Eigen::Vector3d v0 = bottomVertices[i]->getPosition();
			const Eigen::Vector3d v1 = bottomVertices[(i + 1) % 3]->getPosition();
			const Eigen::Vector3d a = v0 - fc;
			const Eigen::Vector3d b = v1 - v0;
			const Eigen::Vector3d n = a.normalized().cross(b.normalized());
			assert(n.dot(Eigen::Vector3d::UnitZ()) > 0);
		}

		// Collect prism vertices in standard topology format
		std::vector<Vertex*> cellVertices(6);
		for (int i = 0; i < 3; i++) {
			const int offset = 3;
			cellVertices[i] = bottomVertices[i];
			cellVertices[i + offset] = topVertices[i];
		}

		// Vertex positions of prism
		const Eigen::Vector3d A = cellVertices[0]->getPosition();
		const Eigen::Vector3d B = cellVertices[1]->getPosition();
		const Eigen::Vector3d C = cellVertices[2]->getPosition();
		const Eigen::Vector3d F = cellVertices[5]->getPosition();

		// Define Piola map: xRef  -->  BK * xRef + bK
		const Eigen::Vector3d bK = C;
		Eigen::Matrix3d BK;
		BK.col(0) = A - C;
		BK.col(1) = B - C;
		BK.col(2) = F - C;
		rt_piolaMatrix.emplace_back(BK);
		rt_piolaVector.emplace_back(bK);
	}

	// Check if all Piola maps are defined at beginning of list
	const int nPiolaMapsDefined = isPiolaMapDefined.count();
	//std::cout << "nPiolaMapsDefined: " << nPiolaMapsDefined << std::endl;
	//std::ofstream file("isPiolaMapDefined.dat");
	//file << isPiolaMapDefined << std::endl;
	assert(isPiolaMapDefined.topRows(nPiolaMapsDefined).all());
	assert(isPiolaMapDefined.bottomRows(nCells - nPiolaMapsDefined).any() == 0);

	// Truncate list of Piola maps
	rt_piolaMatrix.resize(nPiolaMapsDefined);
	rt_piolaVector.resize(nPiolaMapsDefined);
}

void AppmSolver::readParameters(const std::string & filename)
{
	std::ifstream file(filename);
	if (!file.is_open()) {
		std::cout << "File not opened: " << filename;
		exit(-1);
	}

	std::string line;
	const char delim = ':';

	while (std::getline(file, line)) {
		int pos = line.find(delim);
		std::string tag = line.substr(0, pos);

		if (tag == "maxIterations") {
			std::istringstream(line.substr(pos + 1)) >> this->maxIterations;
		}
		if (tag == "maxTime") {
			std::istringstream(line.substr(pos + 1)) >> this->maxTime;
		}
		if (tag == "isFluidEnabled") {
			std::istringstream(line.substr(pos + 1)) >> this->isFluidEnabled;
		}
		if (tag == "isMaxwellEnabled") {
			std::istringstream(line.substr(pos + 1)) >> this->isMaxwellEnabled;
		}
		if (tag == "timestepSize") {
			std::istringstream(line.substr(pos + 1)) >> this->timestepSize;
		}
		if (tag == "lambdaSquare") {
			std::istringstream(line.substr(pos + 1)) >> this->lambdaSquare;
		}
		if (tag == "massFluxScheme") {
			std::string temp;
			std::istringstream(line.substr(pos + 1)) >> temp;
			if (temp == "ImEx") {
				massFluxScheme = MassFluxScheme::IMPLICIT_EXPLICIT;
			}
		}
		if (tag == "initType") {
			std::istringstream(line.substr(pos + 1)) >> initType;
		}
		if (tag == "isElectricLorentzForceActive") {
			std::istringstream(line.substr(pos + 1)) >> isElectricLorentzForceActive;
		}
		if (tag == "isMagneticLorentzForceActive") {
			std::istringstream(line.substr(pos + 1)) >> isMagneticLorentzForceActive;
		}
		//if (tag == "maxwellSolverBCType") {
		//	std::string temp;
		//	std::istringstream(line.substr(pos + 1)) >> temp;
		//	if (temp == "Voltage") {
		//		maxwellSolverBCType = MaxwellSolverBCType::VOLTAGE_BC;
		//	}
		//	if (temp == "Current") {
		//		maxwellSolverBCType = MaxwellSolverBCType::CURRENT_BC;
		//	}
		//}

	}

	std::cout << std::endl;
	std::cout << "Appm Solver parameters:" << std::endl;
	std::cout << "======================="  << std::endl;
	std::cout << "maxIterations:  " << maxIterations << std::endl;
	std::cout << "maxTime:        " << maxTime << std::endl;
	std::cout << "isFluidEnabled: " << isFluidEnabled << std::endl;
	std::cout << "isMaxwellEnabled: " << isMaxwellEnabled << std::endl;
	std::cout << "timestepSize: " << timestepSize << std::endl;
	std::cout << "lambdaSquare: " << lambdaSquare << std::endl;
	std::cout << "massFluxScheme: " << massFluxScheme << std::endl;
	std::cout << "initType: " << initType << std::endl;
	std::cout << "isElectricLorentzForceActive: " << isElectricLorentzForceActive << std::endl;
	std::cout << "isMagneticLorentzForceActive: " << isMagneticLorentzForceActive << std::endl;
	//std::cout << "maxwellSolverBCType: " << maxwellSolverBCType << std::endl;
	std::cout << "=======================" << std::endl;
}

const std::string AppmSolver::xdmf_GridPrimalEdges(const int iteration) const
{
	const std::string datafilename = primalMesh.getMeshDataFilename();
	std::stringstream ss;
	ss << "<Grid Name=\"Primal Edges\">" << std::endl;
	ss << "<Topology TopologyType=\"Polyline\""
		<< " NumberOfElements=\"" << primalMesh.getNumberOfEdges() << "\""
		<< " NodesPerElement=\"2\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << 2 * primalMesh.getNumberOfEdges() << "\" DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << datafilename << ":/edge2vertex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Topology>" << std::endl;

	ss << "<Geometry GeometryType=\"XYZ\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << primalMesh.getNumberOfVertices() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << datafilename << ":/vertexPos" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Geometry>" << std::endl;

	ss << "<Attribute Name=\"Edge index\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << primalMesh.getNumberOfEdges() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << datafilename << ":/edgeIdx" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	if (isWriteEfield) {
		ss << "<Attribute Name=\"Electric field\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << primalMesh.getNumberOfEdges() << " 3\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << "appm-" << iteration << ".h5:/E" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;
	}

	ss << "</Grid>";
	return ss.str();
}
const std::string AppmSolver::xdmf_GridPrimalFaces(const int iteration) const
{
	const std::string datafilename = primalMesh.getMeshDataFilename();
	H5Reader h5reader;
	h5reader = H5Reader(datafilename);
	const int nElements = h5reader.readDataSize("/face2vertex");
	assert(nElements > 0);

	std::stringstream ss;
	ss << "<Grid Name=\"Primal Faces\">" << std::endl;
	ss << "<Topology TopologyType=\"Mixed\""
		<< " NumberOfElements=\"" << primalMesh.getNumberOfFaces() << "\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << nElements << "\" DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << datafilename << ":/face2vertex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Topology>" << std::endl;

	ss << "<Geometry GeometryType=\"XYZ\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << primalMesh.getNumberOfVertices() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << datafilename << ":/vertexPos" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Geometry>" << std::endl;

	ss << "<Attribute Name=\"Face index\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << primalMesh.getNumberOfFaces() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << datafilename << ":/faceIndex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	if (isWriteBfield) {
		ss << "<Attribute Name=\"Magnetic flux\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << primalMesh.getNumberOfFaces() << " 3\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << "appm-" << iteration << ".h5:/B" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;
	}

	ss << "<Attribute Name=\"Magnetic flux interpolated\" AttributeType=\"Vector\" Center=\"Node\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << primalMesh.getNumberOfVertices() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "appm-" << iteration << ".h5:/Bvertex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;
	ss << "</Grid>";
	return ss.str();
}

const std::string AppmSolver::xdmf_GridDualEdges(const int iteration) const
{
	const std::string datafilename = dualMesh.getMeshDataFilename();
	std::stringstream ss;
	ss << "<Grid Name=\"Dual Edges\">" << std::endl;
	ss << "<Topology TopologyType=\"Polyline\""
		<< " NumberOfElements=\"" << dualMesh.getNumberOfEdges() << "\"" 
		<< " NodesPerElement=\"2\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << 2 * dualMesh.getNumberOfEdges() << "\" DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << datafilename << ":/edge2vertex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Topology>" << std::endl;

	ss << "<Geometry GeometryType=\"XYZ\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfVertices() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << datafilename << ":/vertexPos" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Geometry>" << std::endl;

	ss << "<Attribute Name=\"Edge index\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfEdges() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << datafilename << ":/edgeIdx" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	if (isWriteHfield) {
		ss << "<Attribute Name=\"Magnetic field\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfEdges() << " 3\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << "appm-" << iteration << ".h5:/H" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;
	}
	ss << "</Grid>" << std::endl;
	return ss.str();
}

const std::string AppmSolver::xdmf_GridDualFaces(const int iteration) const
{
	const std::string datafilename = dualMesh.getMeshDataFilename();
	H5Reader h5reader;
	h5reader = H5Reader(datafilename);
	const int nElements = h5reader.readDataSize("/face2vertex");
	assert(nElements > 0);

	std::stringstream ss;
	ss << "<Grid Name=\"Dual Faces\">" << std::endl;
	ss << "<Topology TopologyType=\"Mixed\""
		<< " NumberOfElements=\"" << dualMesh.getNumberOfFaces() << "\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << nElements << "\" DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << datafilename << ":/face2vertex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Topology>" << std::endl;

	ss << "<Geometry GeometryType=\"XYZ\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfVertices() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << datafilename << ":/vertexPos" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Geometry>" << std::endl;

	ss << "<Attribute Name=\"Face index\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << datafilename << ":/faceIndex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;
	
	ss << "<Attribute Name=\"Face Fluid Type\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << datafilename << ":/faceFluidType" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"Face Flux Mass\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << "\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "appm-" << iteration << ".h5:/faceFluxMass" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"Face Flux Momentum\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "appm-" << iteration << ".h5:/faceFluxMomentum" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"Face Flux Energy\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << "\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "appm-" << iteration << ".h5:/faceFluxEnergy" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"Current density\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "appm-" << iteration << ".h5:/CurrentDensity" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"J_h_aux_Vector\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "appm-" << iteration << ".h5:/J_h_aux_vector" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "</Grid>" << std::endl;
	return ss.str();
}

const std::string AppmSolver::xdmf_GridDualCells(const int iteration) const
{
	const std::string meshDataFilename = dualMesh.getMeshDataFilename();
	H5Reader h5reader;
	h5reader = H5Reader(meshDataFilename);
	const int nElements = h5reader.readDataSize("/cell2vertex");
	assert(nElements > 0);

	std::stringstream ss;
	ss << "<Grid Name=\"Dual Cells\">" << std::endl;
	ss << "<Topology TopologyType=\"Mixed\""
		<< " NumberOfElements=\"" << dualMesh.getNumberOfCells() << "\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << nElements << "\" DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << meshDataFilename << ":/cell2vertex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Topology>" << std::endl;

	ss << "<Geometry GeometryType=\"XYZ\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfVertices() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << meshDataFilename << ":/vertexPos" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Geometry>" << std::endl;

	ss << "<Attribute Name=\"Cell index\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfCells() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << meshDataFilename << ":/cellIndex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"Cell Type\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfCells() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << meshDataFilename << ":/cellFluidType" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	const std::string datafilename = (std::stringstream() << "appm-" << iteration << ".h5").str();

	ss << "<Attribute Name=\"Magnetic Flux Interpolated\" AttributeType=\"Vector\" Center=\"Node\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfVertices() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << datafilename << ":/Bvertex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"Ecc\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfCells() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << datafilename << ":/Ecc" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;



	//ss << fluidSolver->getXdmfOutput(iteration);
	ss << fluidXdmfOutput(datafilename) << std::endl;
	ss << "</Grid>";
	return ss.str();
}

/** 
* Print XDMF output for all fluids. 
* @param datafilename    Data filename (HDF5)
*/
const std::string AppmSolver::fluidXdmfOutput(const std::string & datafilename) const
{
	std::stringstream ss;
	const int nFluids = this->getNFluids();
	const int nCells = dualMesh.getNumberOfCells();
	for (int k = 0; k < nFluids; k++) {
		const std::string fluidName = (std::stringstream() << "" << k).str();

		ss << "<Attribute Name=\"Density " << fluidName << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << nCells << "\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << datafilename << ":/" << "fluid" << k << "-density" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;

		ss << "<Attribute Name=\"Number Density " << fluidName << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << nCells << "\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << datafilename << ":/" << "fluid" << k << "-numberDensity" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;

		ss << "<Attribute Name=\"Pressure " << fluidName << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << nCells << "\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << datafilename << ":/" << "fluid" << k << "-pressure" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;


		ss << "<Attribute Name=\"Velocity " << fluidName << "\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << nCells << " 3\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << datafilename << ":/" << "fluid" << k << "-velocity" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;

		if (isStateWrittenToOutput) {
			ss << "<Attribute Name=\"stateN " << fluidName << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
			ss << "<DataItem Dimensions=\"" << nCells << "\""
				<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
			ss << datafilename << ":/" << "fluid" << k << "-stateN" << std::endl;
			ss << "</DataItem>" << std::endl;
			ss << "</Attribute>" << std::endl;


			ss << "<Attribute Name=\"stateU " << fluidName << "\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
			ss << "<DataItem Dimensions=\"" << nCells << " 3\""
				<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
			ss << datafilename << ":/" << "fluid" << k << "-stateU" << std::endl;
			ss << "</DataItem>" << std::endl;
			ss << "</Attribute>" << std::endl;

			ss << "<Attribute Name=\"stateE " << fluidName << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
			ss << "<DataItem Dimensions=\"" << nCells << "\""
				<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
			ss << datafilename << ":/" << "fluid" << k << "-stateE" << std::endl;
			ss << "</DataItem>" << std::endl;
			ss << "</Attribute>" << std::endl;
		}
	}
	return ss.str();
}

std::ostream & operator<<(std::ostream & os, const AppmSolver::MassFluxScheme & obj) {
	switch (obj) {
	case AppmSolver::MassFluxScheme::EXPLICIT:
		os << "EXPLICIT";
		break;

	case AppmSolver::MassFluxScheme::IMPLICIT_EXPLICIT:
		os << "IMEX";
	}
	return os;
}

//std::ostream & operator<<(std::ostream & os, const AppmSolver::MaxwellSolverBCType & obj) {
//	switch (obj) {
//	case AppmSolver::MaxwellSolverBCType::CURRENT_BC:
//		os << "CURRENT_BC";
//		break;
//	case AppmSolver::MaxwellSolverBCType::VOLTAGE_BC:
//		os << "VOLTAGE_BC";
//		break;
//	default:
//		os << "Not Implemented";
//		assert(false);
//	}
//	return os;
//}
