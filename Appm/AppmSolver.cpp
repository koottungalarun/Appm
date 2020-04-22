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

	init_maxwellStates();

	B_vertex = Eigen::Matrix3Xd::Zero(3, primalMesh.getNumberOfVertices());
	init_RaviartThomasInterpolation();

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

	if (isMaxwellEnabled) {
		isWriteBfield = true;
		isWriteEfield = true;
		//Eigen::sparseMatrixToFile(Q, "Q.dat");
		//Eigen::sparseMatrixToFile(Meps, "Meps.dat");
		//Eigen::sparseMatrixToFile(Mnu, "Mnu.dat");
		//Eigen::sparseMatrixToFile(C, "C.dat");
		//Eigen::sparseMatrixToFile(M1, "M1.dat");
		//Eigen::sparseMatrixToFile(M2, "M2.dat");
		test_implicitEfieldToCurrent();
	}
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

	
	// Time integration loop
	//double dT = 0.05;
	const int nPrimalTerminalVertices = primalMesh.getMeshInfo().nVerticesTerminal;
	assert(nPrimalTerminalVertices > 0);

	const int nFluids = this->getNFluids();
	const int nFaces = dualMesh.getNumberOfFaces();
	faceFluxes = Eigen::MatrixXd::Zero(5*nFluids, nFaces);
	faceFluxesImExRusanov = Eigen::MatrixXd::Zero(nFaces, nFluids);
	writeOutput(iteration, time);

	while (iteration < maxIterations && time < maxTime) {
		std::cout << "Iteration " << iteration << ",\t time = " << time << std::endl;

		dt_previous = dt;

		// Determine timestep
		if (isFluidEnabled) {
			dt = getNextFluidTimestepSize();
		}
		fluidFluxes.setZero();
		fluidSources.setZero();
		faceFluxes.setZero();
		faceFluxesImExRusanov.setZero();

		// Maxwell equations
		if (isMaxwellEnabled) {
			//maxwellSolver->updateMaxwellState(dt, time);
			Eigen::VectorXd x(maxwellState.size());

			const int nEdges = primalMesh.getNumberOfEdges();
			Eigen::SparseMatrix<double> M_sigma(nEdges, nEdges);

			Eigen::SparseMatrix<double> M;
			M = M1 + pow(dt, 2) * M2;
			M.makeCompressed();
			//Eigen::sparseMatrixToFile(M, "M.dat");

			const int nDirichlet = nPrimalTerminalVertices;
			const int nFree = maxwellState.size() - nDirichlet;

			// Dirichlet conditions
			Eigen::SparseMatrix<double> Md = M.rightCols(nDirichlet);
			Eigen::VectorXd xd = Eigen::VectorXd::Zero(nDirichlet);
			assert(nDirichlet % 2 == 0);
			xd.topRows(xd.size() / 2).array() = 2;
			xd.bottomRows(xd.size() / 2).array() = 1; // voltage condition on terminal B

			Eigen::VectorXd src = Eigen::VectorXd::Zero(maxwellState.size());
			Eigen::VectorXd rhs(x.size());
			rhs.setZero();
			double dt_ratio = dt / dt_previous;
			rhs += M1 * ((1 + dt_ratio) * maxwellState - dt_ratio * maxwellStatePrevious);
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
			Eigen::SparseLU<Eigen::SparseMatrix<double>> maxwellSolver(Mf);
			if (maxwellSolver.info() != Eigen::Success) {
				std::cout << "Maxwell solver setup failed" << std::endl;
			}
			xf = maxwellSolver.solve(rhsFree);
			std::cout << "Maxwell solver finished" << std::endl;
			if (maxwellSolver.info() != Eigen::Success) {
				std::cout << "Maxwell solver solving failed" << std::endl;
			}
			assert(maxwellSolver.info() == Eigen::Success);

			// assemble new state vector
			x.topRows(nFree) = xf;
			x.bottomRows(nDirichlet) = xd;
			//std::ofstream("x.dat") << x << std::endl;

			// Update state vectors for discrete states
			E_h = Q * x;
			B_h -= dt * C * E_h;

			J_h = M_sigma * E_h;

			//std::ofstream("E_h.dat") << E_h << std::endl;

			// Set new state vector
			maxwellStatePrevious = maxwellState;
			maxwellState = x;

			interpolateElectricFieldToCellCenter();

			// Interpolation of B-field to dual cell centers
			interpolateMagneticFluxToPrimalVertices();

			// Consistent formulation of current as mass flux

			// Solve  d^2/dt^2 (E) + curl(curl(E)) + n E = rhs, for new timestep
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
						assert(false); // TODO to be checked
						//assert(faceCells.size() == 1);
						//const Cell * cell = faceCells[0];
						//assert(cell->getFluidType() == Cell::FluidType::FLUID);
						//int idxL = cell->getIndex();
						//fluidFluxes.col(idxL).segment(5 * fluidIdx, 5) += faceFlux3d;
						//assert(false);
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
			if (time < 0) {
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

						// Define geometric position of source region: ball of radius r around reference position
						const Eigen::Vector3d cc = cell->getCenter();
						const Eigen::Vector3d srcRefPos(0, 0, 0);
						const double srcRadius = 0.2;

						if ((cc - srcRefPos).norm() < srcRadius) {
							srcLocal(0) = 0; // mass source
							srcLocal.segment(1, 3) = Eigen::Vector3d::Zero(); // momentum source
							srcLocal(4) = 1; // energy source
						}
						fluidSources.col(cIdx).segment(5 * fluidIdx, 5) = srcLocal;
					}
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
	}
	std::cout << "Final time:      " << time << std::endl;
	std::cout << "Final iteration: " << iteration << std::endl;

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

	electricFieldAtDualCellCenters = Eigen::Matrix3Xd::Zero(3, dualMesh.getNumberOfCells());
	
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
Test to implement Ohm's law j = M_sigma * e, as given by the implicit and consistent formulation of current and mass flux.
*/
void AppmSolver::test_implicitEfieldToCurrent()
{
	std::cout << "Test for Ohms law" << std::endl;
	
	// Set E_h uniform in z-direction
	E_h.setZero();
	assert(E_h.size() == primalMesh.getNumberOfEdges());
	for (int i = 0; i < E_h.size(); i++) {
		const Edge * edge = primalMesh.getEdge(i);
		const Eigen::Vector3d edgeDir = edge->getDirection().normalized();
		const double edgeLength = edge->getLength();
		E_h(i) = edgeDir.dot(Eigen::Vector3d::UnitZ()) * edgeLength;
	}

	// Setup of equation J = M * E + Jb
	assert(J_h.size() > E_h.size());
	Eigen::SparseMatrix<double> M_sigma(J_h.size(), E_h.size());
	Eigen::VectorXd J_h_b(J_h.size());
	J_h_b.setZero();

	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;

	// For each dual face
	const int refIdx = 7565;
	for (int dualFaceIdx = 0; dualFaceIdx < J_h.size(); dualFaceIdx++) {

		//if (dualFaceIdx != refIdx) { continue; }

		const Face * dualFace = dualMesh.getFace(dualFaceIdx);
		const Eigen::Vector3d dualFaceNormal = dualFace->getNormal();
		const bool isDualFaceBoundary = dualFace->isBoundary();

		// TODO remove this limiter
		//if (dualFaceIdx >= primalMesh.getNumberOfEdges()) { break; }

		// TODO only for faces in z-direction
		if (dualFaceNormal.cross(Eigen::Vector3d::UnitZ()).norm() > 1e-8) {
			continue;
		}

		// get list of adjacient cells
		const std::vector<Cell*> adjacientCells = dualFace->getCellList();
		assert(adjacientCells.size() >= 1);
		assert(adjacientCells.size() <= 2);
		assert((adjacientCells.size() == 1 &&  dualFace->isBoundary())
			|| (adjacientCells.size() == 2 && !dualFace->isBoundary()));

		// For each adjacient cell
		for (auto cell : adjacientCells) {
			const double cellVolume = cell->getVolume();
			std::vector<Face*> cellFaces = cell->getFaceList();
			const Eigen::Vector3d cc = cell->getCenter();

			// For each face of that cell
			for (auto face : cellFaces) {
				int faceIdx = face->getIndex();
				const double faceArea = face->getArea();
				const Eigen::Vector3d fn = face->getNormal().normalized();
				const Eigen::Vector3d fc = face->getCenter();
				const Eigen::Vector3d r = fc - cc;

				// Get corresponding primal edge
				const int edgeIdx = dualMesh.getAssociatedPrimalEdgeIndex(faceIdx);
				assert(edgeIdx >= 0);
				assert(edgeIdx < primalMesh.getNumberOfEdges());
				const Edge * edge = primalMesh.getEdge(edgeIdx);
				const double edgeLength = edge->getLength();
				const Eigen::Vector3d edgeDirection = edge->getDirection().normalized();

				// assert that edge direction and face normal are (almost) parallel
				const double tol = 0.2;
				const Eigen::Vector3d e_x_fn = edgeDirection.cross(fn);
				const double e_x_fn_norm = e_x_fn.norm();
				const bool isParallel = e_x_fn_norm < tol;
				assert(isParallel);

				// 
				const int row = dualFaceIdx;
				const int col = edgeIdx;
				const double n = 1; // number density in cell 
				assert(cellVolume > 0);
				const double r_dot_dualFaceNormal = r.dot(dualFaceNormal);
				
				// edgeLength * ((face->isBoundary()) ? 0.5 : 1); // this factor is from the hydrodynamic boundary condition
				
				const Eigen::Vector3d nj = face->getNormal();
				const Eigen::Vector3d L = edge->getDirection().normalized();
				const int incidence = r.dot(fn) > 0 ? 1 : -1;
				const double dL = edge->getLength();
				const double dA = face->getArea();
				const double dV = cell->getVolume();
				double value = n * r_dot_dualFaceNormal * 1. / dL * (L.dot(nj)) * incidence * dA / dV;
				//value *= cellVolume; // TODO 

				value *= (dualFace->isBoundary()) ? 1 : 0.5;
				value *= dualFace->getArea();
				//value *= (face->isBoundary()) ? 2 : 1; // this factor is from the mesh construction: dual boundary cells have half-length of inner cells

				if (row == refIdx && E_h(col) != 0) {
					std::cout << std::scientific << std::setprecision(4);
					std::cout << std::endl;
					std::cout << "row: " << row << "\t";
					std::cout << "col: " << col << "\t";
					std::cout << "value: " << value << std::endl;

					std::cout << "dV:\t" << dV << std::endl;
					std::cout << "dA:\t" << dA << std::endl;
					std::cout << "dA/dV:\t" << dA / dV << std::endl;
					std::cout << "dL:\t" << dL << std::endl;

					std::cout << "fn:\t" << fn.transpose() << std::endl;
					std::cout << "fc:\t" << fc.transpose() << std::endl;
					std::cout << "cc:\t" << cc.transpose() << std::endl;
					std::cout << "ri:\t" << r.transpose() << std::endl;

					std::cout << "E_h(col):\t" << E_h(col) << std::endl;
					std::cout << "L:\t" << L.transpose() << std::endl;
					std::cout << "nj:\t" << nj.transpose() << std::endl;
					std::cout << "L.nj:\t" << L.dot(nj) << std::endl;
					std::cout << "r.dot(fn):\t" << incidence << std::endl;

					std::cout << std::endl;
				}

				triplets.push_back(T(row, col, value));
			}
		}
	}
	std::cout << std::endl;

	M_sigma.setFromTriplets(triplets.begin(), triplets.end());
	M_sigma.makeCompressed();

	// TODO: multiply with electric ionization degree, mass ratio, and timestep size
	//const double dt = 1;
	//const double q = 1;
	//const double massRatio = 1;
	//M_sigma *= dt * q * massRatio;

	J_h = M_sigma * E_h + J_h_b;
	
	std::cout << "J_h(" << refIdx << "): \t";
	std::cout << J_h(refIdx) << std::endl;
	std::cout << "J_h/A: \t";
	std::cout << J_h(refIdx) / dualMesh.getFace(refIdx)->getArea() << std::endl;

	Eigen::sparseMatrixToFile(M_sigma, "Msigma.dat");
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
		assert(faceCells.size() >= 1);
		assert(false);
		exit(-1);
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

/**
* Interpolate electric field from primal edges to dual cell centers with Perot's method.
*/
void AppmSolver::interpolateElectricFieldToCellCenter()
{
	std::cout << "Interpolate electric field to cell center (Perot method)" << std::endl;
	assert(primalMesh.getNumberOfVertices() == dualMesh.getNumberOfCells());
	const int nDualCells = dualMesh.getNumberOfCells();
	const int nPe = primalMesh.getNumberOfEdges();

	for (int k = 0; k < nDualCells; k++) {
		Eigen::Vector3d E_loc;
		E_loc.setZero();

		const Cell * dualCell = dualMesh.getCell(k);
		const std::vector<Face*> cellFaces = dualCell->getFaceList();
		for (auto dualFace : cellFaces) {
			const int i = dualFace->getIndex();
			if (i >= nPe) { continue; }			
			const double faceArea = dualFace->getArea();
			const Eigen::Vector3d r = dualFace->getCenter() - dualCell->getCenter();
			const Eigen::Vector3d fn = dualFace->getNormal();
			const int orientation = getOrientation(dualCell, dualFace);
			const double primalEdgeLength = primalMesh.getEdge(i)->getLength();
			const double eValue = E_h(i) * orientation * 1./primalEdgeLength;
			E_loc += eValue * faceArea * r;
		}
		E_loc *= 1. / dualCell->getVolume();
		electricFieldAtDualCellCenters.col(k) = E_loc;
	}
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
	Eigen::VectorXd B_h(nPrimalFaces); 
	B_h.setZero();

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

	assert(electricFieldAtDualCellCenters.size() > 0);
	writer.writeData(electricFieldAtDualCellCenters, "/EdualCellCenter");

	assert(J_h.size() > 0);
	assert(J_h.size() == dualMesh.getNumberOfFaces());
	//assert(J_h.size() == dualMesh.getNumberOfFaces() - 1);
	Eigen::Matrix3Xd J_consistent(3, dualMesh.getNumberOfFaces());
	J_consistent.setZero();
	for (int i = 0; i < J_h.size(); i++) {
		const Face * face = dualMesh.getFace(i);
		const double fA = face->getArea();
		const Eigen::Vector3d fn = face->getNormal();
		J_consistent.col(i) = J_h(i) / fA * fn;
	}
	writer.writeData(J_consistent, "/J_consistent");
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
	H5Reader h5reader;
	h5reader = H5Reader("dual-mesh.h5");
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
	H5Reader h5reader;
	h5reader = H5Reader("dual-mesh.h5");
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
			(std::stringstream() << dualMesh.getPrefix() << "-mesh.h5:/cell2vertex").str()
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

	// Attribute: Cell index
	XdmfAttribute cellIndexAttribute(
		XdmfAttribute::Tags("Cell index", XdmfAttribute::Type::Scalar, XdmfAttribute::Center::Cell)
	);
	cellIndexAttribute.addChild(
		XdmfDataItem(XdmfDataItem::Tags(
			{ dualMesh.getNumberOfCells() },
			XdmfDataItem::NumberType::Int,
			XdmfDataItem::Format::HDF),
			(std::stringstream() << dualMesh.getPrefix() << "-mesh.h5:/cellIndex").str()
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
	std::cout << "=======================" << std::endl;
}

const std::string AppmSolver::xdmf_GridPrimalEdges(const int iteration) const
{
	std::stringstream ss;
	ss << "<Grid Name=\"Primal Edges\">" << std::endl;
	ss << "<Topology TopologyType=\"Polyline\""
		<< " NumberOfElements=\"" << primalMesh.getNumberOfEdges() << "\""
		<< " NodesPerElement=\"2\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << 2 * primalMesh.getNumberOfEdges() << "\" DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << "primal-mesh.h5:/edge2vertex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Topology>" << std::endl;

	ss << "<Geometry GeometryType=\"XYZ\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << primalMesh.getNumberOfVertices() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "primal-mesh.h5:/vertexPos" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Geometry>" << std::endl;

	ss << "<Attribute Name=\"Edge index\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << primalMesh.getNumberOfEdges() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << "primal-mesh.h5:/edgeIdx" << std::endl;
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
	H5Reader h5reader;
	h5reader = H5Reader("primal-mesh.h5");
	const int nElements = h5reader.readDataSize("/face2vertex");
	assert(nElements > 0);

	std::stringstream ss;
	ss << "<Grid Name=\"Primal Faces\">" << std::endl;
	ss << "<Topology TopologyType=\"Mixed\""
		<< " NumberOfElements=\"" << primalMesh.getNumberOfFaces() << "\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << nElements << "\" DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << "primal-mesh.h5:/face2vertex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Topology>" << std::endl;

	ss << "<Geometry GeometryType=\"XYZ\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << primalMesh.getNumberOfVertices() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "primal-mesh.h5:/vertexPos" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Geometry>" << std::endl;

	ss << "<Attribute Name=\"Face index\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << primalMesh.getNumberOfFaces() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << "primal-mesh.h5:/faceIndex" << std::endl;
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
	std::stringstream ss;
	ss << "<Grid Name=\"Dual Edges\">" << std::endl;
	ss << "<Topology TopologyType=\"Polyline\""
		<< " NumberOfElements=\"" << dualMesh.getNumberOfEdges() << "\"" 
		<< " NodesPerElement=\"2\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << 2 * dualMesh.getNumberOfEdges() << "\" DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << "dual-mesh.h5:/edge2vertex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Topology>" << std::endl;

	ss << "<Geometry GeometryType=\"XYZ\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfVertices() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "dual-mesh.h5:/vertexPos" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Geometry>" << std::endl;

	ss << "<Attribute Name=\"Edge index\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfEdges() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << "dual-mesh.h5:/edgeIdx" << std::endl;
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
	H5Reader h5reader;
	h5reader = H5Reader("dual-mesh.h5");
	const int nElements = h5reader.readDataSize("/face2vertex");
	assert(nElements > 0);

	std::stringstream ss;
	ss << "<Grid Name=\"Dual Faces\">" << std::endl;
	ss << "<Topology TopologyType=\"Mixed\""
		<< " NumberOfElements=\"" << dualMesh.getNumberOfFaces() << "\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << nElements << "\" DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << "dual-mesh.h5:/face2vertex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Topology>" << std::endl;

	ss << "<Geometry GeometryType=\"XYZ\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfVertices() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "dual-mesh.h5:/vertexPos" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Geometry>" << std::endl;

	ss << "<Attribute Name=\"Face index\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << "dual-mesh.h5:/faceIndex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;
	
	ss << "<Attribute Name=\"Face Fluid Type\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << "dualMeshTypes.h5:/faceFluidTypes" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;




	if (isWriteJfield) {
		ss << "<Attribute Name=\"Electric Current\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << " 3\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << "appm-" << iteration << ".h5:/J" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;
	}



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

	ss << "<Attribute Name=\"Current density consistent\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "appm-" << iteration << ".h5:/J_consistent" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "</Grid>" << std::endl;
	return ss.str();
}

const std::string AppmSolver::xdmf_GridDualCells(const int iteration) const
{
	H5Reader h5reader;
	h5reader = H5Reader("dual-mesh.h5");
	const int nElements = h5reader.readDataSize("/cell2vertex");
	assert(nElements > 0);

	std::stringstream ss;
	ss << "<Grid Name=\"Dual Cells\">" << std::endl;
	ss << "<Topology TopologyType=\"Mixed\""
		<< " NumberOfElements=\"" << dualMesh.getNumberOfCells() << "\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << nElements << "\" DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << "dual-mesh.h5:/cell2vertex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Topology>" << std::endl;

	ss << "<Geometry GeometryType=\"XYZ\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfVertices() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "dual-mesh.h5:/vertexPos" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Geometry>" << std::endl;

	ss << "<Attribute Name=\"Cell index\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfCells() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << "dual-mesh.h5:/cellIndex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"Cell Type\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfCells() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << "dualMeshTypes.h5:/cellFluidTypes" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	const std::string datafilename = (std::stringstream() << "appm-" << iteration << ".h5").str();

	ss << "<Attribute Name=\"Magnetic Flux Interpolated\" AttributeType=\"Vector\" Center=\"Node\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfVertices() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << datafilename << ":/Bvertex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"Electric Field Interpolated\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfCells() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << datafilename << ":/EdualCellCenter" << std::endl;
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
