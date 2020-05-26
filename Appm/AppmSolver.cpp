#include "AppmSolver.h"



AppmSolver::AppmSolver() 
	: AppmSolver(PrimalMesh::PrimalMeshParams())
{
}

AppmSolver::AppmSolver(const PrimalMesh::PrimalMeshParams & primalMeshParams)
{
	// Output file for time values at each timestep
	this->timeFile = std::ofstream("timesteps.dat");

	isStateWrittenToOutput = true;
	readParameters("AppmSolverParams.txt");
	init_meshes(primalMeshParams);  // Initialize primal and dual meshes
	if (primalMesh.getNumberOfCells() == 0) {
		return;
	}
	if (dualMesh.getNumberOfCells() == 0) {
		return;
	}

	if (dualMesh.getNumberOfCells() != primalMesh.getNumberOfVertices()) {
		std::cout << "Number of dual cells is not equal to primal vertices" << std::endl;
	}

	init_multiFluid("particleParameters.txt");
	
	switch (initType) {
	case 1:
	{
		std::cout << "Initialize fluid states: " << "Shock tube" << std::endl;
		const double zRef = 0.;
		init_SodShockTube(zRef);
	}
	break;

	case 2:
	{
		std::cout << "Initialize fluid states: " << "Uniform" << std::endl;
		double p = 1.0;
		double n = 1.0;
		Eigen::Vector3d uvec = Eigen::Vector3d::Zero();
		init_Uniformly(n, p, uvec);
	}
		break;

	case 3:
	{
		std::cout << "Initialize fluid states: " << "Explosion" << std::endl;
		const Eigen::Vector3d refPos = Eigen::Vector3d(0, 0, 0);
		const double radius = 0.2;
		init_Explosion(refPos, radius);
	}
		break;

	case 4:
	{
		init_testcase_frictionTerm();
	}
	break;

	default:
		std::cout << "InitType unknown: " << initType << std::endl;
		exit(-1);
	}

	const int nFluids = this->getNFluids();
	const int nFaces = dualMesh.getNumberOfFaces();
	const int nCells = dualMesh.getNumberOfCells();

	for (int i = 0; i < nCells; i++) {
		const Cell * cell = dualMesh.getCell(i);
		if (cell->getType() == Cell::Type::FLUID) {
			assert(fluidStates.col(i).allFinite());
		}
	}


	faceFluxes = Eigen::MatrixXd::Zero(5 * nFluids, nFaces);
	faceFluxesImExRusanov = Eigen::MatrixXd::Zero(nFaces, nFluids);

	LorentzForce_electric = Eigen::MatrixXd::Zero(3 * nFluids, nCells);
	LorentzForce_magnetic = Eigen::MatrixXd::Zero(3 * nFluids, nCells);

	init_maxwellStates();
	B_vertex = Eigen::Matrix3Xd::Zero(3, primalMesh.getNumberOfVertices());
	init_RaviartThomasInterpolation();


	if (appmParams.isMaxwellEnabled) {
		isWriteBfield = true;
		isWriteEfield = true;
	}

	
	//set_Efield_uniform(Eigen::Vector3d::UnitZ());
	//E_cc = getEfieldAtCellCenter();
	//setFluidFaceFluxes();
	//double dt = getNextFluidTimestepSize();
	//dt = 1;
	//double time = 5;
	////E_h.setZero();
	//Eigen::SparseMatrix<double> Msigma;
	//Msigma = get_Msigma_spd(J_h_aux, dt, time);
	//Eigen::sparseMatrixToFile(Msigma, "Msigma.dat");
	//std::cout << "nnz(Msigma): " << Msigma.nonZeros() << std::endl;
	//std::cout << "J_h_aux maxcoeff: " << J_h_aux.cwiseAbs().maxCoeff() << std::endl;
	////const int nEdges = primalMesh.getNumberOfEdges();
	//J_h = Msigma * E_h + J_h_aux;
	//Jcc = getCurrentDensityAtCellCenter();
	//J_h = J_h_aux;
	//assert(J_h.size() == nFaces);
	//assert(J_h_aux.size() == nFaces);

	//set_Bfield_azimuthal();
	//interpolateMagneticFluxToPrimalVertices();

	// write initial data to file (iteration = 0, time = 0)
	int iteration = 0;
	double time = 0;
	writeOutput(iteration, time);

	J_h_aux.setZero();
	J_h.setZero();
	//E_h.setZero();
	//E_cc = getEfieldAtCellCenter();
}


AppmSolver::~AppmSolver()
{
}

void AppmSolver::run()
{
	if (primalMesh.getNumberOfCells() == 0) {
		std::cout << "Primal mesh has no cells" << std::endl;
		return;
	}
	if (dualMesh.getNumberOfCells() == 0) {
		std::cout << "Dual mesh has no cells" << std::endl;
		return;
	}
	if (dualMesh.getNumberOfCells() != primalMesh.getNumberOfVertices()) {
		return;
	}

	double dt = appmParams.timestepSize;
	double dt_previous = dt;
	double time = 0;
	int iteration = 0;


	createStopFile();

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

	const int nFluids = this->getNFluids();
	const int nFaces = dualMesh.getNumberOfFaces();
	const int nCells = dualMesh.getNumberOfCells();

	//debug_checkCellStatus();
	

	auto timer_startAppmSolver = std::chrono::high_resolution_clock::now();

	// Set fluid face fluxes before start of iteration loop. 
	// This ensures that all required data of 'previous' timestep is also available at first iteration.
	setFluidFaceFluxes();


	/*
	* Time integration loop 
	*/
	while (iteration < maxIterations && time < maxTime && !isStopFileActive()) {
		std::cout << std::endl;
		std::cout << "*********************************************" << std::endl;
		std::cout << "* Iteration " << iteration << ",\t time = " << time << std::endl;
		std::cout << "*********************************************" << std::endl;
	
		dt_previous = dt;

		//fluidFluxes.setZero();
		//fluidSources.setZero();
		//faceFluxes.setZero();
		//faceFluxesImExRusanov.setZero();

		// Determine timestep
		if (appmParams.isFluidEnabled) {
			dt = getNextFluidTimestepSize();
			// Limit timestep such that maximum time is reached exactly
			if (time + dt >= maxTime) {
				std::cout << "timestep limited to reach maximum time; value before limiter is applied: dt = " << dt << std::endl;
				dt = maxTime - time;
			}
		}

		// Set explicit fluid source terms
		if (appmParams.isFluidEnabled) {
			if (isMagneticLorentzForceActive) {
				for (int idx = 0; idx < nCells; idx++) {
					const Cell * cell = dualMesh.getCell(idx);
					if (cell->getType() != Cell::Type::FLUID) {
						continue; // Skip cell that are not of type Fluid 
					}
					for (int fluidIdx = 0; fluidIdx < nFluids; fluidIdx++) {
						assert(idx < fluidStates.cols());
						assert(5 * fluidIdx + 3 < fluidStates.rows());
						const Eigen::Vector3d nu = fluidStates.col(idx).segment(5 * fluidIdx + 1, 3);
						assert(nu.allFinite());
						const Eigen::Vector3d B = B_vertex.col(idx);
						assert(B.allFinite());
						const int q = particleParams[fluidIdx].electricCharge;
						const double massRatio = particleParams[fluidIdx].mass;
						const Eigen::Vector3d result = q * 1. / massRatio * nu.cross(B);
						if (!result.allFinite()) {
							std::cout << "result: " << result.transpose() << std::endl;
							assert(result.allFinite());
						}
						LorentzForce_magnetic.col(idx).segment(3 * fluidIdx, 3) = result;
					}
				}
				std::cout << "maxCoeff F_L magnetic: " << LorentzForce_magnetic.cwiseAbs().maxCoeff() << std::endl;
			}
		}

		// Set explicit face fluxes
		if (appmParams.isFluidEnabled) {
			setFluidFaceFluxes();
		}


		// Affine-linear function for implicit-consistent formulation of current density J_h = Msigma * E_h + Jaux
		Eigen::SparseMatrix<double> Msigma;
		Msigma = get_Msigma_spd(J_h_aux, dt, time);

		// Maxwell equations
		if (appmParams.isMaxwellEnabled) {
			solveMaxwellSystem(time, dt, dt_previous, Msigma);
			// E- and B-field have been updated to new timestep
		}
		J_h = Msigma * E_h + J_h_aux;

		// Fluid equations
		if (appmParams.isFluidEnabled) {
			const int nFluids = this->getNFluids();
			const int nCells = dualMesh.getNumberOfCells();
			const int nFaces = dualMesh.getNumberOfFaces();

			// Set source term for electric Lorentz force (implicit)
			if (isElectricLorentzForceActive) {
				for (int idx = 0; idx < nCells; idx++) {
					const Cell * cell = dualMesh.getCell(idx);
					if (cell->getType() != Cell::Type::FLUID) {
						continue; // Skip cells that are not of type Fluid
					}
					for (int fluidIdx = 0; fluidIdx < nFluids; fluidIdx++) {
						const int q = particleParams[fluidIdx].electricCharge;
						const double n = fluidStates(5 * fluidIdx, idx);
						const double massRatio = particleParams[fluidIdx].mass;
						LorentzForce_electric.col(idx).segment(3 * fluidIdx, 3) = q * 1./massRatio * n * E_cc.col(idx);
					}
				}
			}
			std::cout << "Electric Lorentz Force maxCoeff: " << LorentzForce_electric.cwiseAbs().maxCoeff() << std::endl;

			setFluidSourceTerm();

			// Set implicit terms for mass fluxes
			if (appmParams.isMassFluxSchemeImplicit) {
				setImplicitMassFluxTerms(dt);
			}

			setSumOfFaceFluxes();
			updateFluidStates(dt);

			//debug_checkCellStatus();
		}
		std::cout << "dt = " << dt << std::endl;
		iteration++;
		time += dt;
		writeOutput(iteration, time);
	}
	std::cout << std::endl;
	std::cout << "Final time:      " << time << std::endl;
	std::cout << "Final iteration: " << iteration << std::endl;

	auto timer_endAppmSolver = std::chrono::high_resolution_clock::now();
	auto delta_appmSolver = std::chrono::duration<double>(timer_endAppmSolver - timer_startAppmSolver);
	std::cout << "Elapsed time for APPM solver: " << delta_appmSolver.count() << std::endl;

	std::cout << printSolverParameters() << std::endl;

	// test_raviartThomas();

	// Use Paraview (version 5.6.0) to visualize.
	// Light data is given in XDMF files (version 3). 
	// Heavy data is stored in HDF5 files.

	// Note that polyhedrons have been added to XDMF in version 3; version 2 does not support polyhedrons.
	// Moreover, Paraview has 3 readers: XDMF Reader, Xdmf3ReaderS, Xdmf3ReaderT.
	// Xdmf3Reader3: can     read polyhedrons, but not grid-of-grids (i.e., grids with GridType=Tree).
	// XDMF Reader:  can not read polyhedrons, but     grid-of-grids.
	// Therefore, the mesh data for vertices, edges, and faces, are separated from the volume data.
	writeXdmf("appm.xdmf");
	writeXdmfDualVolume("appm-volume.xdmf");
	//writeXdmfDualFaceFluxes();
}

void AppmSolver::debug_checkCellStatus() const
{
	std::cout << "Check cell status" << std::endl;
	const int cellIdx = 52;
	const Cell * cell = dualMesh.getCell(cellIdx);
	auto cellFaces = cell->getFaceList();
	std::cout << "cell state: (id = " << cell->getIndex() << ")" << std::endl;
	Eigen::VectorXd state = fluidStates.col(cell->getIndex());
	for (int fluidIdx = 0; fluidIdx < getNFluids(); fluidIdx++) {
		std::cout << std::scientific << state.segment(5 * fluidIdx, 5).transpose() << std::endl;
	}
	
	Eigen::VectorXd sumFluxes(faceFluxes.rows());
	sumFluxes.setZero();
	for (auto face : cellFaces) {
		const int s_ki = cell->getOrientation(face);
		const double fA = face->getArea();
		Eigen::VectorXd flux = faceFluxes.col(face->getIndex());
		sumFluxes += s_ki * flux * fA;
		std::cout << "face " << face->getIndex() << ", orientation = " << s_ki << std::endl;
		std::cout << "area = " << fA << ", normal = " << face->getNormal().transpose() << std::endl;
		for (int fluidIdx = 0; fluidIdx < getNFluids(); fluidIdx++) {
			std::cout << std::scientific << flux.segment(5 * fluidIdx, 5).transpose() << std::endl;
		}
	}
	std::cout << "sum of fluxes: \t" << std::endl;
	for (int fluidIdx = 0; fluidIdx < getNFluids(); fluidIdx++) {
		std::cout << sumFluxes.segment(5 * fluidIdx, 5).transpose() << std::endl;
	}
	std::cout << std::endl;

}

const int AppmSolver::getNFluids() const
{
	return particleParams.size();
}

std::string AppmSolver::printSolverParameters() const
{
	std::stringstream ss;
	ss << "Appm Solver parameters:" << std::endl;
	ss << "=======================" << std::endl;
	ss << "maxIterations:  " << maxIterations << std::endl;
	ss << "maxTime:        " << maxTime << std::endl;
	ss << "timestepSize:   " << appmParams.timestepSize << std::endl;
	ss << "isFluidEnabled: " << appmParams.isFluidEnabled << std::endl;
	ss << "isMaxwellEnabled: " << appmParams.isMaxwellEnabled << std::endl;
	ss << "lambdaSquare: " << lambdaSquare << std::endl;
	ss << "massFluxScheme: " << appmParams.isMassFluxSchemeImplicit << std::endl;
	ss << "initType: " << initType << std::endl;
	ss << "isElectricLorentzForceActive: " << isElectricLorentzForceActive << std::endl;
	ss << "isMagneticLorentzForceActive: " << isMagneticLorentzForceActive << std::endl;
	ss << "isShowDataWriterOutput:" << isShowDataWriterOutput << std::endl;
	ss << "isMaxwellCurrentDefined: " << appmParams.isMaxwellCurrentDefined << std::endl;
	ss << "isEulerMaxwellCouplingEnabled: " << appmParams.isEulerMaxwellCouplingEnabled << std::endl;
	ss << "=======================" << std::endl;
	return ss.str();
}

void AppmSolver::init_maxwellStates()
{
	std::cout << "Initialize Maxwell states" << std::endl;


	E_h = Eigen::VectorXd::Zero(primalMesh.getNumberOfEdges());
	B_h = Eigen::VectorXd::Zero(primalMesh.getNumberOfFaces());

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
	M_perot = initPerotInterpolationMatrix();
	
	// Electric field at cell centers
	E_cc = Eigen::Matrix3Xd::Zero(3, dualMesh.getNumberOfCells());

	// Current density at cell centers
	Jcc = Eigen::Matrix3Xd::Zero(3, dualMesh.getNumberOfCells());

	// Initialize matrix that maps electric field on primal edges to electric current on dual faces; see Lorentz force in Fluid equations
	//initMsigma();


	const int nFaces = dualMesh.getNumberOfFaces();
	J_h = Eigen::VectorXd::Zero(nFaces);
	J_h_previous = J_h;

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
	fluidSources = Eigen::MatrixXd::Zero(fluidStateLength, nCells);
	sumOfFaceFluxes = Eigen::MatrixXd::Zero(fluidStates.rows(), nCells);


	const int nFaces = dualMesh.getNumberOfFaces();
	faceTypeFluids = Eigen::MatrixXi::Zero(nFaces, nFluids);
	for (int fluidIdx = 0; fluidIdx < nFluids; fluidIdx++) {
		for (int i = 0; i < nFaces; i++) {
			const Face * face = dualMesh.getFace(i);
			Face::Type faceType = getFaceTypeOfFluid(face, fluidIdx);
			faceTypeFluids(i, fluidIdx) = static_cast<int>(faceType);
		}
	}
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
		const Eigen::Vector3d uL = Eigen::Vector3d::Zero();

		const double pR = 0.1;
		const double nR = 0.125;
		const Eigen::Vector3d uR = Eigen::Vector3d::Zero();

		Eigen::VectorXd singleFluidStateLeft;
		singleFluidStateLeft = Physics::primitive2state(epsilon2, nL, pL, uL);
		assert(singleFluidStateLeft.size() == 5);
		//singleFluidStateLeft.setZero();
		//singleFluidStateLeft(0) = nL;
		//singleFluidStateLeft.segment(1, 3) = nL * uL;
		//singleFluidStateLeft(4) = pL / (Physics::gamma - 1) + 0.5 * nL * uL.squaredNorm();

		Eigen::VectorXd singleFluidStateRight;
		singleFluidStateRight = Physics::primitive2state(epsilon2, nR, pR, uR);
		assert(singleFluidStateRight.size() == 5);
		//singleFluidStateRight.setZero();
		//singleFluidStateRight(0) = nR;
		//singleFluidStateRight.segment(1, 3) = nR * uR;
		//singleFluidStateRight(4) = pR / (Physics::gamma - 1) + 0.5 * nR * uR.squaredNorm();

		leftState.segment(5 * k, 5) = singleFluidStateLeft;
		rightState.segment(5 * k, 5) = singleFluidStateRight;
	}

	for (int i = 0; i < nCells; i++) {
		const Cell * cell = dualMesh.getCell(i);
		Eigen::VectorXd cellState(fluidStateLength);
		cellState.setConstant(std::nan("")); // NaN

		if (cell->getType() == Cell::Type::FLUID) {
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

void AppmSolver::init_Uniformly(const double n, const double p, const Eigen::Vector3d uvec)
{
	const int nCells = dualMesh.getNumberOfCells();
	const int nFluids = getNFluids();
	const int fluidStateLength = getFluidStateLength();
	Eigen::VectorXd state(fluidStateLength);
	for (int k = 0; k < nFluids; k++) {
		const double epsilon2 = particleParams[k].mass;
		Eigen::VectorXd singleFluidState;
		singleFluidState = Physics::primitive2state(epsilon2, n, p, uvec);
		//const double e = p / ((Physics::gamma - 1) * epsilon2 * n);
		//const double etot = e + 0.5 * pow(u, 2);
		//singleFluidState(0) = n;
		//singleFluidState.segment(1, 3) = n * u * Eigen::Vector3d::UnitZ();
		//singleFluidState(4) = n * etot;
		state.segment(5 * k, 5) = singleFluidState;
	}
	fluidStates.setConstant(std::nan(""));
	for (int i = 0; i < nCells; i++) {
		const Cell * cell = dualMesh.getCell(i);
		if (cell->getType() == Cell::Type::FLUID) {
			fluidStates.col(i) = state;
		}
	}
}

void AppmSolver::init_Explosion(const Eigen::Vector3d refPos, const double radius)
{

	const int nCells = dualMesh.getNumberOfCells();
	const int nFluids = getNFluids();
	for (int fluidIdx = 0; fluidIdx < nFluids; fluidIdx++) {
		const double massRatio = particleParams[fluidIdx].mass;
		Eigen::VectorXd state_lo = Physics::primitive2state(massRatio, 1, 1, Eigen::Vector3d::Zero());
		Eigen::VectorXd state_hi = Physics::primitive2state(massRatio, 1, 5, Eigen::Vector3d::Zero());
		
		for (int cIdx = 0; cIdx < nCells; cIdx++) {
			const Cell * cell = dualMesh.getCell(cIdx);
			if (cell->getType() != Cell::Type::FLUID) { continue; } // Skip cells that are not Fluid type
			const Eigen::Vector3d cc = cell->getCenter();
			const bool isInside = (cc - refPos).norm() < radius;
			const Eigen::VectorXd state = isInside ? state_hi : state_lo;
			fluidStates.col(cIdx).segment(5 * fluidIdx, 5) = state;
		}
	}
}

void AppmSolver::init_testcase_frictionTerm()
{
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
Eigen::SparseMatrix<double> AppmSolver::initPerotInterpolationMatrix()
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
	Eigen::SparseMatrix<double> M;
	M = Eigen::SparseMatrix<double>(3 * nCells, nEdges);
	M.setFromTriplets(triplets.begin(), triplets.end());
	M.makeCompressed();
	
	std::cout << text << ": DONE" << std::endl;
	return M;
	// Eigen::sparseMatrixToFile(M_perot, "Mperot.dat");
}

const Eigen::Matrix3Xd AppmSolver::getEfieldAtCellCenter()
{
	const int nCells = dualMesh.getNumberOfCells();
	if (M_perot.size() == 0 || M_perot.nonZeros() == 0) {
		std::cout << "Initialize Perot interpolation matrix for Electric Field (edges to cell centers)" << std::endl;
		M_perot = initPerotInterpolationMatrix();
	}
	assert(M_perot.size() > 0);
	assert(M_perot.nonZeros() > 0);

	Eigen::VectorXd E_cc_vectorFormat = M_perot * E_h;

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
	Eigen::Map<Eigen::Matrix3Xd> result(E_cc_vectorFormat.data(), 3, nCells);
	return result;
}

/**
* Interpolate finite face current to current density at cell center. 
* This is especially useful for visualization.
*
* @return current density vector at dual cell centers. 
*/
const Eigen::Matrix3Xd AppmSolver::getCurrentDensityAtCellCenter()
{
	const int nCells = dualMesh.getNumberOfCells();
	const int nFaces = dualMesh.getNumberOfFaces();
	Eigen::SparseMatrix<double> M;

	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;

	for (int cIdx = 0; cIdx < nCells; cIdx++) {
		const Cell * cell = dualMesh.getCell(cIdx);
		if (cell->getType() != Cell::Type::FLUID) { continue; }
		auto cellFaces = cell->getFaceList();
		for (auto face : cellFaces) {
			const int faceIdx = face->getIndex();
			const Eigen::Vector3d cc = cell->getCenter();
			const Eigen::Vector3d fc = face->getCenter(); 
			const Eigen::Vector3d r = fc - cc;
			const double Vk = cell->getVolume();

			const int incidence = getOrientation(cell, face);

			const double value = 1./Vk * incidence;

			for (int i = 0; i < 3; i++) {
				const int row = 3 * cIdx + i;
				const int col = faceIdx;
				triplets.push_back(T(row, col, value * r(i)));
			}
		}
	}
	M = Eigen::SparseMatrix<double>(3 * nCells, nFaces);
	M.setFromTriplets(triplets.begin(), triplets.end());
	M.makeCompressed();

	Eigen::VectorXd Jcc_vectorFormat = M * J_h;
	Eigen::Map<Eigen::Matrix3Xd> result(Jcc_vectorFormat.data(), 3, nCells);
	return result;
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
	if (!isMomentumFluxActive) {
		std::cout << "Warning: momentum flux disabled in implicit-consistent current model" << std::endl;
	}

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
				if (cell->getType() != Cell::Type::FLUID) {
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
						const bool isFaceTypeOpening = dualFace->getType() == Face::Type::OPENING;
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

	//assert(Msigma.nonZeros() > 0);
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
			if (cell->getType() != Cell::Type::FLUID) {
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

/**
* @return Fluid state in given cell and fluid, projected in direction of a unit normal vector. 
* The normal vector may point in any direction. 
*/
const Eigen::Vector3d AppmSolver::getFluidState(const int cellIdx, const int fluidIdx, const Eigen::Vector3d & faceNormal) const
{
	assert(cellIdx >= 0);
	assert(cellIdx < fluidStates.cols());
	assert(fluidIdx >= 0);
	assert(fluidIdx < this->getNFluids());
	const double tol = 4 * std::numeric_limits<double>::epsilon();
	assert(abs(faceNormal.norm() - 1) <= tol); // normal vector should be of unit length
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

/**
* Get the face flux with explicit Rusanov scheme. 
*/
const Eigen::Vector3d AppmSolver::getRusanovFluxExplicit(const int faceIdx, const int fluidIdx) const
{
	assert(false); // Do not use this anymore

	const bool showOutput = false;
	Eigen::Vector3d qL, qR;
	Eigen::Vector3d faceFlux;
	faceFlux.setZero();

	const Face * face = dualMesh.getFace(faceIdx);
	const Face::Type faceFluidType = getFaceTypeOfFluid(face, fluidIdx);

	switch (faceFluidType) {
	case Face::Type::INTERIOR:
		getAdjacientCellStates(face, fluidIdx, qL, qR);
		faceFlux = Physics::getRusanovFlux(qL, qR, showOutput);
		break;


	}

	return faceFlux; 
}

const Eigen::Vector3d AppmSolver::getRusanovFluxImEx(const int faceIdx, const int fluidIdx, const double dt) 
{
	const Face * face = dualMesh.getFace(faceIdx);
	const Eigen::Vector3d faceNormal = face->getNormal();
	//const Face::FluidType faceFluidType = face->getFluidType();
	const Face::Type faceFluidType = getFaceTypeOfFluid(face, fluidIdx);

	Eigen::Vector3d qL, qR;
	const std::pair<int, int> adjacientCellIdx = getAdjacientCellStates(face, fluidIdx, qL, qR);

	double implicitExtraTerms = 0;
	int cellIndex;

	if (faceFluidType != Face::Type::INTERIOR) {
		assert(adjacientCellIdx.first == -1 || adjacientCellIdx.second == -1);
	}
	const bool isInteriorFace = faceFluidType == Face::Type::INTERIOR;
	const bool isBoundaryFace = !isInteriorFace;
	if (appmParams.isMassFluxSchemeImplicit) {
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

const std::pair<int,int> AppmSolver::getAdjacientCellStates(const Face * face, const int fluidIdx, Eigen::Vector3d & qL, Eigen::Vector3d & qR) const
{
	assert(face != nullptr);
	const std::vector<Cell*> faceCells = face->getCellList();
	assert(faceCells.size() >= 1);

	const Eigen::Vector3d faceNormal = face->getNormal();
	const int orientation = getOrientation(faceCells[0], face);
	qL.setZero();
	qR.setZero();
	int idxL = -1;
	int idxR = -1;

	//const Face::FluidType faceFluidType = face->getFluidType();
	const Face::Type faceFluidType = getFaceTypeOfFluid(face, fluidIdx);
	switch (faceFluidType) {
	case Face::Type::INTERIOR:
		assert(faceCells.size() == 2);
		assert(faceCells[0]->getType() == Cell::Type::FLUID);
		assert(faceCells[1]->getType() == Cell::Type::FLUID);
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

	case Face::Type::OPENING:
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

	case Face::Type::WALL:
		assert(faceCells.size() >= 1 && faceCells.size() <= 2);
		if (faceCells[0]->getType() == Cell::Type::FLUID) {
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
			assert(faceCells[1]->getType() == Cell::Type::FLUID);
			assert(false); // Not yet implemented
		}
		break;

	case Face::Type::TERMINAL:
		std::cout << "FaceFluidType TERMINAL should not be used directly; call getFaceTypeOfFluid(face,fluidIdx) instead" << std::endl;
		assert(false); // Not yet implemented
		break;

	default:
		std::cout << "Face Fluid Type not implemented: " << faceFluidType << std::endl;
		assert(false);
	}

	return std::pair<int, int>(idxL, idxR);
}

//const double AppmSolver::getMomentumUpdate(const int k, const Eigen::Vector3d & nvec, const int fluidIdx) const
//{
//	double result = 0;
//	const Cell * cell = dualMesh.getCell(k);
//	const std::vector<Face*> cellFaces = cell->getFaceList();
//	for (auto face : cellFaces) {
//		const double faceArea = face->getArea();
//		const Eigen::Vector3d faceNormal = face->getNormal();
//		const Eigen::Vector3d faceFlux = getRusanovFluxExplicit(face->getIndex(), fluidIdx);
//		const double momentumFlux = faceFlux(1);
//		double localResult = momentumFlux * faceArea * faceNormal.dot(nvec);
//		result += localResult;
//	}
//	result /= cell->getVolume();
//	return result;
//}

/**
* Create stop file with default value.
*/
void AppmSolver::createStopFile()
{
	std::ofstream(stopFilename) << 0 << std::endl;
}

/**
* @return if value in stop file is positive; this indicates to stop the solver iteration loop. 
*/
bool AppmSolver::isStopFileActive()
{
	int stopValue = 0;
	std::ifstream(stopFilename) >> stopValue;
	std::cout << "Check for value in stop file: " << stopValue << std::endl;
	if (stopValue > 0) {
		std::cout << "!!! Stop !!!" << std::endl;
		return true;
	}
	return false;
}

/**
* Set fluid face fluxes at each dual face.
*/
void AppmSolver::setFluidFaceFluxes()
{
	const int nFaces = dualMesh.getNumberOfFaces();
	const int nFluids = getNFluids();

	for (int fidx = 0; fidx < nFaces; fidx++) {
		const Face * face = dualMesh.getFace(fidx);
		const Eigen::Vector3d faceNormal = face->getNormal();

		// skip faces that have no adjacient fluid cell
		if (!face->hasFluidCells()) { continue; }

		assert(face->getType() != Face::Type::DEFAULT);

		for (int fluidIdx = 0; fluidIdx < nFluids; fluidIdx++) {
			Eigen::VectorXd faceFlux3d(5);
			Eigen::Vector3d flux = getSpeciesFaceFlux(face, fluidIdx);
			faceFlux3d(0) = flux(0);
			faceFlux3d.segment(1, 3) = flux(1) * faceNormal;
			faceFlux3d(4) = flux(2);
			faceFluxes.col(fidx).segment(5 * fluidIdx, 5) = faceFlux3d;
		}
	}
}

/**
* Collect sum of face fluxes for each cell (i.e., summation of the divergence term)
*/
void AppmSolver::setSumOfFaceFluxes() {
	const int nCells = dualMesh.getNumberOfCells();
	sumOfFaceFluxes.setZero();
	for (int k = 0; k < nCells; k++) {
		const Cell * cell = dualMesh.getCell(k);
		// Skip non-fluid cells
		if (cell->getType() != Cell::Type::FLUID) { continue; }

		double cellVolume = cell->getVolume();

		// get sum of face fluxes
		const std::vector<Face*> cellFaces = cell->getFaceList();
		for (auto face : cellFaces) {
			const double faceArea = face->getArea();
			const int orientation = getOrientation(cell, face);
			sumOfFaceFluxes.col(k) += orientation * faceFluxes.col(face->getIndex()) * faceArea / cellVolume;
		}
	}
}

/**
* Given the explicit Rusanov fluxes at each fluid face, this function 
* adds the extra terms in mass flux that stem from the implicit mass flux formulation.
*/
void AppmSolver::setImplicitMassFluxTerms(const double dt)
{
	assert(dt > 0);
	const int nFaces = dualMesh.getNumberOfFaces();
	const int nFluids = getNFluids();

	// For each face ...
	for (int faceIdx = 0; faceIdx < nFaces; faceIdx++) {
		const Face * face = dualMesh.getFace(faceIdx);

		// skip faces that have no adjacient fluid cell
		if (!face->hasFluidCells()) { continue; }
		assert(face->getType() != Face::Type::DEFAULT);

		const Eigen::Vector3d ni = face->getNormal();
		const Face::Type faceType = face->getType();
		const double numSchemeFactor = (faceType == Face::Type::INTERIOR) ? 0.5 : 1;

		// ... for each fluid ...
		for (int fluidIdx = 0; fluidIdx < nFluids; fluidIdx++) {
			// ... for each adjacient cell of that face ...
			auto adjacientCells = face->getCellList();
			double cellSum = 0;
			for (auto cell : adjacientCells) {
				// Skip non-fluid cells
				if (cell->getType() != Cell::Type::FLUID) {
					continue;
				}
				const double Vk = cell->getVolume();

				// Implicit terms due to face fluxes
				auto cellFaces = cell->getFaceList();
				double faceSum = 0;
				for (auto cellFace : cellFaces) {
					const int j = cellFace->getIndex();
					const double Aj = cellFace->getArea();
					const int s_kj = cell->getOrientation(cellFace);
					const Eigen::Vector3d nj = cellFace->getNormal();
					double ni_dot_nj = ni.dot(nj);
					if (abs(ni_dot_nj) < 1e-10) {
						ni_dot_nj = 0; // truncate small values
					}
					const double fj = faceFluxes.col(j).segment(5 * fluidIdx + 1, 3).dot(nj);
					faceSum += s_kj * fj * Aj * ni_dot_nj;
				}
				cellSum -= numSchemeFactor / Vk * faceSum;

				// Implicit term due to fluid momentum source 				
				Eigen::Vector3d Sk; // Momentum source term for this fluid
				Sk = fluidSources.col(cell->getIndex()).segment(5 * fluidIdx + 1, 3);
				cellSum += Sk.dot(ni);
			}
			cellSum *= dt;

			// Update mass flux with implicit terms
			faceFluxes(5 * fluidIdx, faceIdx) += cellSum;
		}
	}
}

Eigen::Vector3d AppmSolver::getSpeciesFaceFlux(const Face * face, const int fluidIdx)
{
	assert(face != nullptr);
	assert(fluidIdx >= 0);
	assert(fluidIdx < getNFluids());
	const int faceIdx = face->getIndex();

	Eigen::Vector3d flux;
	flux.setZero();

	Face::Type faceFluidType = face->getType();
	const bool isElectronFluid = particleParams[fluidIdx].electricCharge < 0;
	const bool isCathode = faceFluidType == Face::Type::TERMINAL && face->getCenter()(2) < 0; // TODO

	if (false && isCathode && isElectronFluid) {
		flux = getSpeciesFaceFluxAtCathode(face, fluidIdx);
	}
	else {
		Eigen::Vector3d qL, qR;
		std::pair<int,int> cellIndices = getAdjacientCellStates(face, fluidIdx, qL, qR);
		flux = Physics::getRusanovFlux(qL, qR);
	}
	return flux;
}

const Eigen::Vector3d AppmSolver::getSpeciesFaceFluxAtCathode(const Face * face, const int fluidIdx)
{
	const double Ts = 3695; // melting point of Wolfram, Kelvin
	const double workFunction = 4.55; // units: eV
	const double j_em = Physics::thermionicEmissionCurrentDensity(Ts, workFunction); 

	// Check if face is at cathode terminal and we have an electron fluid
	const bool isElectronFluid = particleParams[fluidIdx].electricCharge < 0;
	assert(isElectronFluid);
	Face::Type faceFluidType = face->getType();
	assert(faceFluidType == Face::Type::TERMINAL);

	// Set a default value for testing
	Eigen::Vector3d flux;
	flux(0) = 1;
	flux(1) = 2;
	flux(2) = 1;

	flux(1) *= face->getNormal().dot(Eigen::Vector3d::UnitZ());
	flux *= face->getNormal().dot(Eigen::Vector3d::UnitZ());
	return flux;
}

void AppmSolver::setFluidSourceTerm()
{
	const int nCells = dualMesh.getNumberOfCells();
	const int nFluids = getNFluids();

	// Set source term
	for (int cIdx = 0; cIdx < nCells; cIdx++) {
		const Cell * cell = dualMesh.getCell(cIdx);

		// Skip cell that ore not of type Fluid
		if (cell->getType() != Cell::Type::FLUID) {
			continue;
		}
		// For all fluids ... 
		for (int fluidIdx = 0; fluidIdx < nFluids; fluidIdx++) {
			Eigen::VectorXd srcLocal(5);
			srcLocal.setZero();

			//srcLocal = testcase_001_FluidSourceTerm(time, cell, fluidIdx);

			if (isElectricLorentzForceActive) {
				srcLocal.segment(1, 3) += LorentzForce_electric.col(cIdx).segment(3 * fluidIdx, 3);
			}
			if (isMagneticLorentzForceActive) {
				srcLocal.segment(1, 3) += LorentzForce_magnetic.col(cIdx).segment(3 * fluidIdx, 3);
			}

			// TODO: reaction source terms

			fluidSources.col(cIdx).segment(5 * fluidIdx, 5) = srcLocal;
		}
	}
	std::cout << "fluid sources maxCoeff: " << fluidSources.cwiseAbs().maxCoeff() << std::endl;
}

/** 
* Update to next timestep: U(m+1) = U(m) - dt / volume * sum(fluxes) + dt * source.
*/
void AppmSolver::updateFluidStates(const double dt)
{
	const int nCells = dualMesh.getNumberOfCells();
	for (int k = 0; k < nCells; k++) {
		fluidStates.col(k) += -dt * sumOfFaceFluxes.col(k) + dt * fluidSources.col(k);
	}
}

void AppmSolver::solveMaxwellSystem(const double time, const double dt, const double dt_previous, const Eigen::SparseMatrix<double> & Msigma)
{
	// Number of primal edges
	const int nEdges = primalMesh.getNumberOfEdges();

	// Number of primal vertices on terminals
	const int nPrimalTerminalVertices = primalMesh.getMeshInfo().nVerticesTerminal;
	assert(nPrimalTerminalVertices > 0);


	// Setup of implicit system of equations for the reformulated Ampere equation, 
	//     d_tt (Meps * e) + C' * Mnu * C * e = -d_t (j), 
	// with j = Msigma * e  (+ j_aux)  

	// New state vector
	Eigen::VectorXd x(maxwellState.size());

	// Get matrix such that J_h = Msigma * E_h + J_h_aux
	//Eigen::SparseMatrix<double> Msigma;
	//Msigma = get_Msigma_spd(J_h_aux, dt);

	// Get interior matrix 
	Eigen::SparseMatrix<double> Msigma_inner;
	Msigma_inner = Msigma.topLeftCorner(nEdges, nEdges);

	// System matrix on left hand side, i.e., M*x = rhs, to be solved for x
	Eigen::SparseMatrix<double> M;
	assert(M1.size() > 0);
	assert(M1.nonZeros() > 0);
	assert(M2.size() > 0);
	assert(M2.nonZeros() > 0);

	M = M1 + pow(dt, 2) * M2;
	M += dt * Q.transpose() * Msigma_inner * Q;
	M.makeCompressed();
	//Eigen::sparseMatrixToFile(M, "M.dat");

	double dt_ratio = dt / dt_previous;

	// Data vector for right hand side
	Eigen::VectorXd rhs(x.size());
	rhs.setZero();
	rhs = M1 * (1 + dt_ratio) * maxwellState;
	rhs -= dt_ratio * M1 * maxwellStatePrevious;
	//rhs += dt * Q.transpose() * Msigma_inner * Q * maxwellState;
	rhs -= dt * Q.transpose() * (J_h_aux - J_h).topRows(nEdges);


	// The system of equations has fixed and free values; 
	// - fixed values: electric potential at terminals (Dirichlet boundary condition)
	// - free  values: electric potential at non-terminal vertices, and electric voltages at non-boundary edges
	int nDirichlet = nPrimalTerminalVertices;
	if (appmParams.isMaxwellCurrentDefined) {
		// If we define a current, then only one electrode is grounded 
		// and the other is on a free-floating potential.
		//assert(nPrimalTerminalVertices % 2 == 0);
		nDirichlet = nPrimalTerminalVertices;
	} else {
		// If no current is defined, then we specify the electric potential 
		// on both electrodes as fixed values.
		nDirichlet = nPrimalTerminalVertices;
	}
	int nFree = maxwellState.size() - nDirichlet;

	// The vector of degrees of freedom (DoF) is sorted such that free values are in front of fixed values:
	// x = [freeValues, fixedValues]
	// Therefore, the system of free DoF is given as: (without considering proper array sizes)
	// M_free * x_free = -M_fixed * x_fixed + rhs

	// Dirichlet conditions
	Eigen::SparseMatrix<double> Md = M.rightCols(nDirichlet);

	// Data vector for degrees of freedom
	Eigen::VectorXd xf(nFree); // Vector of free coefficients
	Eigen::VectorXd xd = setVoltageBoundaryConditions(time);
	assert(Md.cols() == xd.rows());

	// Subtract Dirichlet values from left side
	rhs -= Md * xd;

	// Load vector for free coefficients
	Eigen::VectorXd rhsFree = rhs.topRows(nFree);

	// Matrix of free coefficients
	Eigen::SparseMatrix<double> Mf = M.topLeftCorner(nFree, nFree);
	Mf.makeCompressed();

	// Solve system
	std::cout << "Setup Maxwell solver" << std::endl;

	// Note: the matrix Mf is not symmetric, neither positive definite!
	//xf = solveMaxwell_sparseLU(Mf, rhsFree);
	//xf = solveMaxwell_BiCGStab(Mf, rhsFree); // Slow solver
	//xf = solveMaxwell_LSCG(Mf, rhsFree); // error: iterative solver base is not initialized
	//xf = solveMaxwell_CG(Mf, rhsFree);
	xf = solveMaxwell_PardisoLU(Mf, rhsFree);

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

	std::cout << std::endl;
	std::cout << "Init dual mesh" << std::endl;
	dualMesh = DualMesh();
	dualMesh.init_dualMesh(primalMesh, primalParams.getElectrodeRadius());
	dualMesh.writeToFile();
	dualMesh.writeXdmf();

	std::cout << "Dual mesh has " << dualMesh.getNumberOfVertices() << " vertices" << std::endl;
}

/**
* Write XDMF output file for edges and faces of primal and dual meshes.
*/
void AppmSolver::writeXdmf(const std::string & filename)
{
	std::cout << "Write XDMF output file: " << filename << std::endl;
	
	// Check that filename has correct extension
	std::size_t pos = filename.rfind(".xdmf");
	assert(pos != std::string::npos);
	assert(pos > 0);

	const int nTimesteps = timeStamps.size();
	
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

/**
* Write XDMF output file for cells of dual mesh. 
*/
void AppmSolver::writeXdmfDualVolume(const std::string & filename)
{
	std::cout << "Write XDMF output file: " << filename << std::endl;

	// Check that filename has correct extension
	std::size_t pos = filename.rfind(".xdmf");
	assert(pos != std::string::npos);
	assert(pos > 0);

	const int nTimesteps = timeStamps.size();

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

/**
* Write data to output file.
*/
void AppmSolver::writeOutput(const int iteration, const double time)
{
	timeStamps.push_back(time);
	timeFile << iteration << "," << time << std::endl;

	std::cout << "Write output at iteration " << iteration << ", time = " << time << std::endl;
	const std::string filename = (std::stringstream() << "appm-" << iteration << ".h5").str();

	const int nPrimalFaces = primalMesh.getNumberOfFaces();
	const int nDualEdges = dualMesh.getNumberOfEdges();
	const int nDualFaces = dualMesh.getNumberOfFaces();
	const int nDualCells = dualMesh.getNumberOfCells();

	H5Writer h5writer(filename);
	h5writer.setShowOutput(this->isShowDataWriterOutput);

	// Fluid states
	writeFluidStates(h5writer);

	// Maxwell states
	writeMaxwellStates(h5writer);

	Eigen::VectorXd timeVec(1);
	timeVec(0) = time;
	h5writer.writeData(timeVec, "/time");
	Eigen::VectorXi iterVec(1);
	iterVec(0) = iteration;
	h5writer.writeData(iterVec, "/iteration");
}

void AppmSolver::writeFluidStates(H5Writer & writer)
{
	const int nFluids = this->getNFluids();
	const int nCells = dualMesh.getNumberOfCells();
	const int nFaces = dualMesh.getNumberOfFaces();

	for (int fluidIdx = 0; fluidIdx < nFluids; fluidIdx++) {
		const std::string fluidName = particleParams[fluidIdx].name;
		const std::string fluidTag = (std::stringstream() << "/" << fluidName << "-").str();
		const std::string pressureTag = fluidTag + "pressure";
		const std::string velocityTag = fluidTag + "velocity";
		const std::string densityTag = fluidTag + "density";
		const std::string numberDensityTag = fluidTag + "numberDensity";
		Eigen::VectorXd numberDensity = fluidStates.row(5 * fluidIdx);
		Eigen::VectorXd density(nCells);
		Eigen::MatrixXd velocity(3, nCells);
		Eigen::VectorXd pressure(nCells);
		Eigen::VectorXd sumMassFluxes = sumOfFaceFluxes.row(5 * fluidIdx);
		Eigen::MatrixXd sumMomentumFluxes = sumOfFaceFluxes.block(5 * fluidIdx + 1, 0, 3, nCells);
		Eigen::VectorXd sumEnergyFluxes = sumOfFaceFluxes.row(5*fluidIdx + 4);

		const std::string stateN = fluidTag + "stateN";
		const std::string stateU = fluidTag + "stateU";
		const std::string stateE = fluidTag + "stateE";
		Eigen::VectorXd qN(nCells);
		Eigen::MatrixXd qU(3, nCells);
		Eigen::VectorXd qE(nCells);

		// Initialize data vectors with NaN (not-a-number)
		density.setConstant(std::nan("")); 
		velocity.setConstant(std::nan(""));
		pressure.setConstant(std::nan(""));
		qN.setConstant(std::nan("")); 
		qU.setConstant(std::nan("")); 
		qE.setConstant(std::nan("")); 

		const double epsilon2 = particleParams[fluidIdx].mass;
		for (int i = 0; i < nCells; i++) {
			const Cell * cell = dualMesh.getCell(i);
			if (cell->getType() != Cell::Type::FLUID) { 
				continue; // Skip non-fluid cells
			}
			const Eigen::VectorXd state = fluidStates.block(5 * fluidIdx, i, 5, 1);
			if (isStateWrittenToOutput) {
				qN(i) = state(0);
				qU.col(i) = state.segment(1, 3);
				qE(i) = state(4);
			}
			double n = 0;
			double p = 0;
			Eigen::Vector3d u;
			u.setZero();
			Physics::state2primitive(epsilon2, state, n, p, u);
			density(i) = epsilon2 * n;
			velocity.col(i) = u;
			pressure(i) = p;
		}

		writer.writeData(density, densityTag);
		writer.writeData(numberDensity, numberDensityTag);
		writer.writeData(pressure, pressureTag);
		writer.writeData(velocity, velocityTag);
		writer.writeData(sumMassFluxes, (std::stringstream() << fluidTag << "sumMassFlux").str());
		writer.writeData(sumMomentumFluxes, (std::stringstream() << fluidTag << "sumMomentumFlux").str());
		writer.writeData(sumEnergyFluxes, (std::stringstream() << fluidTag << "sumEnergyFlux").str());


		Eigen::Matrix3Xd el_source = LorentzForce_electric.block(3 * fluidIdx, 0, 3, nCells);
		writer.writeData(el_source, fluidTag + "LorentzForceEl");
		Eigen::Matrix3Xd mag_source = LorentzForce_magnetic.block(3 * fluidIdx, 0, 3, nCells);
		writer.writeData(mag_source, fluidTag + "LorentzForceMag");


		if (isStateWrittenToOutput) {
			writer.writeData(qN, stateN);
			writer.writeData(qU, stateU);
			writer.writeData(qE, stateE);
		}
	
		assert(faceFluxes.cols() == nFaces);
		assert(faceFluxes.rows() == 5 * getNFluids());
		{
			Eigen::VectorXd faceFluxMass;
			faceFluxMass = faceFluxes.row(5 * fluidIdx + 0);
			assert(faceFluxMass.size() == nFaces);
			writer.writeData(faceFluxMass, (std::stringstream() << fluidTag << "massFlux").str());
		}
		{
			Eigen::MatrixXd faceFluxMomentum;
			faceFluxMomentum = faceFluxes.block(5 * fluidIdx + 1, 0, 3, nFaces);
			assert(faceFluxMomentum.rows() == 3);
			assert(faceFluxMomentum.cols() == nFaces);
			writer.writeData(faceFluxMomentum, (std::stringstream() << fluidTag << "momentumFlux").str());
		}
		{
			Eigen::VectorXd faceFluxEnergy;
			faceFluxEnergy = faceFluxes.row(5 * fluidIdx + 4);
			assert(faceFluxEnergy.size() == nFaces);
			writer.writeData(faceFluxEnergy, (std::stringstream() << fluidTag << "energyFlux").str());
		}
		{
			const Eigen::VectorXi faceTypeFluid = faceTypeFluids.col(fluidIdx);
			writer.writeData(faceTypeFluid, (std::stringstream() << fluidTag << "faceType").str());
		}
	}

	writer.writeData(faceFluxesImExRusanov, "/faceFluxesImExRusanov");
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
		E.col(i) = E_h(i) / edge->getLength() * edge->getDirection().normalized();
	}
	writer.writeData(E, "/E");

	writer.writeData(E_cc, "/Ecc");

	assert(J_h.size() > 0);
	assert(J_h.size() == dualMesh.getNumberOfFaces());
	assert(J_h.allFinite());
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
	assert(J_h_aux.allFinite());
	Eigen::Matrix3Xd J_h_aux_vector(3, J_h_aux.size());
	for (int i = 0; i < J_h_aux.size(); i++) {
		const Face * face = dualMesh.getFace(i);
		const double fA = face->getArea();
		const Eigen::Vector3d fn = face->getNormal();
		J_h_aux_vector.col(i) = J_h_aux(i) / fA * fn;
	}
	writer.writeData(J_h_aux_vector, "/j_h_aux_vector");
	writer.writeData(J_h_aux, "/J_h_aux");
	writer.writeData(J_h, "/J_h");

	Jcc = getCurrentDensityAtCellCenter();
	writer.writeData(Jcc, "/Jcc");

	// Interpolated values of B-field to primal vertices
	writer.writeData(B_vertex, "/Bcc");
	const int nCells = dualMesh.getNumberOfCells();
	Eigen::VectorXd Bcc_mag(nCells);
	for (int i = 0; i < nCells; i++) {
		Bcc_mag(i) = B_vertex.col(i).norm();
	}
	writer.writeData(Bcc_mag, "/Bcc_mag");

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
				(std::stringstream() << dataFilename << ":/Bcc").str()
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
				(std::stringstream() << dataFilename << ":/Bcc").str()
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
	appmParams = AppmParameters();

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
		if (tag == "timestepSize") {
			std::istringstream(line.substr(pos + 1)) >> appmParams.timestepSize;
		}
		if (tag == "isFluidEnabled") {
			std::istringstream(line.substr(pos + 1)) >> appmParams.isFluidEnabled;
		}
		if (tag == "isMaxwellEnabled") {
			std::istringstream(line.substr(pos + 1)) >> appmParams.isMaxwellEnabled;
		}
		if (tag == "lambdaSquare") {
			std::istringstream(line.substr(pos + 1)) >> this->lambdaSquare;
		}
		if (tag == "isMassFluxSchemeImplicit") {
			std::istringstream(line.substr(pos + 1)) >> appmParams.isMassFluxSchemeImplicit;
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
		if (tag == "isShowDataWriterOutput") {
			std::istringstream(line.substr(pos + 1)) >> isShowDataWriterOutput;
		}
		if (tag == "isEulerMaxwellCouplingEnabled") {
			std::istringstream(line.substr(pos + 1)) >> appmParams.isEulerMaxwellCouplingEnabled;
		}
		if (tag == "isMaxwellCurrentDefined") {
			std::istringstream(line.substr(pos + 1)) >> appmParams.isMaxwellCurrentDefined;
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
	std::cout << printSolverParameters() << std::endl;
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

	ss << "<Attribute Name=\"Vertex type\" AttributeType=\"Scalar\" Center=\"Node\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << primalMesh.getNumberOfVertices() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << datafilename << ":/vertexType" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;


	if (isWriteBfield) {
		ss << "<Attribute Name=\"Magnetic Flux\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << primalMesh.getNumberOfFaces() << " 3\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << "appm-" << iteration << ".h5:/B" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;
	}

	ss << "<Attribute Name=\"Magnetic Flux Vertex Centered\" AttributeType=\"Vector\" Center=\"Node\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << primalMesh.getNumberOfVertices() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "appm-" << iteration << ".h5:/Bcc" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"Magnetic Flux Vertex Centered Mag\" AttributeType=\"Scalar\" Center=\"Node\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfCells() << "\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "appm-" << iteration << ".h5" << ":/Bcc_mag" << std::endl;
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

	ss << "<Attribute Name=\"" << "Face Normal" << "\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << datafilename << ":/faceNormal" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;


	for (int fidx = 0; fidx < getNFluids(); fidx++) {
		std::string attributeName = (std::stringstream() << particleParams[fidx].name <<" Face Type").str();
		std::string dataName = (std::stringstream() << particleParams[fidx].name <<"-faceType").str();
		ss << "<Attribute Name=\"" << attributeName << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << "\""
			<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
		ss << "appm-" << iteration << ".h5:/" << dataName << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;

		attributeName = (std::stringstream() << particleParams[fidx].name << " Face Mass Flux").str();
		ss << "<Attribute Name=\"" << attributeName << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << "\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << "appm-" << iteration << ".h5:/" << particleParams[fidx].name << "-massFlux" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;

		attributeName = (std::stringstream() << particleParams[fidx].name << " Face Momentum Flux").str();
		ss << "<Attribute Name=\"" << attributeName << "\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << " 3\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << "appm-" << iteration << ".h5:/" << particleParams[fidx].name << "-momentumFlux" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;

		attributeName = (std::stringstream() << particleParams[fidx].name << " Face Energy Flux").str();
		ss << "<Attribute Name=\"" << attributeName << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << "\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << "appm-" << iteration << ".h5:/" << particleParams[fidx].name << "-energyFlux" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;
	}

	ss << "<Attribute Name=\"Current density\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "appm-" << iteration << ".h5:/CurrentDensity" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;
	

	ss << "<Attribute Name=\"Current density aux\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "appm-" << iteration << ".h5:/j_h_aux_vector" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"J h aux\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << "\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << "appm-" << iteration << ".h5:/J_h_aux" << std::endl;
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

	ss << "<Attribute Name=\"Magnetic Flux Cell Centered\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfCells() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << datafilename << ":/Bcc" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"Magnetic Flux Cell Centered Mag\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfCells() << "\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << datafilename << ":/Bcc_mag" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"Electric Field Cell Centered\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfCells() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << datafilename << ":/Ecc" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"Current Density Cell Centered\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfCells() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << datafilename << ":/Jcc" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

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
		const std::string fluidName = (std::stringstream() << particleParams[k].name).str();

		ss << "<Attribute Name=\"" << fluidName << " density" << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << nCells << "\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << datafilename << ":/" << fluidName << "-density" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;

		ss << "<Attribute Name=\"" << fluidName << " number density" << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << nCells << "\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << datafilename << ":/" << fluidName << "-numberDensity" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;

		ss << "<Attribute Name=\"" << fluidName << " pressure" << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << nCells << "\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << datafilename << ":/" << fluidName << "-pressure" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;


		ss << "<Attribute Name=\"" << fluidName << " velocity" << "\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << nCells << " 3\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << datafilename << ":/" << fluidName << "-velocity" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;

		ss << "<Attribute Name=\"Electric Lorentz Force\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfCells() << " 3\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << datafilename << ":/" << fluidName << "-LorentzForceEl" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;

		ss << "<Attribute Name=\"Magnetic Lorentz Force\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfCells() << " 3\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << datafilename << ":/" << fluidName << "-LorentzForceMag" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;

		ss << "<Attribute Name=\"" << fluidName << " SumMassFlux" << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << nCells << "\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << datafilename << ":/" << fluidName << "-sumMassFlux" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;

		ss << "<Attribute Name=\"" << fluidName << " SumMomentumFlux" << "\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfCells() << " 3\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << datafilename << ":/" << fluidName << "-sumMomentumFlux" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;

		ss << "<Attribute Name=\"" << fluidName << " SumEnergyFlux" << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << nCells << "\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << datafilename << ":/" << fluidName << "-sumEnergyFlux" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;


		if (isStateWrittenToOutput) {
			ss << "<Attribute Name=\"" << fluidName << " state N" << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
			ss << "<DataItem Dimensions=\"" << nCells << "\""
				<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
			ss << datafilename << ":/" << fluidName << "-stateN" << std::endl;
			ss << "</DataItem>" << std::endl;
			ss << "</Attribute>" << std::endl;


			ss << "<Attribute Name=\"" << fluidName << " stateU\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
			ss << "<DataItem Dimensions=\"" << nCells << " 3\""
				<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
			ss << datafilename << ":/" << fluidName << "-stateU" << std::endl;
			ss << "</DataItem>" << std::endl;
			ss << "</Attribute>" << std::endl;

			ss << "<Attribute Name=\"" << fluidName << " stateE\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
			ss << "<DataItem Dimensions=\"" << nCells << "\""
				<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
			ss << datafilename << ":/" << fluidName << "-stateE" << std::endl;
			ss << "</DataItem>" << std::endl;
			ss << "</Attribute>" << std::endl;
		}
	}
	return ss.str();
}

const Face::Type AppmSolver::getFaceTypeOfFluid(const Face * face, const int fluidIdx) const
{
	Face::Type faceTypeOfFluid = face->getType();
	if (faceTypeOfFluid == Face::Type::TERMINAL) {
		//faceTypeOfFluid = (particleParams[fluidIdx].electricCharge >= 0) ? Face::Type::WALL : Face::Type::OPENING;
		faceTypeOfFluid = Face::Type::OPENING; // for testing
	}
	return faceTypeOfFluid;
}

const Eigen::VectorXd AppmSolver::solveMaxwell_PardisoLU(Eigen::SparseMatrix<double>& Mf, Eigen::VectorXd & rhs)
{
	std::cout << "Setup solver for Maxwell system: PardisoLU" << std::endl;
	Eigen::VectorXd xf(rhs.size());
	Eigen::PardisoLU<Eigen::SparseMatrix<double>> maxwellSolver;

	auto timer_start = std::chrono::high_resolution_clock::now();
	maxwellSolver.compute(Mf);
	auto timer_afterCompute = std::chrono::high_resolution_clock::now();

	auto delta_afterCompute = std::chrono::duration<double>(timer_afterCompute - timer_start);
	std::cout << "Solver time for compute step: " << delta_afterCompute.count() << std::endl;

	if (maxwellSolver.info() != Eigen::Success) {
		std::cout << "Maxwell solver setup failed" << std::endl;
	}

	auto timer_startSolve = std::chrono::high_resolution_clock::now();
	xf = maxwellSolver.solve(rhs);
	auto timer_endSolve = std::chrono::high_resolution_clock::now();

	std::cout << "Maxwell solver finished" << std::endl;

	auto delta_solve = std::chrono::duration<double>(timer_endSolve - timer_startSolve);
	std::cout << "Solver time for solve step: " << delta_solve.count() << std::endl;

	Eigen::VectorXd residual = rhs - Mf * xf;
	double rms = std::sqrt(residual.cwiseAbs2().sum() / residual.size());
	std::cout << "maxRes: \t" << residual.cwiseAbs().maxCoeff() << std::endl;
	std::cout << "rms:    \t" << rms << std::endl;
	if (maxwellSolver.info() != Eigen::Success) {
		std::cout << "Maxwell solver solving failed" << std::endl;
	}

	assert(maxwellSolver.info() == Eigen::Success);
	return xf;
}

const Eigen::VectorXd AppmSolver::solveMaxwell_sparseLU(Eigen::SparseMatrix<double>& Mf, Eigen::VectorXd & rhs)
{
	std::cout << "Setup solver for Maxwell system: SparseLU" << std::endl;
	Eigen::VectorXd xf(rhs.size());
	Eigen::SparseLU<Eigen::SparseMatrix<double>> maxwellSolver; 

	auto timer_start = std::chrono::high_resolution_clock::now();
	maxwellSolver.compute(Mf);
	auto timer_afterCompute = std::chrono::high_resolution_clock::now();

	auto delta_afterCompute = std::chrono::duration<double>(timer_afterCompute - timer_start);
	std::cout << "Solver time for compute step: " << delta_afterCompute.count() << std::endl;
	
	if (maxwellSolver.info() != Eigen::Success) {
		std::cout << "Maxwell solver setup failed" << std::endl;
	}

	auto timer_startSolve = std::chrono::high_resolution_clock::now();
	xf = maxwellSolver.solve(rhs);
	auto timer_endSolve = std::chrono::high_resolution_clock::now();
	
	std::cout << "Maxwell solver finished" << std::endl;
	
	auto delta_solve = std::chrono::duration<double>(timer_endSolve - timer_startSolve);
	std::cout << "Solver time for solve step: " << delta_solve.count() << std::endl;
	
	if (maxwellSolver.info() != Eigen::Success) {
		std::cout << "Maxwell solver solving failed" << std::endl;
	}
	
	assert(maxwellSolver.info() == Eigen::Success);
	return xf;
}

const Eigen::VectorXd AppmSolver::solveMaxwell_BiCGStab(Eigen::SparseMatrix<double>& Mf, Eigen::VectorXd & rhs)
{
	std::cout << "Setup solver for Maxwell system: BiCGStab" << std::endl;
	Eigen::VectorXd xf(rhs.size());
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> maxwellSolver;
	maxwellSolver.preconditioner().setDroptol(0.1);
	maxwellSolver.preconditioner().setFillfactor(0.1);

	auto timer_start = std::chrono::high_resolution_clock::now();
	maxwellSolver.compute(Mf);
	auto timer_afterCompute = std::chrono::high_resolution_clock::now();

	auto delta_afterCompute = std::chrono::duration<double>(timer_afterCompute - timer_start);
	std::cout << "Solver time for compute step: " << delta_afterCompute.count() << std::endl;

	if (maxwellSolver.info() != Eigen::Success) {
		std::cout << "Maxwell solver setup failed" << std::endl;
	}

	auto timer_startSolve = std::chrono::high_resolution_clock::now();
	xf = maxwellSolver.solve(rhs);
	auto timer_endSolve = std::chrono::high_resolution_clock::now();

	std::cout << "Maxwell solver finished" << std::endl;

	auto delta_solve = std::chrono::duration<double>(timer_endSolve - timer_startSolve);
	std::cout << "Solver time for solve step: " << delta_solve.count() << std::endl;

	if (maxwellSolver.info() != Eigen::Success) {
		std::cout << "Maxwell solver solving failed" << std::endl;
	}

	std::cout << "# Iterations: " << maxwellSolver.iterations() << std::endl;
	std::cout << "error:        " << maxwellSolver.error() << std::endl;
	
	assert(maxwellSolver.info() == Eigen::Success);
	return xf;
}

const Eigen::VectorXd AppmSolver::solveMaxwell_LSCG(Eigen::SparseMatrix<double>& Mf, Eigen::VectorXd & rhs)
{
	std::cout << "Setup solver for Maxwell system: LeastSquaresCG" << std::endl;
	Eigen::VectorXd xf(rhs.size());
	Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> maxwellSolver;
	maxwellSolver.setTolerance(Eigen::NumTraits<double>::epsilon() * 1024);

	auto timer_start = std::chrono::high_resolution_clock::now();
	maxwellSolver.compute(Mf);
	auto timer_afterCompute = std::chrono::high_resolution_clock::now();

	auto delta_afterCompute = std::chrono::duration<double>(timer_afterCompute - timer_start);
	std::cout << "Solver time for compute step: " << delta_afterCompute.count() << std::endl;

	if (maxwellSolver.info() != Eigen::Success) {
		std::cout << "Maxwell solver setup failed" << std::endl;
	}

	auto timer_startSolve = std::chrono::high_resolution_clock::now();
	xf = maxwellSolver.solve(rhs);
	auto timer_endSolve = std::chrono::high_resolution_clock::now();

	std::cout << "Maxwell solver finished" << std::endl;

	auto delta_solve = std::chrono::duration<double>(timer_endSolve - timer_startSolve);
	std::cout << "Solver time for solve step: " << delta_solve.count() << std::endl;

	if (maxwellSolver.info() != Eigen::Success) {
		std::cout << "Maxwell solver solving failed" << std::endl;
	}

	std::cout << "# Iterations: " << maxwellSolver.iterations() << std::endl;
	std::cout << "error:        " << maxwellSolver.error() << std::endl;

	assert(maxwellSolver.info() == Eigen::Success);
	return xf;
}

const Eigen::VectorXd AppmSolver::solveMaxwell_CG(Eigen::SparseMatrix<double>& Mf, Eigen::VectorXd & rhs)
{
	std::cout << "Setup solver for Maxwell system: ConjugateGradient" << std::endl;
	Eigen::VectorXd xf(rhs.size());
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> maxwellSolver;
	maxwellSolver.setTolerance(Eigen::NumTraits<double>::epsilon() * 1024);

	auto timer_start = std::chrono::high_resolution_clock::now();
	maxwellSolver.compute(Mf);
	auto timer_afterCompute = std::chrono::high_resolution_clock::now();

	auto delta_afterCompute = std::chrono::duration<double>(timer_afterCompute - timer_start);
	std::cout << "Solver time for compute step: " << delta_afterCompute.count() << std::endl;

	if (maxwellSolver.info() != Eigen::Success) {
		std::cout << "Maxwell solver setup failed" << std::endl;
	}

	auto timer_startSolve = std::chrono::high_resolution_clock::now();
	xf = maxwellSolver.solve(rhs);
	auto timer_endSolve = std::chrono::high_resolution_clock::now();

	std::cout << "Maxwell solver finished" << std::endl;

	auto delta_solve = std::chrono::duration<double>(timer_endSolve - timer_startSolve);
	std::cout << "Solver time for solve step: " << delta_solve.count() << std::endl;

	if (maxwellSolver.info() != Eigen::Success) {
		std::cout << "Maxwell solver solving failed" << std::endl;
	}

	std::cout << "# Iterations: " << maxwellSolver.iterations() << std::endl;
	std::cout << "error:        " << maxwellSolver.error() << std::endl;

	assert(maxwellSolver.info() == Eigen::Success);
	return xf;
}

/** 
* @return matrix Meps such that: J_h = Meps * E_h + J_aux. 
*/
Eigen::SparseMatrix<double> AppmSolver::get_Msigma_spd(Eigen::VectorXd & Jaux, const double dt, const double time)
{
	const int nDualFaces = dualMesh.getNumberOfFaces();
	const int nEdges = primalMesh.getNumberOfEdges();
	Eigen::SparseMatrix<double> Msigma(nDualFaces, nEdges);
	Jaux.setZero();

	if (appmParams.isMaxwellCurrentDefined) {
		// Define a uniform current density for those faces being normal to z-axis
		const double t0 = 2;
		const double tscale = 0.5;
		const double currentDensity = 0.5 * (1 + tanh((time + dt - t0) / tscale)); // current density at implicit timestep t + dt
		const double radius = 0.5;
		const double tol = 2 * std::numeric_limits<double>::epsilon();

		for (int i = 0; i < nEdges; i++) {
			const Face * dualFace = dualMesh.getFace(i);
			const Eigen::Vector3d fn = dualFace->getNormal();
			const double fA = dualFace->getArea();
			const Eigen::Vector3d fc = dualFace->getCenter();

			const bool isZnormalFace = std::abs(fn(2)) > tol; // face is normal to z-axis 
			const bool isInsideRadius = fc.segment(0, 2).norm() < radius; // distance to z-axis is small enough
			if (isZnormalFace && isInsideRadius) { 
				Jaux(i) = currentDensity * fA * fn.dot(Eigen::Vector3d::UnitZ());
			}
		}
		//return Msigma;
	}

	// Check if primal edges and dual face normals are oriented in same direction
	const bool debug = false;
	if (debug) {
		Eigen::VectorXd orientation(nEdges);
		orientation.setZero();
		for (int i = 0; i < nEdges; i++) {
			const Edge * edge = primalMesh.getEdge(i);
			const Eigen::Vector3d Lhat = edge->getDirection().normalized();

			const Face * face = dualMesh.getFace(i);
			const Eigen::Vector3d nhat = face->getNormal();

			orientation(i) = Lhat.dot(nhat); // this value should be equal to 1.000
		}
		orientation.array() -= 1;
		const double tol = 4 * std::numeric_limits<double>::epsilon();
		assert((orientation.cwiseAbs().array() <= tol).all());
	}

	std::cout << "isEulerMaxwellCouplingEnabled: " << appmParams.isEulerMaxwellCouplingEnabled << std::endl;
	if (appmParams.isEulerMaxwellCouplingEnabled) {
		const int nFluids = getNFluids();

		typedef Eigen::Triplet<double> T;
		std::vector<T> triplets;
		std::vector<T> geomTriplets;

		// Loop over all dual faces (where current density lives)
		for (int i = 0; i < nDualFaces; i++) {
			//std::cout << "i = " << i << std::endl;
			const Face * dualFace = dualMesh.getFace(i);
			const Eigen::Vector3d ni = dualFace->getNormal();
			const double Ai = dualFace->getArea();
			const Face::Type faceType = dualFace->getType();
			const double numSchemeFactor = (faceType == Face::Type::INTERIOR) ? 0.5 : 1;

			// Do not consider dual boundary faces (those coincident with domain boundary) ...
			// ... this leaves us with inner faces of the dual mesh 
			if (i >= nEdges) { 
				assert(!dualFace->hasFluidCells());
				continue; 
			}

			for (int fluidIdx = 0; fluidIdx < nFluids; fluidIdx++) {
				const int q = particleParams[fluidIdx].electricCharge;
				const double eps2 = particleParams[fluidIdx].mass;
				const Face::Type faceFluidType = getFaceTypeOfFluid(dualFace, fluidIdx);

				// Term due to explicit mass flux
				assert(i < faceFluxes.cols());
				const double massFlux = faceFluxes(5 * fluidIdx + 0, i);
				Jaux(i) += q * Ai * massFlux;


				// Add extra terms due to implicit formulation
				if (appmParams.isMassFluxSchemeImplicit) {
					auto adjacientCells = dualFace->getCellList();
					for (auto cell : adjacientCells) {
						// skip non-fluid cells
						if (cell->getType() != Cell::Type::FLUID) {
							continue;
						}
						assert(cell->getType() == Cell::Type::FLUID);
						const double nk = fluidStates(5 * fluidIdx + 0, cell->getIndex()); // number density in cell k
						const double Vk = cell->getVolume(); // volume of cell k

						// Extra term by explicit magnetic Lorentz force (magnetic flux)
						if (appmParams.isLorentzForceMagneticEnabled) {
							const Eigen::Vector3d nu_x_B = LorentzForce_magnetic.col(cell->getIndex()).segment(3 * fluidIdx, 3);
							Jaux(i) += numSchemeFactor * pow(q, 2) * 1. / eps2 * Ai * nu_x_B.dot(ni);
						}

						auto cellFaces = cell->getFaceList();
						for (auto face : cellFaces) {
							const double Aj = face->getArea();
							const Eigen::Vector3d nj = face->getNormal();
							const int s_kj = cell->getOrientation(face);
							const int j = face->getIndex();
							double nj_dot_ni = nj.dot(ni); // angle between face i and face j
							if (abs(nj_dot_ni) < 1e-10) { 
								nj_dot_ni = 0; // Truncate small angles
							}

							// Extra term by explicit momentum flux
							const double fj = faceFluxes.col(j).segment(5 * fluidIdx + 1, 3).dot(nj); // momentum flux at face j
							Jaux(i) += numSchemeFactor * q * Ai * -dt * 1 / Vk * s_kj * Aj * fj * nj_dot_ni;
							
							// Extra term by implicit electric Lorentz force (electric field)
							const double geomFactor = nj_dot_ni * (Ai * Aj) / Vk; // geometric factor
							const double value = numSchemeFactor * 0.5 * dt * pow(q, 2) / eps2 * nk / Vk * Ai * Aj * nj_dot_ni;
							triplets.push_back(T(i, j, value));
							if (fluidIdx == 0) {
								geomTriplets.push_back(T(i, j, geomFactor));
							}
						}
					}
				} // end if isMassFluxSchemeImplicit
			}
		}
		Msigma.setFromTriplets(triplets.begin(), triplets.end());
		Msigma.makeCompressed();

		//Eigen::SparseMatrix<double> Msigma_geom(nDualFaces, nEdges);
		//Msigma_geom.setFromTriplets(geomTriplets.begin(), geomTriplets.end());
		//Msigma_geom.makeCompressed();
		//std::cout << "Msigma_geom triplets size: " << geomTriplets.size() << std::endl;
		//std::cout << "Msigma_geom nnz: " << Msigma_geom.nonZeros() << std::endl;
		//Eigen::sparseMatrixToFile(Msigma_geom, "Msigma_geom.dat");
	}

	return Msigma;
}

/**
* Testcase: apply energy source if cell center is inside a given radius and time smaller than a threshold.
*/
const Eigen::VectorXd AppmSolver::testcase_001_FluidSourceTerm(const double time, const Cell * cell, const int fluidIdx) const
{
	Eigen::VectorXd srcVector(5);
	srcVector.setZero();
	const double radius = 0.2;
	const double t0 = 1;

	if (cell->getCenter().norm() < radius && time < t0) {
		srcVector(4) = 10;
	}
	return srcVector;
}

const Eigen::VectorXd AppmSolver::setVoltageBoundaryConditions(const double time) const
{
	const int nPrimalTerminalVertices = primalMesh.getMeshInfo().nVerticesTerminal;
	assert(nPrimalTerminalVertices > 0);

	int nDirichlet = nPrimalTerminalVertices;
	if (appmParams.isMaxwellCurrentDefined) {
		nDirichlet = nPrimalTerminalVertices;
	}

	Eigen::VectorXd xd = Eigen::VectorXd::Zero(nDirichlet);

	if (appmParams.isMaxwellCurrentDefined) {
		xd.topRows(nPrimalTerminalVertices / 2).array() = 0;
	}
	else {
		// Set voltage values at terminal A (stressed electrode)
		xd.topRows(nPrimalTerminalVertices / 2).array() = 0.5 * (1 + tanh((time - 1) / 0.1));
	}
	// Set voltage condition to zero at terminal B (grounded electrode)
	xd.bottomRows(nPrimalTerminalVertices / 2).array() = 0;

	return xd;
}

