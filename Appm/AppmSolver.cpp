#include "AppmSolver.h"



AppmSolver::AppmSolver() 
{
}


AppmSolver::~AppmSolver()
{
	for (int i = elasticCollisions.size() - 1; i >= 0; i--) {
		delete elasticCollisions[i];
		elasticCollisions[i] = nullptr;
	}
}

/**
* @return get message how data should be visualized.
*/
const std::string AppmSolver::message_howToVisualizeData(const std::string & outFilename_volume, const std::string & outFilename_surface) const
{
	std::stringstream ss;
	ss << "Note on data output: \n";
	ss << " - heavy data is stored in HDF5 files. \n";
	ss << " - light data is stored in XDMF files (version 3).  \n";
	ss << " (Polyhedrons have been added to XDMF in version 3; this was not supported in version 2.)  \n";
	ss << "  \n";
	ss << "Use Paraview (version 5.6.0) to visualize. It comes with three readers: \n";
	ss << " - XDMF Reader:  CANNOT read polyhedrons, but  CAN   read grid-of-grids. \n";
	ss << " - Xdmf3ReaderS:  CAN   read polyhedrons, but CANNOT grid-of-grids (i.e., grids with GridType=Tree). \n";
	ss << " - Xdmf3ReaderT: ??? \n";
	ss << "  \n";
	ss << "Therefore output data has been splitted into two *.xdmf files:  \n";
	ss << " - " << outFilename_surface << ":   contains data with geometric dimensionality <= 2 (i.e., vertices, edges, faces) \n";
	ss << " - " << outFilename_volume << ": contains volumetric data \n";
	ss << "  \n";
	ss << "To visualize volumetric data: \n";
	ss << " - Start paraview, open " << outFilename_volume << ", and read it with Xdmf3ReaderS. \n";
	ss << "  \n";
	ss << "To visualize point, edge, or surface data:  \n";
	ss << " - Start paraview, open " << outFilename_surface << ", and read it with XDMF Reader. \n";
	ss << "  \n";
	return ss.str();
}


void AppmSolver::init()
{
	// Show solver parameters
	std::cout << solverParams << std::endl;

	// Show mesh parameters
	std::cout << primalMeshParams << std::endl;
	std::cout << "Species list:" << std::endl;
	for (int i = 0; i < speciesList.size(); i++) {
		std::cout << "Species " << i << ":" << std::endl;
		std::cout << speciesList[i] << std::endl;
		std::cout << std::endl;
	}

	// Output file for time values at each timestep
	this->timeFile = std::ofstream("timesteps.dat");

	isStateWrittenToOutput = true;
	//readParameters("AppmSolverParams.txt");
	init_meshes(this->primalMeshParams);  // Initialize primal and dual meshes
	if (primalMesh.getNumberOfCells() == 0) {
		return;
	}
	if (dualMesh.getNumberOfCells() == 0) {
		return;
	}

	if (dualMesh.getNumberOfCells() != primalMesh.getNumberOfVertices()) {
		std::cout << "Number of dual cells is not equal to primal vertices" << std::endl;
	}

	init_multiFluid();

	
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


	if (solverParams.getMaxwellEnabled()) {
		isWriteBfield = true;
		isWriteEfield = true;
	}

	int iteration = 0;
	double time = 0;
	double dt = 1;

	//set_Efield_uniform(Eigen::Vector3d::UnitZ());
	//E_cc = getEfieldAtCellCenter();
	//setFluidFaceFluxes();
	//dt = 1;
	//dt = getNextFluidTimestepSize();
	//Eigen::SparseMatrix<double> Msigma;
	//Msigma = get_Msigma_spd(J_h_aux, dt, time);
	//Eigen::sparseMatrixToFile(Msigma, "Msigma.dat");

	//Eigen::isSymmetricPositiveDefinite(Msigma);
	//std::cout << "nnz(Msigma): " << Msigma.nonZeros() << std::endl;
	//J_h = Msigma * E_h + J_h_aux;
	//Jcc = getCurrentDensityAtCellCenter();

	//set_Bfield_azimuthal();
	//interpolateMagneticFluxToPrimalVertices();

	// write initial data to file (iteration = 0, time = 0)

	//iteration = 0;
	//time = 0;
	//writeOutput(iteration, time);

	//J_h_aux.setZero();
	//J_h.setZero();
	//E_h.setZero();
	//E_cc = getEfieldAtCellCenter();
}

std::string AppmSolver::getIterationHeader(const int iter, const double time, const double dt) const
{
	std::stringstream ss;
	ss << "*********************************************" << std::endl;
	ss << "* Iteration " << iter;
	ss << ",\t time = " << time;
	ss << ",\t dt = " << dt;
	ss << std::endl;
	ss << "*********************************************";
	return ss.str();
}



void AppmSolver::run()
{
	init();
	
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

	double dt = solverParams.getMaxTimestepSize();
	double dt_previous = dt;
	double time = 0;
	int iteration = 0;


	createStopFile(0);

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
	if (solverParams.getMaxIterations() > 0) {
		setFluidFaceFluxes();
	}


	/*
	* Time integration loop 
	*/
	const int maxIterations = solverParams.getMaxIterations();
	const double maxTime = solverParams.getMaxTime();
	const int outputFrequency = solverParams.getOutputFrequency();
	std::cout << getIterationHeader(iteration, time, dt) << std::endl;

	writeOutput(iteration, time);
	
	try {
		while (iteration < maxIterations && time < maxTime) {
			// Initialize data
			fluidSources.setZero();
			faceFluxes.setZero();

			// Determine timestep
			dt_previous = dt;
			if (solverParams.getFluidEnabled()) {
				dt = getNextFluidTimestepSize();
				// Limit timestep such that maximum time is reached exactly
				if (time + dt >= maxTime) {
					std::cout << "timestep limited to reach maximum time; value before limiter is applied: dt = " << dt << std::endl;
					dt = maxTime - time;
				}
			}

			time += dt;
			iteration += 1;

			std::cout << std::endl;
			std::cout << getIterationHeader(iteration, time, dt) << std::endl;

			// Set explicit fluid source terms
			setMagneticLorentzForceSourceTerms();
			//setFrictionSourceTerms();
			setElasticCollisionSourceTerms();
			//setInelasticCollisionSourceTerms();
			//setRadiationSource(); // <<<----  TODO

			// Set explicit face fluxes
			if (solverParams.getFluidEnabled()) {
				setFluidFaceFluxes();
			}

			// Affine-linear function for implicit-consistent formulation of current density J_h = Msigma * E_h + Jaux
			Eigen::SparseMatrix<double> Msigma;
			Msigma = get_Msigma_spd(J_h_aux, dt, time);

			// Maxwell equations
			if (solverParams.getMaxwellEnabled()) {
				solveMaxwellSystem(time, dt, dt_previous, Msigma);
				// E- and B-field have been updated to new timestep
			}
			J_h = Msigma * E_h + J_h_aux;

			// Fluid equations
			if (solverParams.getFluidEnabled()) {
				const int nFluids = this->getNFluids();
				const int nCells = dualMesh.getNumberOfCells();
				const int nFaces = dualMesh.getNumberOfFaces();

				// Set source term for electric Lorentz force (implicit)
				if (solverParams.getLorentzForceEnabled()) {
					for (int i = 0; i < nCells; i++) {
						const Cell * cell = dualMesh.getCell(i);
						if (cell->getType() != Cell::Type::FLUID) {
							continue; // Skip cells that are not of type Fluid
						}
						for (int fluidIdx = 0; fluidIdx < nFluids; fluidIdx++) {
							//const int q = particleParams[fluidIdx].electricCharge;
							const int q = getSpecies(fluidIdx).getCharge();
							const double n = fluidStates(5 * fluidIdx, i);
							const double massRatio = getSpecies(fluidIdx).getMassRatio();
							LorentzForce_electric.col(i).segment(3 * fluidIdx, 3) = q * 1. / massRatio * n * E_cc.col(i);
						}
					}

					// Update momentum source terms ...
					for (int fluidIdx = 0; fluidIdx < nFluids; fluidIdx++) {
						fluidSources.block(5 * fluidIdx + 1, 0, 3, nCells) += LorentzForce_electric.block(3 * fluidIdx, 0, 3, nCells);
					}
					// Update energy source terms ...
					for (int fluidIdx = 0; fluidIdx < nFluids; fluidIdx++) {
						int q = getSpecies(fluidIdx).getCharge();
						for (int i = 0; i < nCells; i++) {
							Eigen::VectorXd state = fluidStates.col(i);
							Eigen::Vector3d nu = state.segment(5 * fluidIdx + 1, 3); // mass flux in cell i at timestep m; 
							Eigen::Vector3d E = E_cc.col(i); // electric field in cell i at timestep m+1
							double ohmicSource = nu.dot(E);
							ohmicSource *= q;
							fluidSources(5 * fluidIdx + 4, i) += ohmicSource;
						}
					}
				}
				std::cout << "Electric Lorentz Force maxCoeff: " << LorentzForce_electric.cwiseAbs().maxCoeff() << std::endl;

				//setFluidSourceTerm();

				// Set implicit terms for mass fluxes
				if (solverParams.getMassfluxSchemeImplicit()) {
					setImplicitMassFluxTerms(dt);
				}

				setSumOfFaceFluxes();
				const bool isImplicitSources = solverParams.getEulerSourcesImplicit();
				updateFluidStates(dt, isImplicitSources);

				//debug_checkCellStatus();
			}

			// Add this time value to list of timesteps
			timeFile << iteration << ", " << time << ", " << dt << std::endl;


			// Criteria to stop loop iterations
			const bool isStopFile = isStopFileActive();
			const bool isMaxIters = iteration >= maxIterations;
			const bool isMaxTime = time >= maxTime;
			const bool isStop = isStopFile || isMaxIters || isMaxTime;

			if (!isStop) {
				if (outputFrequency > 0) {
					if (iteration % outputFrequency == 0) {
						writeOutput(iteration, time);
					}
					// Write intermediate output file in case that the simulation crashes
					//if (iteration % (10 * outputFrequency) == 0) {
					//	writeXdmf("appm.xdmf");
					//	writeXdmfDualVolume("appm-volume.xdmf");
					//}
				}
			}
			else {
				writeOutput(iteration, time);
				writeXdmf("appm.xdmf");
				writeXdmfDualVolume("appm-volume.xdmf");
				break; // <<<--- stop loop iterations
			}
		}
	}
	catch (const std::exception & e) {
		std::cout << "**********************************************************" << std::endl;
		std::cout << "**********************************************************" << std::endl;
		std::cout << "**********************************************************" << std::endl;
		std::cout << "AppmSolver failed; error stack is: " << std::endl;
		std::cout << e.what() << std::endl;
		std::cout << "**********************************************************" << std::endl;
		std::cout << "**********************************************************" << std::endl;
		std::cout << "**********************************************************" << std::endl;
		writeOutput(iteration, time);
	}
	const std::string filename_surfaceData = "appm.xdmf";
	const std::string filename_volumeData = "appm-volume.xdmf";
	writeXdmf(filename_surfaceData);
	writeXdmfDualVolume(filename_volumeData);

	std::cout << std::endl;
	std::cout << "Final time:      " << time << std::endl;
	std::cout << "Final iteration: " << iteration << std::endl;

	auto timer_endAppmSolver = std::chrono::high_resolution_clock::now();
	auto delta_appmSolver = std::chrono::duration<double>(timer_endAppmSolver - timer_startAppmSolver);
	std::cout << "Elapsed time for APPM solver: " << delta_appmSolver.count() << std::endl;

	std::cout << solverParams << std::endl;

	std::cout << message_howToVisualizeData(filename_volumeData, filename_surfaceData) << std::endl;
}


void AppmSolver::setElasticCollisionSourceTerms()
{
	std::cout << "Elastic collisions " << std::endl;
	if (!solverParams.getFluidEnabled()) {
		return;
	}
	// Number of elastic collisions
	auto nElColls = elasticCollisions.size();
	if (nElColls <= 0) {
		return;
	}
	
	// Number of dual cells
	auto nCells = dualMesh.getNumberOfCells();

	// Check size and reset source terms
	const int nFluids = getNFluids();
	Eigen::MatrixXd elCollMomSourceTerm(3*nFluids, nCells);
	elCollMomSourceTerm.setZero();
	Eigen::MatrixXd elCollEnergySourceTerm(nFluids, nCells);
	elCollEnergySourceTerm.setZero();

	// Reduced temperature for each collision and each fluid cell
	Eigen::MatrixXd Tab(nElColls, nCells);
	//Tab.setZero();
	Tab.setConstant(std::nan(""));

	// Get reduced temperature for each fluid cell ...
	for (int i = 0; i < nCells; i++) {
		const Cell * cell = dualMesh.getCell(i);
		if (cell->getType() != Cell::Type::FLUID) {
			// skip non-fluid cells
			continue;
		}
		// ... and each collision
		for (int collisionIdx = 0; collisionIdx < nElColls; collisionIdx++) {
			ElasticCollision * coll = elasticCollisions[collisionIdx];
			const int idxA = coll->getFluidIdxA();
			const int idxB = coll->getFluidIdxB();
			const Eigen::VectorXd stateA = getState(i, idxA);
			const Eigen::VectorXd stateB = getState(i, idxB);

			const double ma = getSpecies(idxA).getMassRatio();
			const double mb = getSpecies(idxB).getMassRatio();
			const double Ta = Physics::getTemperature(stateA, ma);
			const double Tb = Physics::getTemperature(stateB, mb);

			const double Treduced = (ma * Tb + mb * Ta) / (ma + mb);
			Tab(collisionIdx, i) = Treduced;
		}
	}

	// Compute rate of change in momentum and energy due to collisions 
	for (int collIdx = 0; collIdx < nElColls; collIdx++) {
		ElasticCollision * coll = elasticCollisions[collIdx];
		Eigen::VectorXd T = Tab.row(collIdx);
		Eigen::VectorXd Qbar11 = coll->getAvgMomCrossSection(T);

		// get reduced mass of this collision
		const int idxA = coll->getFluidIdxA();
		const int idxB = coll->getFluidIdxB();
		const double ma = getSpecies(idxA).getMassRatio();
		const double mb = getSpecies(idxB).getMassRatio();
		const double mab = ma * mb / (ma + mb);

		// number densities of fluid A and B
		Eigen::VectorXd na = fluidStates.row(5 * idxA);
		Eigen::VectorXd nb = fluidStates.row(5 * idxB);

		// velocity of fluid A and B
		Eigen::MatrixXd ua = fluidStates.block(5 * idxA + 1, 0, 3, nCells);
		Eigen::MatrixXd ub = fluidStates.block(5 * idxB + 1, 0, 3, nCells);


		// factor in source term
		Eigen::VectorXd meanRelVelocity;
		meanRelVelocity = (1. / M_PI * 8. / mab * T.array()).sqrt();
		
		// Momentum rate of change due to elastic collisions;
		// temp_ab = na * mab * tau_ab^-1
		Eigen::VectorXd temp_ab;
		temp_ab = na.array() * nb.array() * (mab * 4. / 3. * meanRelVelocity.array() * Qbar11.array());

		// check data size
		assert(temp_ab.cols() == 1);
		assert(temp_ab.size() == nCells);

		std::cout << "Set default value: 1" << std::endl;
		temp_ab.setOnes(); // <<<----------------------------------------  TODO: default value! ---------------
		assert(false);


		// for each fluid cell ...
		for (int i = 0; i < nCells; i++) {
			if (dualMesh.getCell(i)->getType() == Cell::Type::FLUID) {
				const Eigen::VectorXd stateA = getState(i, idxA);
				const Eigen::VectorXd stateB = getState(i, idxB);

				// temperature of fluid A and B
				double Ta = Physics::getTemperature(stateA, ma);
				double Tb = Physics::getTemperature(stateB, mb);

				// ... apply momentum source term
				const Eigen::Vector3d R_local = (-1) * temp_ab(i) * (ua.col(i) - ub.col(i));
				elCollMomSourceTerm.col(i).segment(3 * idxA, 3) += R_local;
				elCollMomSourceTerm.col(i).segment(3 * idxB, 3) -= R_local;

				// ... apply energy source term 
				const double Q_local = (-1) * 3. / (ma + mb) * temp_ab(i) * (Ta - Tb);
				elCollEnergySourceTerm(idxA, i) += Q_local;
				elCollEnergySourceTerm(idxB, i) -= Q_local;
			}
		}
	}

	// Apply elastic collision sources to fluid source term
	for (int i = 0; i < nCells; i++) {
		for (int fluidx = 0; fluidx < nFluids; fluidx++) {
			fluidSources.col(i).segment(5 * fluidx + 1, 3) += elCollMomSourceTerm.col(i).segment(3 * fluidx, 3);
			fluidSources(5 * fluidx + 4, i) += elCollEnergySourceTerm(fluidx, i);
		}
	}
	std::cout << "max elastic momentum source term: " << elCollMomSourceTerm.array().abs().maxCoeff() << std::endl;
	std::cout << "max elastic energy   source term: " << elCollEnergySourceTerm.array().abs().maxCoeff() << std::endl;
	std::cout << "Elastic collisions DONE" << std::endl;
}


void AppmSolver::setSolverParameters(const AppmSolver::SolverParameters & solverParams)
{
	this->solverParams = solverParams;
}

void AppmSolver::setMeshParameters(const PrimalMesh::PrimalMeshParams & meshParams)
{
	this->primalMeshParams = meshParams;
}

void AppmSolver::setSpecies(const std::vector<Species>& speciesList)
{
	this->speciesList = speciesList;
}

void AppmSolver::setElasticCollisions(const std::vector<std::string> & list)
{
	const double Tscale = scalingParameters.getTemperatureScale();
	const double nScale = scalingParameters.getNumberDensityScale();
	const double xScale = scalingParameters.getLengthScale();
	assert(nScale > 0);
	assert(xScale > 0);
	const double crossSectionScale = scalingParameters.getCrossSectionsScale();
	assert(crossSectionScale > 0);
	std::cout << "Set elastic collisions" << std::endl;
	std::cout << "Scaling variables:" << std::endl;
	std::cout << "  Tscale: " << Tscale << std::endl;
	std::cout << "  nScale: " << nScale << std::endl;
	std::cout << "  xScale: " << xScale << std::endl;
	std::cout << "  Qscale: " << crossSectionScale << std::endl;


	this->elasticCollisions = std::vector<ElasticCollision*>();
	for (auto tag : list) {
		int pos = tag.find('-');
		const std::string tagA = tag.substr(0, pos);
		const std::string tagB = tag.substr(pos + 1);
		int idxA = getSpeciesIndex(tagA);
		int idxB = getSpeciesIndex(tagB);

		std::stringstream ss;
		ss << "collisions/elastic/" << tag << ".dat";
		const std::string filename = ss.str();

		ElasticCollision * elasticCollision = new ElasticCollision(filename, idxA, idxB, Tscale, crossSectionScale);
		{
			std::stringstream ss;
			ss << "elastic-" << tag << ".dat";
			std::string filename = ss.str();
			Eigen::MatrixXd data = elasticCollision->getData();
			std::ofstream file(filename);
			std::cout << "Write data file: " << filename << std::endl;
			file << data << std::endl;
		}
		elasticCollisions.push_back(elasticCollision);
	}


}

/**
* Set list of inelastic collisions.
*
* Inelastic collision are restricted to ionization (and its reverse, recombination) reaction. 
* The parameter file is identified by the name of electron and atom, e.g., e-Ar.dat
* The electron must be in first place. We assume positive ions.
*
* The file has two columns: electron temperature (Te) and ionization rate (ki).
* 
*/
void AppmSolver::setInelasticCollisions(const std::vector<std::string>& list)
{
	std::cout << "Set inelastic collisions" << std::endl;

	this->inelasticCollisions = std::vector<InelasticCollision*>();
	for (auto tag : list) {
		int pos = tag.find('-');
		const std::string tagElectron = tag.substr(0, pos);
		const std::string tagAtom = tag.substr(pos + 1);
		const std::string tagIon = tagAtom + "+";
		assert(tagElectron == "e");
		std::cout << "Inelastic collision: " << tagElectron << ", " << tagAtom << ", " << tagIon << std::endl;

		const int idxE = getSpeciesIndex(tagElectron);
		const int idxA = getSpeciesIndex(tagAtom);
		const int idxI = getSpeciesIndex(tagIon);

		const std::string folderPath = "collisions/inelastic/";
		InelasticCollision * inelasticCollision = new InelasticCollision(folderPath, idxE, idxA, idxI, scalingParameters);
		const double electronMassRatio = getSpecies(idxE).getMassRatio();
		inelasticCollisions.push_back(inelasticCollision);
	}
	std::cout << "Number of inelastic collisions: " << inelasticCollisions.size() << std::endl;
}

void AppmSolver::setScalingParameters(const std::string & filename)
{
	// Read file with scaling parameters
	scalingParameters = ScalingParameters(filename);
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
	return speciesList.size();
}

const Eigen::VectorXd AppmSolver::getState(const int cIdx, const int fluidx) const
{
	assert(cIdx >= 0 && cIdx < dualMesh.getNumberOfCells());
	assert(fluidx >= 0 && fluidx < getNFluids());
	return fluidStates.col(cIdx).segment(5 * fluidx, 5);
}

//std::string AppmSolver::printSolverParameters() const
//{
//	std::stringstream ss;
//	ss << "Appm Solver parameters:" << std::endl;
//	ss << "=======================" << std::endl;
//	ss << "maxIterations:  " << appmParams.maxIterations << std::endl;
//	ss << "maxTime:        " << appmParams.maxTime << std::endl;
//	ss << "timestepSize:   " << appmParams.timestepSize << std::endl;
//	ss << "isFluidEnabled: " << appmParams.isFluidEnabled << std::endl;
//	ss << "isMaxwellEnabled: " << appmParams.isMaxwellEnabled << std::endl;
//	ss << "MaxwellSolverType: " << appmParams.maxwellSolverType << std::endl;
//	ss << "lambdaSquare: " << lambdaSquare << std::endl;
//	ss << "massFluxScheme: " << appmParams.isMassFluxSchemeImplicit << std::endl;
//	ss << "initType: " << initType << std::endl;
//	ss << "isElectricLorentzForceActive: " << appmParams.isLorentzForceElectricEnabled << std::endl;
//	ss << "isMagneticLorentzForceActive: " << appmParams.isLorentzForceMagneticEnabled << std::endl;
//	ss << "isShowDataWriterOutput: " << isShowDataWriterOutput << std::endl;
//	ss << "isMaxwellCurrentDefined: " << appmParams.isMaxwellCurrentDefined << std::endl;
//	ss << "isEulerMaxwellCouplingEnabled: " << appmParams.isEulerMaxwellCouplingEnabled << std::endl;
//	ss << "isFrictionActive: " << appmParams.isFrictionActive << std::endl;
//	ss << "=======================" << std::endl;
//	return ss.str();
//}

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

	//Eigen::sparseMatrixToFile(Q, "Q.dat");
	H5Writer h5writer("discreteMaps.h5");
	Eigen::sparseMatrixToFile(Q, "/Q", h5writer);

	// Material laws
	this->Meps = getElectricPermittivityOperator();
	assert(Meps.nonZeros() > 0);
	//Eigen::sparseMatrixToFile(Meps, "Meps.dat");
	Eigen::sparseMatrixToFile(Meps, "/Meps", h5writer);

	this->Mnu = getMagneticPermeabilityOperator();
	assert(Mnu.nonZeros() > 0);
	//Eigen::sparseMatrixToFile(Mnu, "Mnu.dat");
	Eigen::sparseMatrixToFile(Mnu, "/Mnu", h5writer);

	// Curl operator on all primal faces and all primal edges
	this->C = primalMesh.get_f2eMap().cast<double>();
	assert(C.nonZeros() > 0);

	// Curl operator on inner primal faces and inner primal edges
	const Eigen::SparseMatrix<double> C_inner = C.topLeftCorner(nPfi, nPei);
	//Eigen::sparseMatrixToFile(C_inner, "Ci.dat");
	Eigen::sparseMatrixToFile(C_inner, "/Ci_", h5writer);

	Eigen::SparseMatrix<double> P(nPfi, nDof);
	assert(P.cols() == nPei + nPvb);
	P.leftCols(nPei) = C_inner;
	//Eigen::sparseMatrixToFile(P, "P.dat");
	Eigen::sparseMatrixToFile(P, "/P", h5writer);

	const double lambdaSq = solverParams.getApParameter();
	M1 = lambdaSq * Q.transpose() * Meps * Q;
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
	Jaux_cc = Eigen::Matrix3Xd::Zero(3, dualMesh.getNumberOfCells());


	// Initialize matrix that maps electric field on primal edges to electric current on dual faces; see Lorentz force in Fluid equations
	//initMsigma();


	const int nFaces = dualMesh.getNumberOfFaces();
	J_h = Eigen::VectorXd::Zero(nFaces);
	J_h_previous = J_h;

	J_h_aux = Eigen::VectorXd(nFaces);
	J_h_aux.setZero();

	J_h_aux_mm1 = Eigen::VectorXd(nFaces);
	J_h_aux_mm1.setZero();

	// Initialize electric field 
	const Eigen::Vector3d initEfield = solverParams.getInitEfield();
	assert(initEfield.allFinite());
	const int nEdges = primalMesh.getNumberOfEdges();
	for (int i = 0; i < nEdges; i++) {
		const Edge * edge = primalMesh.getEdge(i);
		const Eigen::Vector3d edgeDir = edge->getDirection();
		E_h(i) = edgeDir.dot(initEfield);
	}
	E_cc = getEfieldAtCellCenter();

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

void AppmSolver::init_multiFluid()
{
	assert(speciesList.size() > 0);
	const int nFluids = speciesList.size();

	// Initialize data
	const int nCells = dualMesh.getNumberOfCells();
	const int nFaces = dualMesh.getNumberOfFaces();
	const int fluidStateLength = getFluidStateLength();

	fluidStates = Eigen::MatrixXd::Zero(fluidStateLength, nCells);
	fluidSources = Eigen::MatrixXd::Zero(fluidStateLength, nCells);
	sumOfFaceFluxes = Eigen::MatrixXd::Zero(fluidStates.rows(), nCells);
	frictionForceSourceTerm = Eigen::MatrixXd::Zero(3 * nFluids, nCells);
	frictionEnergySourceTerm = Eigen::MatrixXd::Zero(nFluids, nCells);
	diffusionVelocity = Eigen::MatrixXd::Zero(3 * nFluids, nCells);
	bulkVelocity = Eigen::Matrix3Xd::Zero(3, nCells);
	massFluxImplicitTerm = Eigen::MatrixXd::Zero(nFluids, nFaces);

	faceTypeFluids = Eigen::MatrixXi::Zero(nFaces, nFluids);
	for (int fluidIdx = 0; fluidIdx < nFluids; fluidIdx++) {
		for (int i = 0; i < nFaces; i++) {
			const Face * face = dualMesh.getFace(i);
			Face::Type faceType = getFaceTypeOfFluid(face, fluidIdx);
			faceTypeFluids(i, fluidIdx) = static_cast<int>(faceType);
		}
	}

	applyFluidInitializationType();
}

/**
* Define initial state of fluid.
*/
void AppmSolver::applyFluidInitializationType()
{
	// set default value for all fluid states
	fluidStates.setConstant(std::nan(""));

	const FluidInitType initType = solverParams.getFluidInitType();
	switch (initType) {
	case FluidInitType::SHOCKTUBE:
	{
		std::cout << "Initialize fluid states: " << "Shock tube" << std::endl;
		const double zRef = 0.;
		init_SodShockTube(zRef);
	}
	break;

	case FluidInitType::UNIFORM:
	{
		std::cout << "Initialize fluid states: " << "Uniform" << std::endl;
		double p = 1.0;
		double n = 1.0;
		Eigen::Vector3d uvec = Eigen::Vector3d::Zero();
		init_Uniformly(n, p, uvec);
	}
	break;

	case FluidInitType::TEST_FRICTION:
	{
		std::cout << "Initialize fluid states: " << "Test Friction" << std::endl;
		init_fluid_frictonTest();
	}
	break;

	case FluidInitType::TEST_FRICTION_TEMPERATURE:
	{
		std::cout << "Initialize fluid states: " << initType << std::endl;
		init_fluid_frictionTest_temperature();
	}
	break;

	case FluidInitType::TEST_FRICTION_NUMBERDENSITY:
	{
		std::cout << "Initialize fluid states: " << initType << std::endl;
		init_fluid_frictionTest_numberDensity();
	}
	break;

	case FluidInitType::TEST_FRICTION_ELECTRONS_NONZERO_VELOCITY:
	{
		std::cout << "Initialize fluid states: " << initType << std::endl;
		init_fluid_frictionTest_electrons_nonzero_velocity();
	}
	break;

	case FluidInitType::INIT_FILE:
	{
		const std::string filename = "initFluid.txt";
		std::cout << "Initialize from init file: " << filename << std::endl;
		init_multiFluid_readFromFile(filename);
	}
	break;


	//case FluidInitType::DEFAULT:
	//{
	//	std::cout << "Initialize fluid states: " << "Explosion" << std::endl;
	//	const Eigen::Vector3d refPos = Eigen::Vector3d(0, 0, 0);
	//	const double radius = 0.2;
	//	init_Explosion(refPos, radius);
	//}
	//break;

	//case 4:
	//{
	//	init_testcase_frictionTerm();
	//}
	//break;

	//case 5:
	//{
	//	init_ignitionWire();
	//}
	//break;

	default:
		std::cout << "FluidInitType not implemented: " << initType << std::endl;
		exit(-1);
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
		const double epsilon2 = getSpecies(k).getMassRatio();

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
		const double epsilon2 = getSpecies(k).getMassRatio();
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
		const double massRatio = getSpecies(fluidIdx).getMassRatio();
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

/**
* Initialize fluids such that first fluid is at velocity u = 1 and all others are quiescent (u = 0).
*/
void AppmSolver::init_testcase_frictionTerm()
{
	const int nCells = dualMesh.getNumberOfCells();
	const int nFluids = getNFluids();
	assert(nFluids >= 2);
	for (int fluidIdx = 0; fluidIdx < nFluids; fluidIdx++) {
		const double massRatio = getSpecies(fluidIdx).getMassRatio();
		const double n = 1;
		const double p = 1;
		Eigen::Vector3d uvec;
		Eigen::VectorXd state;
		if (fluidIdx == 0) {
			uvec = Eigen::Vector3d::UnitZ();
			state = Physics::primitive2state(massRatio, n, p, uvec);
		}
		else {
			uvec = Eigen::Vector3d::Zero();
			state = Physics::primitive2state(massRatio, n, p, uvec);
		}
		for (int i = 0; i < nCells; i++) {
			fluidStates.col(i).segment(5 * fluidIdx, 5) = state;
		}
	}
}

/**
* Initialize fluid states like an ignition wire, i.e., a quasineutral, hot, thin channel of gas.
*/
void AppmSolver::init_ignitionWire()
{
	const int nCells = dualMesh.getNumberOfCells();
	for (int cellIdx = 0; cellIdx < nCells; cellIdx++) {
		const Cell * cell = dualMesh.getCell(cellIdx);
		if (cell->getType() != Cell::Type::FLUID) {
			// Skip non-fluid cells
			continue;
		}
		const Eigen::Vector3d & pos = cell->getCenter();
		const Eigen::Vector2d pos2d = pos.segment(0, 2);
		for (int fluidx = 0; fluidx < getNFluids(); fluidx++) {
			const double massRatio = getSpecies(fluidx).getMassRatio();
			const int q = getSpecies(fluidx).getCharge();
			const Eigen::Vector3d u = Eigen::Vector3d::Zero();
			double n = 1;
			double p = 1;
			Eigen::VectorXd state;
			if (q == 0) {
				n = 1;
				p = 1;
			}
			else {
				if (pos2d.norm() < 0.1) {
					n = 1;
					p = 2;
				}
				else {
					n = 0.1;
					p = 1;
				}
			}
			state = Physics::primitive2state(massRatio, n, p, u);
			fluidStates.col(cellIdx).segment(5 * fluidx, 5) = state;
		}
	}
}

/**
* Initialize all fluids with zero velocity, but the first charge-neutral fluid with velocity u = 1.
*/
void AppmSolver::init_fluid_frictonTest()
{
	const int nCells = dualMesh.getNumberOfCells();
	const int nFluids = getNFluids();

	// index of first fluid that is charge-neutral
	int idxN = -1;
	for (int fluidx = 0; fluidx < nFluids; fluidx++) {
		if (getSpecies(fluidx).getCharge() == 0) {
			idxN = fluidx;
			break;
		}
	}
	assert(idxN >= 0);
	std::cout << "Fluid with nonzero velocity: " << getSpecies(idxN).getName() << std::endl;

	const double n = 1;
	const double p = 1;
	const Eigen::Vector3d uZero = Eigen::Vector3d::Zero();
	const Eigen::Vector3d uNonZero = Eigen::Vector3d::UnitZ();

	for (int fluidx = 0; fluidx < nFluids; fluidx++) {
		const double massRatio = getSpecies(fluidx).getMassRatio();
		Eigen::VectorXd stateZero = Physics::primitive2state(massRatio, n, p, uZero);
		Eigen::VectorXd stateNonZero = Physics::primitive2state(massRatio, n, p, uNonZero);

		for (int i = 0; i < nCells; i++) {
			if (fluidx == idxN) {
				fluidStates.col(i).segment(5 * fluidx, 5) = stateNonZero;
			}
			else {
				fluidStates.col(i).segment(5 * fluidx, 5) = stateZero;
			}
		}
	}

	std::cout << "Number of elastic collisions: " << elasticCollisions.size() << std::endl;
}

void AppmSolver::init_fluid_frictionTest_temperature()
{
	const int nCells = dualMesh.getNumberOfCells();
	const int nFluids = getNFluids();

	// index of first fluid that is charge-neutral
	int idxN = -1;
	for (int fluidx = 0; fluidx < nFluids; fluidx++) {
		if (getSpecies(fluidx).getCharge() == 0) {
			idxN = fluidx;
			break;
		}
	}
	assert(idxN >= 0);

	const Eigen::Vector3d uZero = Eigen::Vector3d::Zero();
	const double n = 1;
	const double p0 = 1;
	const double p1 = 2;

	for (int fluidx = 0; fluidx < nFluids; fluidx++) {
		const double massRatio = getSpecies(fluidx).getMassRatio();
		Eigen::VectorXd stateZero = Physics::primitive2state(massRatio, n, p0, uZero);
		Eigen::VectorXd stateNonZero = Physics::primitive2state(massRatio, n, p1, uZero);

		for (int i = 0; i < nCells; i++) {
			if (fluidx == idxN) {
				fluidStates.col(i).segment(5 * fluidx, 5) = stateNonZero;
			}
			else {
				fluidStates.col(i).segment(5 * fluidx, 5) = stateZero;
			}
		}
	}
}

/**
* Initialize charge-neutral fluid
*/
void AppmSolver::init_fluid_frictionTest_numberDensity()
{
	const int nCells = dualMesh.getNumberOfCells();
	const int nFluids = getNFluids();
	const int nFluidCells = dualMesh.getNumberFluidCells();

	const double n0 = 1;
	const double n1 = 2;
	const double p0 = 1;
	const Eigen::Vector3d u0 = Eigen::Vector3d::Zero();

	// assume that we have three fluids: neutrals, electrons (negative), ions (positive)
	assert(nFluids == 3);

	// find index of neutral species, electrons, and ions
	int idxN = -1;
	int idxE = -1;
	int idxI = -1;
	assert(nFluids >= 3);
	for (int fluidx = 0; fluidx < nFluids; fluidx++) {
		if (idxN == -1 && getSpecies(fluidx).getCharge() == 0) {
			idxN = fluidx;
		}
		if (idxE == -1 && getSpecies(fluidx).getCharge() == -1) {
			idxE = fluidx;
		}
		if (idxI == -1 && getSpecies(fluidx).getCharge() == +1) {
			idxI = fluidx;
		}
	}
	assert(idxN >= 0);
	assert(idxE >= 0);
	assert(idxI >= 0);
	
	// define quiescent, charge-neutral plasma with number density n > n0 for neutral species
	for (int fluidx = 0; fluidx < nFluids; fluidx++) {
		const double massRatio = getSpecies(fluidx).getMassRatio();
		const Eigen::VectorXd state0 = Physics::primitive2state(massRatio, n0, p0, u0);
		const Eigen::VectorXd state1 = Physics::primitive2state(massRatio, n1, p0, u0);
		const Eigen::VectorXd state = (fluidx != idxN) ? state0 : state1;
		for (int k = 0; k < nFluidCells; k++) {
			fluidStates.col(k).segment(5 * fluidx, 5) = state;
		}
	}
}

void AppmSolver::init_fluid_frictionTest_electrons_nonzero_velocity()
{
	const int nCells = dualMesh.getNumberOfCells();
	const int nFluids = getNFluids();
	const int nFluidCells = dualMesh.getNumberFluidCells();

	const double n0 = 1;
	const double p0 = 1;
	const Eigen::Vector3d u0 = Eigen::Vector3d::Zero();
	const Eigen::Vector3d u1 = Eigen::Vector3d::UnitZ();

	// assume that we have three fluids: neutrals, electrons (negative), ions (positive)
	assert(nFluids == 3);

	// find index of neutral species, electrons, and ions
	int idxN = -1;
	int idxE = -1;
	int idxI = -1;
	assert(nFluids >= 3);
	for (int fluidx = 0; fluidx < nFluids; fluidx++) {
		if (idxN == -1 && getSpecies(fluidx).getCharge() == 0) {
			idxN = fluidx;
		}
		if (idxE == -1 && getSpecies(fluidx).getCharge() == -1) {
			idxE = fluidx;
		}
		if (idxI == -1 && getSpecies(fluidx).getCharge() == +1) {
			idxI = fluidx;
		}
	}
	assert(idxN >= 0);
	assert(idxE >= 0);
	assert(idxI >= 0);

	// define quiescent, charge-neutral plasma with nonzero electron velocity
	for (int fluidx = 0; fluidx < nFluids; fluidx++) {
		const double massRatio = getSpecies(fluidx).getMassRatio();
		const Eigen::VectorXd state0 = Physics::primitive2state(massRatio, n0, p0, u0);
		const Eigen::VectorXd stateNonZero = Physics::primitive2state(massRatio, n0, p0, u1);
		const Eigen::VectorXd state = (fluidx != idxE) ? state0 : stateNonZero;
		for (int k = 0; k < nFluidCells; k++) {
			fluidStates.col(k).segment(5 * fluidx, 5) = state;
		}
	}

}

/**
* Read initial conditions for multi-fluid from a text file. 
* Assume that it is structured as follows:
* - Fluids are identified by their name (the name is given in the input file)
* - Initial state is given by number density (n), temperature (T), and velocity z-component (uz)
* - Format example: 
*     Ion
*     n: 1
*     T: 1
*     uz: 0
* - Empty lines and comment lines (start with #) are skipped
*/
void AppmSolver::init_multiFluid_readFromFile(const std::string & filename)
{
	std::vector<InitDataStruct> initDataVec;
	InitDataStruct initData;


	// Read input file and store data in a struct
	std::ifstream file(filename);
	std::string line;
	while (std::getline(file, line)) {
		if (line.empty()) { 
			continue; // skip empty lines
		}
		if (line.front() == '#') {
			continue; // skip comment lines
		}
		std::stringstream ss(line);
		std::string varName;
		double value;
		ss >> varName >> value;

		if (varName.size() == 2) {
			std::string temp = varName.substr(0,2);
			if (temp == "n:") {
				initData.n = value;
			}
			if (temp == "T:") {
				initData.T = value;
			}
			if (temp == "u:") {
				initData.uz = value;
			}
		}
		else {
			// If fluidName is already defined
			if (initData.fluidName.size() > 0) {
				// push existing init data to vector of initData 
				assert(initData.n > 0 && initData.T > 0 && initData.fluidName.size() > 0);
				initDataVec.push_back(initData);

				// create an empty struct for next fluid
				initData = InitDataStruct();
			}
			initData.fluidName = varName.substr(0, varName.size() - 1);
		}
	}
	// Add last fluid init data
	assert(initData.n > 0 && initData.T > 0 && initData.fluidName.size() > 0);
	initDataVec.push_back(initData);

	// Show init data that has been read
	std::cout << "Number of fluid init data read: " << initDataVec.size() << std::endl;
	std::cout << "**************" << std::endl;
	for (int i = 0; i < initDataVec.size(); i++) {
		std::cout << initDataVec[i] << std::endl;
		std::cout << "**************" << std::endl;
	}

	// Apply data to fluids
	for (int i = 0; i < initDataVec.size(); i++) {
		// find species index 
		int fidx = -1;
		double massRatio = 0;
		for (int j = 0; j < getNFluids(); j++) {
			if (speciesList[j].getName() == initDataVec[i].fluidName) {
				fidx = j;
				massRatio = speciesList[j].getMassRatio();
			}
		}
		assert(fidx >= 0);
		assert(fidx < getNFluids());

		// Define fluid state from init data
		double n = initDataVec[fidx].n;
		double T = initDataVec[fidx].T;
		double uz = initDataVec[fidx].uz;
		Eigen::Vector3d u(0, 0, uz);
		double p = n * T;
		Eigen::VectorXd state = Physics::primitive2state(massRatio, n, p, u);

		// check if state is computed consistently:
		{
			double Ts = Physics::getTemperature(state, massRatio);
			std::cout << "Check temperature conversion: " << std::endl;
			std::cout << "T  = " << T << std::endl;
			std::cout << "Ts = " << Ts << std::endl;
			const double relErr = std::abs(Ts / T - 1.); // relative error between T and Ts
			assert(relErr < 1e-4); // check if Ts and T are almost equal
		}

		// Set fluid state for all fluid cells
		for (int k = 0; k < dualMesh.getNumberFluidCells(); k++) {
			fluidStates.col(k).segment(5 * fidx, 5) = state;
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
		//B_h(idx) = thetaVec.dot(fn) * fA;  // Projection of azimuthal vector onto face normal
		B_h(idx) = thetaVec.dot(fn);  // Projection of azimuthal vector onto face normal
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
			const Eigen::Vector3d n = face->getNormal().normalized();
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

	// Test: guard against trunctation error
	//const double EmaxCoeff = E_h.cwiseAbs().maxCoeff();
	//const double scale = 1e6;
	//E_h.array() += scale * EmaxCoeff;
	//E_h.array() -= scale * EmaxCoeff;

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

	assert(J_h_aux.allFinite());
	Eigen::VectorXd Jaux_cc_vectorFormat = M * J_h_aux;
	Jaux_cc = Eigen::Map<Eigen::Matrix3Xd>(Jaux_cc_vectorFormat.data(), 3, nCells);

	Eigen::VectorXd Jcc_vectorFormat = M * J_h;
	Eigen::Map<Eigen::Matrix3Xd> result(Jcc_vectorFormat.data(), 3, nCells);

	return result;
}

/**
* Set source term in energy equation due to radiation.
*/
void AppmSolver::setRadiationSource()
{
	assert(false);
	// Extend implementation before using it:
	// - apply to electron fluid
	// - include and interpolate from material data 

	int nFluids = getNFluids();
	assert(nFluids >= 1);
	nFluids = std::max(1, nFluids);	 // Consider only the first fluid species
	const int nCells = dualMesh.getNumberOfCells();
	for (int i = 0; i < nCells; i++) {
		const Cell * cell = dualMesh.getCell(i);
		if (cell->getType() != Cell::Type::FLUID) { continue; }
		for (int fluidx = 0; fluidx < nFluids; fluidx++) {
			const Eigen::VectorXd state = fluidStates.col(i).segment(5 * fluidx, 5);
			assert(state.size() == 5);
			const double massRatio = getSpecies(fluidx).getMassRatio();
			double T = 0;
			T = Physics::getTemperature(state, massRatio);
			const double Qrad = -pow(T, 4);	// use T^4-law
			fluidSources(5 * 0 + 4, i) += Qrad;
		}
	}
}

/**
* Apply friction force to all fluids. 
*
* Momentum source for species a: R_a = (-1) * sum_b( n_a * n_b * (u_a - u_b) * z_ab ), 
* where b stands for all other fluids (b != a), and z_ab is a function that describes further parameters 
* of the interaction between species a and b, e.g. temperature, mass ratio, and collision cross section.
*
* This also results in an energy source term Q_a = u_a * R_a.
*/
//void AppmSolver::setFrictionSourceTerms()
//{
//	auto nCells = dualMesh.getNumberOfCells();
//	auto nFluids = getNFluids();
//
//	if (solverParams.getFluidEnabled()) {
//		// For each fluid cell ...
//		for (int i = 0; i < nCells; i++) {
//			const Cell * cell = dualMesh.getCell(i);
//			if (cell->getType() != Cell::Type::FLUID) {
//				continue; // Skip non-fluid cells
//			}
//
//			double rho_avg = 0;    // bulk density
//			Eigen::Vector3d u_avg = Eigen::Vector3d::Zero(); // bulk velocity (mass-averaged velocity)
//
//			for (int alpha = 0; alpha < nFluids; alpha++) {
//				auto massRatio = getSpecies(alpha).getMassRatio();
//				auto state = fluidStates.col(i).segment(5 * alpha, 5);
//				assert(state.size() == 5);
//				auto n_u = state.segment(1, 3);
//				u_avg += massRatio * n_u;
//				rho_avg += massRatio * state(0);
//			}
//			assert(rho_avg > 0);
//			u_avg /= rho_avg;
//			bulkVelocity.col(i) = u_avg;
//
//			// For each species alpha ...
//			for (int alpha = 0; alpha < nFluids; alpha++) {
//				Eigen::Vector3d R_a = Eigen::Vector3d::Zero();
//				auto stateA = fluidStates.col(i).segment(5 * alpha, 5); // state vector of fluid A in cell i
//				auto n_a = stateA(0); // number density
//				auto u_a = stateA.segment(1, 3) / n_a; // velocity vector
//
//				// ... for each species beta != alpha ...
//				for (int beta = 0; beta < nFluids; beta++) {
//					auto stateB = fluidStates.col(i).segment(5 * beta, 5);
//					auto n_b = stateB(0); // number density
//					auto u_b = stateB.segment(1, 3) / n_b; // velocity vector
//
//					if (beta == alpha) {
//						continue; // skip if indices are equal, since they contribute nothing
//					}
//					// species interaction term
//					//auto z_ab = 1; // assume a constant value for testing purpose
//					// R_a -= n_a * n_b * (u_a - u_b) * z_ab;
//
//					// ... get friction force due to interaction of species Alpha and Beta
//					auto m_ab = getReducedMass(alpha, beta);
//					auto nu_ab = getCollisionFrequency(alpha, beta, i);
//					R_a -= n_a * m_ab * nu_ab * (u_a - u_b);
//				} // end for each species B
//				auto w_a = u_a - u_avg; // Diffusion velocity w_a 
//				auto Q_a = u_a.dot(R_a); // Source term for energy equation due to friction
//
//				// Save local data for post-processing
//				diffusionVelocity.col(i).segment(3 * alpha, 3) = w_a;
//				if (solverParams.getFrictionActive()) {
//					frictionForceSourceTerm.col(i).segment(3 * alpha, 3) = R_a;
//					frictionEnergySourceTerm(alpha, i) = Q_a;
//
//					fluidSources.col(i).segment(5 * alpha + 1, 3) += R_a; // set momentum source of fluid A
//					fluidSources(5 * alpha + 4, i) += Q_a;                // set   energy source of fluid A
//				} 
//			} // end for each species A
//		} // end for each fluid cell i
//	} // end if isFluidEnabled
//}

/**
* Evaluate magnetic Lorentz force (i.e., F_* = q_* * eps_*^(-2) * (nu)_* x B) 
* and add this source term to the fluid sources.
*/
void AppmSolver::setMagneticLorentzForceSourceTerms()
{
	const int nFluids = this->getNFluids();
	const int nCells = dualMesh.getNumberOfCells();

	for (int i = 0; i < nCells; i++) {
		const Cell * cell = dualMesh.getCell(i);
		if (cell->getType() != Cell::Type::FLUID) {
			continue; // Skip cell that are not of type Fluid 
		}
		for (int fluidIdx = 0; fluidIdx < nFluids; fluidIdx++) {
			assert(i < fluidStates.cols());
			assert(5 * fluidIdx + 3 < fluidStates.rows());
			const Eigen::Vector3d nu = fluidStates.col(i).segment(5 * fluidIdx + 1, 3);
			assert(nu.allFinite());
			const Eigen::Vector3d B = B_vertex.col(i);
			assert(B.allFinite());
			const int q = getSpecies(fluidIdx).getCharge();
			const double massRatio = getSpecies(fluidIdx).getMassRatio();
			const Eigen::Vector3d result = q * 1. / massRatio * nu.cross(B);
			if (!result.allFinite()) {
				std::cout << "result: " << result.transpose() << std::endl;
				assert(result.allFinite());
			}
			LorentzForce_magnetic.col(i).segment(3 * fluidIdx, 3) = result;
		}
	}
	//std::cout << "maxCoeff F_L magnetic: " << LorentzForce_magnetic.cwiseAbs().maxCoeff() << std::endl;

	// Update of momentum source term
	if (solverParams.getFluidEnabled() && solverParams.getLorentzForceEnabled()) {
		for (int fluidx = 0; fluidx < nFluids; fluidx++) {
			fluidSources.block(5 * fluidx + 1, 0, 3, nCells) += LorentzForce_magnetic.block(3 * fluidx, 0, 3, nCells);
		}
	}
}


/** 
* @return timestep size for next iteration, based on maximum local wavespeed 
*         and projected cell size at each face.
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
		const Eigen::Vector3d faceNormal = face->getNormal().normalized();
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
	double dt = dt_faces.minCoeff();
	assert(dt > 1e-12);

	dt *= solverParams.getTimestepSizeFactor();
	return dt;
}


/**
* @return Fluid state in given cell and fluid, projected in direction of 
*         a unit normal vector. The normal vector may point in any direction. 
*/
const Eigen::Vector3d AppmSolver::getFluidState(const int cellIdx, const int fluidIdx, const Eigen::Vector3d & faceNormal) const
{
	assert(cellIdx >= 0);
	assert(cellIdx < fluidStates.cols());
	assert(fluidIdx >= 0);
	assert(fluidIdx < this->getNFluids());
	Eigen::Vector3d fn = faceNormal.normalized();
	const double tol = 4 * std::numeric_limits<double>::epsilon();
	assert(abs(fn.norm() - 1) <= tol); // normal vector should be of unit length
	const Eigen::VectorXd state = fluidStates.col(cellIdx).segment(5 * fluidIdx, 5);
	return Eigen::Vector3d(state(0), state.segment(1,3).dot(fn), state(4));
}

/**
* Get orientation of cell and face as indicated by the face normal vector. 
* @return 1 if face normal has same orientation as the position vector of face center with respect to cell center; otherwise -1.
*/
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
* Get adjacient fluid cell states of a given face.
* @param face       (input) face between adjacient fluid cells.
* @param fluidIdx   (input) fluid index of which the states are retrieved.
* @param qL
* @param qR         (output) fluid states on the left- and right-side of face, as indicated by the face normal vector. 
* @return   cell index of the left- and right-adjacient cells.
*/
const std::pair<int,int> AppmSolver::getAdjacientCellStates(const Face * face, const int fluidIdx, Eigen::Vector3d & qL, Eigen::Vector3d & qR) const
{
	assert(face != nullptr);
	const std::vector<Cell*> faceCells = face->getCellList();
	assert(faceCells.size() >= 1);

	const Eigen::Vector3d faceNormal = face->getNormal().normalized();
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
			std::cout << "Boundary condition not implemented: " << std::endl;
			std::cout << "Face type: WALL" << std::endl;
			std::cout << "Cell types: SOLID and FLUID (presumably)" << std::endl;
			assert(false); // Not yet implemented
			exit(-1);
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

/**
* Create stop file with default value.
*/
void AppmSolver::createStopFile(const int value)
{
	std::ofstream(stopFilename) << value << std::endl;
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
	//std::cout << "face fluxes: " << std::endl;
	for (int fidx = 0; fidx < nFaces; fidx++) {
		const Face * face = dualMesh.getFace(fidx);
		// skip faces that have no adjacient fluid cell
		if (!face->hasFluidCells()) { continue; }
		assert(face->getType() != Face::Type::DEFAULT);

		const double fA = face->getArea();
		//const Eigen::Vector3d faceNormal = face->getNormal().normalized();

		// Face normal vector, with length equal to face area
		const Eigen::Vector3d faceNormal = face->getNormal();

		for (int fluidIdx = 0; fluidIdx < nFluids; fluidIdx++) {
			Eigen::VectorXd faceFlux3d(5);
			Eigen::Vector3d flux = getSpeciesFaceFlux(face, fluidIdx);
			//std::cout << "face idx: " << face->getIndex() << ", fluid " << fluidIdx << ": " << flux.transpose() << std::endl;
			faceFlux3d(0) = flux(0) * fA;
			faceFlux3d.segment(1, 3) = flux(1) * faceNormal;
			faceFlux3d(4) = flux(2) * fA;
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

	const int kRef = -1;
	for (int k = 0; k < nCells; k++) {
		const Cell * cell = dualMesh.getCell(k);
		// Skip non-fluid cells
		if (cell->getType() != Cell::Type::FLUID) { continue; }

		double cellVolume = cell->getVolume();

		// get sum of face fluxes
		const std::vector<Face*> cellFaces = cell->getFaceList();
		Eigen::VectorXd sumFluxes = Eigen::VectorXd::Zero(faceFluxes.rows()); 
		if (k == kRef) {
			std::cout << "Face fluxes: " << std::endl;
		}

		/*
		* The following loop is a summation over the face fluxes. 
		* Without further precautions, we will observe numerical cancellation effects 
		* because of adding and subtracting 'large' numbers (e.g., order of 1), and 
		* we will be left by truncation errors (e.g., order of 1e-16 = machine precision, 
		* times number of summands). 
		*
		* For guarding against such truncation errors, we also sum the absolute values 
		* and add/subtract it after the loop. This removes the truncation errors.
		*/

		// Accumulator for absolute values
		Eigen::VectorXd c = Eigen::VectorXd::Zero(sumFluxes.size());
		

		for (auto face : cellFaces) {
			//double faceArea = face->getArea();
			double orientation = getOrientation(cell, face);
			Eigen::VectorXd ff = faceFluxes.col(face->getIndex());
			//temp += orientation * ff * faceArea / cellVolume;
			Eigen::VectorXd temp = orientation * ff;
			if (k == kRef) {
				std::cout << temp.transpose() << std::endl;
			}
			c += temp.cwiseAbs();
			sumFluxes += temp;
		}
		// Guard against truncation errors
		c *= this->fluxTruncationErrorGuardScale; // scaling factor for effectively removing truncation errors in the low bits
		sumFluxes += c;
		sumFluxes -= c;

		if (k == kRef) {
			std::cout << "Sum of face fluxes: " << std::endl;
			std::cout << sumFluxes.transpose() << std::endl;
		}
		sumOfFaceFluxes.col(k) = sumFluxes / cellVolume;

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

		const double Ai = face->getArea();
		const Eigen::Vector3d ni = face->getNormal().normalized();

		// ... for each fluid ...
		for (int fluidx = 0; fluidx < nFluids; fluidx++) {
			const double numSchemeFactor = getNumericalSchemeFactor(face, fluidx);

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
					const Eigen::Vector3d nj = cellFace->getNormal().normalized();
					double ni_dot_nj = ni.dot(nj);
					if (abs(ni_dot_nj) < 1e-10) {
						ni_dot_nj = 0; // truncate small values
					}
					bool isValid = abs(ni_dot_nj) <= (1 + 2 * std::numeric_limits<double>::epsilon());
					if (!isValid) {
						std::cout << "ni_dot_nj: " << ni_dot_nj << std::endl;
						std::cout << std::endl;
					}
					assert(isValid);
					//const double fj = faceFluxes.col(j).segment(5 * fluidIdx + 1, 3).dot(nj);
					//faceSum += s_kj * fj * Aj * ni_dot_nj;
					const Eigen::Vector3d fj = faceFluxes.col(j).segment(5 * fluidx + 1, 3);
					faceSum += s_kj * fj.dot(ni) * Ai;
				}
				cellSum -= numSchemeFactor / Vk * faceSum;

				// Implicit term due to fluid momentum source 				
				Eigen::Vector3d Sk; // Momentum source term for this fluid
				Sk = fluidSources.col(cell->getIndex()).segment(5 * fluidx + 1, 3);
				cellSum += numSchemeFactor * Sk.dot(ni) * Ai;
			}
			cellSum *= dt;

			// Update mass flux with implicit terms
			massFluxImplicitTerm(fluidx, faceIdx) = cellSum;
			faceFluxes(5 * fluidx, faceIdx) += cellSum;
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
	const bool isElectronFluid = getSpecies(fluidIdx).getCharge() < 0;
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
	assert(false); // Check implementation before usage
	const double Ts = 3695; // melting point of Wolfram, Kelvin
	const double workFunction = 4.55; // units: eV
	const double j_em = Physics::thermionicEmissionCurrentDensity(Ts, workFunction); 

	// Check if face is at cathode terminal and we have an electron fluid
	const bool isElectronFluid = getSpecies(fluidIdx).getCharge() < 0;
	assert(isElectronFluid);
	Face::Type faceFluidType = face->getType();
	assert(faceFluidType == Face::Type::TERMINAL);

	// Set a default value for testing
	Eigen::Vector3d flux;
	flux(0) = 1;
	flux(1) = 2;
	flux(2) = 1;

	Eigen::Vector3d fn = face->getNormal().normalized();
	flux(1) *= fn.dot(Eigen::Vector3d::UnitZ());
	flux *= fn.dot(Eigen::Vector3d::UnitZ());
	return flux;
}


/** 
* Update to next timestep: U(m+1) = U(m) - dt / volume * sum(fluxes) + dt * source.
*/
void AppmSolver::updateFluidStates(const double dt, const bool isImplicitSources)
{
	const int nFluidCells = dualMesh.getNumberFluidCells();
	const int nFluids = getNFluids();
	const int n = 5 * nFluidCells * nFluids;

	fluidSources.setZero();
	Eigen::SparseMatrix<double> M_elastic;
	M_elastic = getJacobianEulerSourceElasticCollisions();
	assert(M_elastic.rows() == n);
	assert(M_elastic.cols() == n);

	Eigen::SparseMatrix<double> M_inelastic;
	Eigen::VectorXd rhs_inelastic(n);
	rhs_inelastic.setZero();

	M_inelastic = getJacobianEulerSourceInelasticCollisions(rhs_inelastic);
	assert(M_inelastic.rows() == n);
	assert(M_inelastic.cols() == n);

	Eigen::SparseMatrix<double> M;
	M = M_elastic + M_inelastic;

	if (isImplicitSources) {
		// implicit scheme
		Eigen::SparseMatrix<double> A(n, n);
		A.setIdentity();
		A -= dt * M;

		// data mapper from matrix to vector format
		Eigen::MatrixXd states = fluidStates.leftCols(nFluidCells);
		assert(states.size() == n);
		Eigen::Map<Eigen::VectorXd> statesVec(states.data(), states.size());
		Eigen::MatrixXd sumOfFluxes = sumOfFaceFluxes.leftCols(nFluidCells);
		assert(sumOfFluxes.size() == n);
		Eigen::Map<Eigen::VectorXd> sumOfFluxesVec(sumOfFluxes.data(), sumOfFluxes.size());

		// solve sparse system A*x = rhs
		assert(statesVec.size() == sumOfFluxesVec.size());
		Eigen::VectorXd rhs = statesVec - dt * sumOfFluxesVec - dt * rhs_inelastic;
		Eigen::VectorXd x(rhs.size());
		Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
		solver.compute(A);
		x = solver.solve(rhs);

		// Show output statistics of linear solver
		const int iters = solver.iterations();
		const double err = solver.error();
		const bool isConverged = solver.info() == Eigen::ComputationInfo::Success;
		std::cout << "Linear solver converged: " << isConverged << std::endl;
		std::cout << "iterations: " << iters << std::endl;
		std::cout << "error:      " << err << std::endl;
		assert(isConverged);

		// data mapper from vector format to matrix
		Eigen::Map<Eigen::MatrixXd> newStates(x.data(), fluidStates.rows(), nFluidCells);

		// set new fluid states
		fluidStates.leftCols(nFluidCells) = newStates;
	}
	else {
		/* Explicit scheme */

		// Get states of fluid cells in vector format
		Eigen::MatrixXd states = fluidStates.leftCols(nFluidCells);
		Eigen::Map<Eigen::VectorXd> statesVec(states.data(), states.size());

		// Get sources in vector format
		assert(M.cols() == statesVec.size());
		assert(statesVec.allFinite());
		Eigen::VectorXd srcVec = M * statesVec; 
		assert(srcVec.allFinite());

		// Reshape source vector to matrix format
		Eigen::Map<Eigen::MatrixXd> src_as_matrix(srcVec.data(), fluidSources.rows(), nFluidCells);
		fluidSources.leftCols(nFluidCells) = src_as_matrix;

		// Update states
		fluidStates.leftCols(nFluidCells) += 
			- dt * sumOfFaceFluxes.leftCols(nFluidCells) 
			+ dt * fluidSources.leftCols(nFluidCells);
	}
}

/**
* @return Jacobian matrix of Euler source terms due to elastic collisions.
*/
Eigen::SparseMatrix<double> AppmSolver::getJacobianEulerSourceElasticCollisions() const
{
	const bool isShowDebugMessage = false;
	const int nFluidCells = dualMesh.getNumberFluidCells();
	const int nFluids = getNFluids();

	const double gamma = Physics::gamma; 	// heat capacity ratio
	assert(gamma > 1);

	// Triplets for creating Jacobian matrix
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;

	const int nElasticCollisions = elasticCollisions.size();
	if (isShowDebugMessage) {
		std::cout << "Number of elastic collisions: " << nElasticCollisions << std::endl;
	}

	// For all elastic collisions ... 
	for (int collIdx = 0; collIdx < nElasticCollisions; collIdx++) {
		const ElasticCollision * collision = elasticCollisions[collIdx];
		// get fluid indices of species A and B
		const int fidxA = collision->getFluidIdxA();
		const int fidxB = collision->getFluidIdxB();
		assert(fidxA >= 0 && fidxA < getNFluids());
		assert(fidxB >= 0 && fidxB < getNFluids());

		// get states of fluid A and fluid B
		const Eigen::MatrixXd & statesA = getStates(fidxA, nFluidCells);
		const Eigen::MatrixXd & statesB = getStates(fidxB, nFluidCells);

		// get reduced mass ratio m_ab = m_a * m_b / (m_a + m_b)
		const double ma = getSpecies(fidxA).getMassRatio();
		const double mb = getSpecies(fidxB).getMassRatio();
		const double mab = getReducedMass(fidxA, fidxB);
		if (isShowDebugMessage) {
			std::cout << "mab: " << mab << std::endl;
		}

		// get reduced temperature T_ab = (ma * Tb + mb * Ta) / (ma + mb)
		const Eigen::VectorXd Ta = Physics::getTemperature(statesA, ma);
		const Eigen::VectorXd Tb = Physics::getTemperature(statesB, mb);
		const Eigen::VectorXd Tab = (1. / (ma + mb)) * (ma * Tb + mb * Ta);

		// get avg collision cross section Q_ab
		Eigen::VectorXd Qab = collision->getAvgMomCrossSection(Tab);

		// get thermal velocity v_th = sqrt(8/pi * T_ab / m_ab)
		const double sqrt_8_pi = M_2_SQRTPI * M_SQRT2; // sqrt(8/pi) = 2/sqrt(pi) * sqrt(2)
		Eigen::VectorXd v_th = (sqrt_8_pi / sqrt(mab)) * Tab.array().sqrt();

		const bool isDefaultDataUsed = true;
		if (isDefaultDataUsed) {
			std::cout << "Warning: use default data for elastic collisions" << std::endl;
			Qab.setConstant(1);
			v_th.setOnes();
		}

		// get number densities
		const Eigen::VectorXd & na = statesA.row(0);
		const Eigen::VectorXd & nb = statesB.row(0);
		assert(na.size() == statesA.cols() && na.size() > 0);
		assert(nb.size() == statesB.cols() && nb.size() > 0);

		// Define source terms in Jacobian
		int i, j;
		for (int k = 0; k < nFluidCells; k++) {
			// Gas velocity
			const Eigen::Vector3d ua = statesA.col(k).segment(1, 3) / na(k);
			const Eigen::Vector3d ub = statesB.col(k).segment(1, 3) / nb(k);

			// Factor for momentum and energy source
			const double factor_m = -4. / 3. * mab * v_th(k) * Qab(k);
			const double factor_e = factor_m / (ma + mb);

			// Index of fluid states of cell k in state vector
			const int idxA = 5 * (nFluids * k + fidxA);
			const int idxB = 5 * (nFluids * k + fidxB);

			/*
			* Momentum source
			*/
			const double faa = +factor_m / ma * nb(k);
			const double fab = -factor_m / ma * na(k);
			const double fba = -factor_m / mb * nb(k);
			const double fbb = +factor_m / mb * na(k);

			i = idxA;
			j = idxA;
			triplets.push_back(T(i + 1, j + 1, faa));
			triplets.push_back(T(i + 2, j + 2, faa));
			triplets.push_back(T(i + 3, j + 3, faa));

			i = idxA;
			j = idxB;
			triplets.push_back(T(i + 1, j + 1, fab));
			triplets.push_back(T(i + 2, j + 2, fab));
			triplets.push_back(T(i + 3, j + 3, fab));

			i = idxB;
			j = idxA;
			triplets.push_back(T(i + 1, j + 1, fba));
			triplets.push_back(T(i + 2, j + 2, fba));
			triplets.push_back(T(i + 3, j + 3, fba));

			i = idxB;
			j = idxB;
			triplets.push_back(T(i + 1, j + 1, fbb));
			triplets.push_back(T(i + 2, j + 2, fbb));
			triplets.push_back(T(i + 3, j + 3, fbb));

			/*
			* Kinetic energy source
			*/
			const Eigen::Vector3d uColl = ma * ua + mb * ub;
			const Eigen::Vector3d qMaa = +factor_e / ma * nb(k) * uColl;
			const Eigen::Vector3d qMab = -factor_e / ma * na(k) * uColl;
			const Eigen::Vector3d qMba = -factor_e / mb * nb(k) * uColl;
			const Eigen::Vector3d qMbb = +factor_e / mb * na(k) * uColl;

			i = idxA;
			j = idxA;
			triplets.push_back(T(i + 4, j + 1, qMaa(0)));
			triplets.push_back(T(i + 4, j + 2, qMaa(1)));
			triplets.push_back(T(i + 4, j + 3, qMaa(2)));

			i = idxA;
			j = idxB;
			triplets.push_back(T(i + 4, j + 1, qMab(0)));
			triplets.push_back(T(i + 4, j + 2, qMab(1)));
			triplets.push_back(T(i + 4, j + 3, qMab(2)));

			i = idxB;
			j = idxA;
			triplets.push_back(T(i + 4, j + 1, qMba(0)));
			triplets.push_back(T(i + 4, j + 2, qMba(1)));
			triplets.push_back(T(i + 4, j + 3, qMba(2)));

			i = idxB;
			j = idxB;
			triplets.push_back(T(i + 4, j + 1, qMbb(0)));
			triplets.push_back(T(i + 4, j + 2, qMbb(1)));
			triplets.push_back(T(i + 4, j + 3, qMbb(2)));

			/* 
			* Thermal energy source
			*/
			const double qTaa = +factor_e / ma * 3 * (gamma - 1) * ma * nb(k);
			const double qTab = -factor_e / ma * 3 * (gamma - 1) * mb * na(k);
			const double qTba = -factor_e / mb * 3 * (gamma - 1) * ma * nb(k);
			const double qTbb = +factor_e / mb * 3 * (gamma - 1) * mb * na(k);

			i = idxA;
			j = idxA;
			triplets.push_back(T(i + 4, j + 1, qTaa * (-0.5) * ua(0)));
			triplets.push_back(T(i + 4, j + 2, qTaa * (-0.5) * ua(1)));
			triplets.push_back(T(i + 4, j + 3, qTaa * (-0.5) * ua(2)));
			triplets.push_back(T(i + 4, j + 4, qTaa));

			i = idxA;
			j = idxB;
			triplets.push_back(T(i + 4, j + 1, qTab * (-0.5) * ub(0)));
			triplets.push_back(T(i + 4, j + 2, qTab * (-0.5) * ub(1)));
			triplets.push_back(T(i + 4, j + 3, qTab * (-0.5) * ub(2)));
			triplets.push_back(T(i + 4, j + 4, qTab));

			i = idxB;
			j = idxA;
			triplets.push_back(T(i + 4, j + 1, qTba * (-0.5) * ua(0)));
			triplets.push_back(T(i + 4, j + 2, qTba * (-0.5) * ua(1)));
			triplets.push_back(T(i + 4, j + 3, qTba * (-0.5) * ua(2)));
			triplets.push_back(T(i + 4, j + 4, qTba));

			i = idxB;
			j = idxB;
			triplets.push_back(T(i + 4, j + 1, qTbb * (-0.5) * ub(0)));
			triplets.push_back(T(i + 4, j + 2, qTbb * (-0.5) * ub(1)));
			triplets.push_back(T(i + 4, j + 3, qTbb * (-0.5) * ub(2)));
			triplets.push_back(T(i + 4, j + 4, qTbb));
		}
	}

	// Create Jacobian matrix from triplets
	const int n = 5 * nFluidCells * nFluids;
	Eigen::SparseMatrix<double> M(n, n);
	M.setFromTriplets(triplets.begin(), triplets.end());
	M.makeCompressed();
	return M;
}

/**
* Get linearized, implicit system of equations to solve for source terms S in fluid model 
* due to inelastic collisions. 
*
* The implicit system is given by A*x - b = S, where A is the Jacobian, 
* x is the vector of fluid states at new timestep m+1, and b is the right-hand-side vector with 
* data at timestep m (e.g., due to the linearization process).
*
* @param rhs   right-hand side vector
* @return Jacobian of inelastic collisions in fluid model.
*/
Eigen::SparseMatrix<double> AppmSolver::getJacobianEulerSourceInelasticCollisions(Eigen::VectorXd & rhs) const
{
	const int n = dualMesh.getNumberFluidCells();
	const int nFluids = getNFluids();
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;

	rhs.setZero();

	// Ratio of heat capacities for all fluids, assuming that all fluids are described as an ideal gas
	const double gamma = Physics::gamma;

	// Auxiliary vector with unit values
	const Eigen::VectorXd ones = Eigen::VectorXd::Ones(n);

	// Number of collisions
	const int nCollisions = inelasticCollisions.size(); 

	// For all collisions ...
	for (int collIdx = 0; collIdx < nCollisions; collIdx++) {
		const InelasticCollision * collision = inelasticCollisions[collIdx];

		// Ionization energy, dimensionless
		const double E_ion = collision->getIonizationEnergyScaled();

		// Fluid index for atoms, electrons, ions
		const int fidxA = collision->getAtomFluidx();
		const int fidxE = collision->getElectronFluidx();
		const int fidxI = collision->getIonFluidx();

		// Species mass ratios
		const double mA = getSpecies(fidxA).getMassRatio();
		const double mE = getSpecies(fidxE).getMassRatio();
		const double mI = getSpecies(fidxI).getMassRatio();

		// Species states
		const Eigen::MatrixXd & statesA = getStates(fidxA, n);
		const Eigen::MatrixXd & statesE = getStates(fidxE, n);
		const Eigen::MatrixXd & statesI = getStates(fidxI, n);

		// Species number density
		const Eigen::VectorXd nA = statesA.row(0);
		const Eigen::VectorXd nE = statesE.row(0);
		const Eigen::VectorXd nI = statesI.row(0);

		// Electron temperature
		const Eigen::VectorXd Ta = Physics::getTemperature(statesA, mA);
		const Eigen::VectorXd Te = Physics::getTemperature(statesE, mE);
		const Eigen::VectorXd Ti = Physics::getTemperature(statesI, mI);
		
		// Thermal velocity in electron fluid
		const Eigen::VectorXd vthE = ((8. / M_PI) * 1. / mE * Te).array().sqrt();

		// Fluid velocities 
		const Eigen::MatrixXd uAmat = Physics::getVelocity(statesA);
		const Eigen::MatrixXd uEmat = Physics::getVelocity(statesE);
		const Eigen::MatrixXd uImat = Physics::getVelocity(statesI);

		// Relative velocities in collision processes ...
		const Eigen::Matrix3Xd w0mat = uEmat - uAmat;  // ... in ionization process 
		const Eigen::Matrix3Xd w1mat = mI * (uEmat - uImat); // .. in recombination process

		// Ratio of kinetic energies to temperature in collision processes
		Eigen::VectorXd lambdaIon = Eigen::VectorXd::Zero(nA.size());
		Eigen::VectorXd lambdaRec = Eigen::VectorXd::Zero(nA.size());
		for (int i = 0; i < n; i++) {
			const Eigen::VectorXd w0 = w0mat.col(i);
			const Eigen::VectorXd w1 = w1mat.col(i);
			lambdaIon(i) = 0.5 * mE / Te(i) * w0.squaredNorm();
			lambdaRec(i) = 0.5 * mE / Te(i) * w1.squaredNorm();
		}
		// Ratio of ionization energy to temperature, dimensionless, for each fluid cell
		const Eigen::VectorXd xStar = E_ion * Te.array().inverse();

		// Inelastic collision coefficients as given by Le & Cambier (2016), with number densities set to 1
		const Eigen::VectorXd psi_Gion = collision->getGion(ones, ones, vthE, Te, lambdaIon);
		const Eigen::VectorXd psi_Grec = collision->getGrec(ones, ones, vthE, xStar, mE, Te, lambdaRec);
		const Eigen::VectorXd psi_R0ion = collision->getR0ion(ones, ones, vthE, Te, lambdaIon);
		const Eigen::VectorXd psi_R1rec = collision->getR1rec(ones, ones, vthE, xStar, Te, lambdaRec);
		const Eigen::VectorXd psi_R2rec = collision->getR2rec(ones, ones, vthE, xStar, Te, lambdaRec);
		const Eigen::VectorXd psi_J00ion = collision->getJ00ion(ones, ones, vthE, Te, lambdaIon);
		const Eigen::VectorXd psi_J11rec = collision->getJ11rec(ones, ones, vthE, xStar, Te, lambdaRec);
		const Eigen::VectorXd psi_J22rec = collision->getJ22rec(ones, ones, vthE, xStar, Te, lambdaRec);
		const Eigen::VectorXd psi_J12rec = collision->getJ12rec(ones, ones, vthE, xStar, Te, lambdaRec);

		// Coefficients in system of equations as given by Le & Cambier (2016)
		const Eigen::VectorXd psi_Rion = psi_R0ion;
		const Eigen::VectorXd psi_Rrec = psi_R1rec + psi_R2rec;
		const Eigen::VectorXd psi_Kion = psi_Gion - psi_R0ion;
		const Eigen::VectorXd psi_Krec = 2 * psi_Grec - psi_R1rec - psi_R2rec;
		const Eigen::VectorXd psi_Wion = psi_J00ion.array() - 2 * lambdaIon.array() * psi_R0ion.array() + lambdaIon.array() * psi_Gion.array();
		const Eigen::VectorXd psi_Wrec = psi_J11rec.array() + psi_J22rec.array() + 2 * psi_J12rec.array() + 4 * lambdaRec.array() * (psi_Gion.array() - psi_R1rec.array() - psi_R2rec.array());
		const Eigen::VectorXd psi_Jion = psi_J00ion.array() - lambdaIon.array() * psi_R0ion.array();
		const Eigen::VectorXd psi_Jrec = psi_J11rec.array() + psi_J22rec.array() + 2 * psi_J12rec.array() - 2 * lambdaRec.array() * (psi_R1rec.array() + psi_R2rec.array());

		// Define Jacobian matrix for inelastic sources.
		// For each fluid cell k ...
		for (int k = 0; k < n; k++) {
			const int idxE = getLinearIndexInJacobian(fidxE, k);
			const int idxI = getLinearIndexInJacobian(fidxI, k);
			const int idxN = getLinearIndexInJacobian(fidxA, k);

			// Species source
			{
				// Net species sources Gnet = Gion - Grec
				const double cNet_e = +psi_Gion(k) * nA(k) - 2 * psi_Grec(k) * nE(k) * nI(k); // factors in front of n_e^{m+1}
				const double cNet_i = -psi_Grec(k) * pow(nE(k), 2); // factors in front of n_i^{m+1}
				const double cNet_n = +psi_Gion(k) * nE(k); // factors in front of n_n^{m+1}
				const double cNet_rhs = -psi_Gion(k) * nE(k) * nA(k) + 2 * psi_Grec(k) * pow(nE(k), 2) * nI(k);

				triplets.push_back(T(idxE, idxE, +cNet_e)); // add Gnet to electron fluid
				triplets.push_back(T(idxE, idxN, +cNet_n)); // add Gnet to electron fluid
				triplets.push_back(T(idxE, idxI, +cNet_i)); // add Gnet to electron fluid
				
				triplets.push_back(T(idxI, idxE, +cNet_e)); // add Gnet to ion fluid
				triplets.push_back(T(idxI, idxN, +cNet_n)); // add Gnet to ion fluid
				triplets.push_back(T(idxI, idxI, +cNet_i)); // add Gnet to ion fluid

				triplets.push_back(T(idxN, idxE, -cNet_e)); // subtract Gnet from neutral fluid
				triplets.push_back(T(idxN, idxN, -cNet_n)); // subtract Gnet from neutral fluid
				triplets.push_back(T(idxN, idxI, -cNet_i)); // subtract Gnet from neutral fluid

				rhs(idxE) -= cNet_rhs; // add terms to rhs-vector
				rhs(idxI) -= cNet_rhs;
				rhs(idxN) += cNet_rhs;
			}
			// Momentum source
			{
				// (1 + mE) * (Gion*U0 - Grec*U1)
				const double c = 1. + mE;
				const double ce = +c * (psi_Gion(k) * nA(k) - psi_Grec(k) * 2 * nE(k) * nI(k)) / (1 + 1./mE); // factor in front of (nu)_e^{m+1}
				const double ci = -c * psi_Grec(k) * mI / (1. + mE) * pow(nE(k), 2); // factor in front of (nu)_i^{m+1}
				const double cn = +c * psi_Gion(k) * 1. / (1. + mE) * nE(k); // factor in front of (nu)_n^{m+1}
				for (int j = 1; j < 4; j++) {
					triplets.push_back(T(idxI + j, idxE + j, +mI * ce)); // add to ion fluid
					triplets.push_back(T(idxI + j, idxI + j, +mI * ci)); // add to ion fluid
					triplets.push_back(T(idxI + j, idxN + j, +mI * cn)); // add to ion fluid
					triplets.push_back(T(idxN + j, idxE + j, -ce)); // subtract from neutral fluid
					triplets.push_back(T(idxN + j, idxI + j, -ci)); // subtract from neutral fluid
					triplets.push_back(T(idxN + j, idxN + j, -cn)); // subtract from neutral fluid
				}
			}
			{
				// mE * Rion * w0 
				const double c = mE;
				const double ce = +c * psi_Rion(k) * nA(k); 
				const double cn = -c * psi_Rion(k) * nE(k);
				for (int j = 1; j < 4; j++) {
					triplets.push_back(T(idxE + j, idxE + j, -mE * ce)); // subtract from electron fluid
					triplets.push_back(T(idxE + j, idxN + j, -mE * cn)); // subtract from electron fluid
					triplets.push_back(T(idxN + j, idxE + j, +ce)); // add to neutral fluid
					triplets.push_back(T(idxN + j, idxN + j, -cn)); // add to neutral fluid
				}

			}
			{
				// mE * Rrec * w1
				const double c = mE;
				const double ce = +c * mI * psi_Rrec(k) * nE(k) * nI(k);
				const double ci = -c * psi_Rrec(k) * mI * pow(nE(k), 2);
				for (int j = 1; j < 4; j++) {
					triplets.push_back(T(idxE + j, idxE + j, -mE * ce)); // subtract from electron fluid
					triplets.push_back(T(idxE + j, idxI + j, -mE * ci)); // subtract from electron fluid
					triplets.push_back(T(idxI + j, idxE + j, +mI * ce)); // add to ion fluid
					triplets.push_back(T(idxI + j, idxI + j, +mI * ci)); // add to ion fluid
				}
			}
			{
				// mE * (Tn / Te - 1) * Kion * w0
				const double c = mE * psi_Kion(k) * (Ta(k) / Te(k) - 1.);
				const double ce = +c * nA(k);
				const double cn = -c * nE(k);
				for (int j = 1; j < 4; j++) {
					triplets.push_back(T(idxI + j, idxE + j, +mI * ce)); // add to ion fluid
					triplets.push_back(T(idxI + j, idxN + j, +mI * cn)); // add to ion fluid
					triplets.push_back(T(idxN + j, idxE + j, -ce)); // subtract from neutral fluid
					triplets.push_back(T(idxN + j, idxN + j, -cn)); // subtract from neutral fluid
				}
			}
			{
				// mE * (Ti / Te - 1) * Krec * w1
				const double c = mE * psi_Krec(k) * mI * (Ti(k) / Te(k) - 1.);
				const double ce = +c * nE(k) * nI(k);
				const double ci = -c * pow(nE(k), 2);
				for (int j = 1; j < 4; j++) {
					triplets.push_back(T(idxN + j, idxE + j, +ce)); // add to neutral fluid
					triplets.push_back(T(idxN + j, idxI + j, +ci)); // add to neutral fluid
					triplets.push_back(T(idxI + j, idxE + j, -mI * ce)); // subtract from ion fluid
					triplets.push_back(T(idxI + j, idxI + j, -mI * ci)); // subtract from ion fluid
				}
			}
			// Energy source
			const int idxEE = getLinearIndexInJacobian(fidxE, k) + 4; // row index for Energy source in Electron fluid
			const int idxEI = getLinearIndexInJacobian(fidxI, k) + 4; // row index for Energy sourcce in Ion fluid
			const int idxEA = getLinearIndexInJacobian(fidxA, k) + 4; // row index for Energy source in Atom (Neutral) fluid
			{
				// ionization energy source 
				// (Grec - Gion) * E_ion
				const double cNet_e = +E_ion * psi_Gion(k) * nA(k) - 2 * psi_Grec(k) * nE(k) * nI(k); // factors in front of n_e^{m+1}
				const double cNet_i = -E_ion * psi_Grec(k) * pow(nE(k), 2); // factors in front of n_i^{m+1}
				const double cNet_n = +E_ion * psi_Gion(k) * nE(k); // factors in front of n_n^{m+1}
				const double cNet_rhs = -E_ion * psi_Gion(k) * nE(k) * nA(k) + 2 * psi_Grec(k) * pow(nE(k), 2) * nI(k);
				triplets.push_back(T(idxEE, idxE, +mE * cNet_e)); // add to electron fluid
				triplets.push_back(T(idxEE, idxN, +mE * cNet_n)); // add to electron fluid
				triplets.push_back(T(idxEE, idxI, +mE * cNet_i)); // add to electron fluid
				
				rhs(idxE + 4) -= mE * cNet_rhs; // add terms to rhs vector
			}
			{
				// E_net = Gion*Eion - Grec*Erec
				{
					// kinetic energy term
					const double c = (1. + mE) / 2.;
					const Eigen::Vector3d ven = psi_Gion(k) * nA(k) / (1. + 1. / mE) * statesE.col(k).segment(1,3) + nE(k) / (1 + mE) * statesA.col(k).segment(1,3);
					const Eigen::Vector3d vei = psi_Grec(k) * 2 * nE(k) * nI(k) / (1. + 1./mE) * statesE.col(k).segment(1,3) + mI * pow(nE(k), 2) / (1. + mE) * statesI.col(k).segment(1,3);
					const Eigen::Vector3d ve = +c * (nA(k) / (1. + 1./mE) * ven - 2 * nE(k) * nI(k) / (1. + 1./mE) * vei); // vector in dot-product with (nu)_e^{m+1}
					const Eigen::Vector3d vi = -c * psi_Grec(k) * mI / (1. + mE) * pow(nE(k),2) * vei; // vector in dot-product with (nu)_i^{m+1}
					const Eigen::Vector3d vn = +c * psi_Gion(k) * nE(k) / (1 + mE) * ven; // vector in dot-product with (nu)_n^{m+1}
					for (int j = 1; j < 4; j++) {
						triplets.push_back(T(idxI + 4, idxE + j, +mI * ve(j - 1))); // add to ion fluid
						triplets.push_back(T(idxN + 4, idxE + j, -ve(j - 1))); // subtract from neutral fluid
						triplets.push_back(T(idxI + 4, idxI + j, +mI * vi(j - 1))); // add to ion fluid
						triplets.push_back(T(idxN + 4, idxI + j, -vi(j - 1))); // subtract from neutral fluid
						triplets.push_back(T(idxI + 4, idxN + j, +mI * vn(j - 1))); // add to ion fluid
						triplets.push_back(T(idxN + 4, idxN + j, -vn(j - 1))); // subtract from neutral fluid
					}
				}
				{
					// thermal energy term
					const double c_ion = +3. / 2. * psi_Gion(k) * nE(k) * (gamma - 1.) * mA;
					const double c_rec = -3. / 2. * psi_Grec(k) * pow(nE(k), 2) * (gamma - 1.) * mI;
					const Eigen::Vector3d vn = c_ion * -0.5 / nA(k) * statesA.col(k).segment(1, 3); // vector in dot-product with (nu)_n^{m+1}
					const Eigen::Vector3d vi = c_rec * -0.5 / nI(k) * statesI.col(k).segment(1, 3); // vector in dot-product with (nu)_i^{m+1}
					for (int j = 1; j < 5; j++) {
						if (j < 4) {
							triplets.push_back(T(idxI + 4, idxN + j, +mI * vn(j - 1))); // add to ion fluid
							triplets.push_back(T(idxI + 4, idxI + j, +mI * vi(j - 1))); // add to ion fluid
							triplets.push_back(T(idxN + 4, idxN + j, -vn(j - 1))); // subtract from neutral fluid
							triplets.push_back(T(idxN + 4, idxI + j, -vi(j - 1))); // subtract from neutral fluid
						}
						else {
							triplets.push_back(T(idxI + 4, idxN + j, +mI * c_ion)); // add to ion fluid
							triplets.push_back(T(idxI + 4, idxI + j, +mI * c_rec)); // add to ion fluid
							triplets.push_back(T(idxN + 4, idxN + j, -c_ion)); // subtract from neutral fluid
							triplets.push_back(T(idxN + 4, idxI + j, -c_rec)); // subtract from neutral fluid
						}
					}
				}
			}
			{
				// 1. / (1. + 1./mE) * (Tn - Te)^2 / Te * Wion
				const double c = 1. / (1. + 1. / mE) * psi_Wion(k) * (Ta(k) / Te(k) - 1.);
				const double ce = -c * nA(k) * (gamma - 1) * mE; // factor in front of implicit term (nT)_e^{m+1}
				const double cn = +c * nE(k) * (gamma - 1) * mA; // factor in front of implicit term (nT)_n^{m+1}
				for (int j = 1; j < 5; j++) {
					if (j < 4) {
						const double valueE = ce * -0.5 / nE(k) * statesE(j, k);
						const double valueA = cn * -0.5 / nA(k) * statesA(j, k);
						triplets.push_back(T(idxEA, idxE + j, +valueE)); // add to neutral fluid
						triplets.push_back(T(idxEI, idxE + j, -mI * valueE)); // subtract from ion fluid
						triplets.push_back(T(idxEA, idxN + j, +valueA)); // add to neutral fluid
						triplets.push_back(T(idxEI, idxN + j, -mI * valueA)); // subtract from ion fluid
					}
					else {
						triplets.push_back(T(idxEA, idxE + j, +ce)); // add to neutral fluid
						triplets.push_back(T(idxEI, idxE + j, -mI * ce)); // subtract from ion fluid
						triplets.push_back(T(idxEA, idxN + j, +cn)); // add to neutral fluid
						triplets.push_back(T(idxEI, idxN + j, -mI * cn)); // subtract from ion fluid
					}
				}
			}
			{
				// 1 / ( 1 + 1 / mE) * (Ti - Te)^2 / Te * Wrec
				const double c = 1. / (1 + 1 / mE) * psi_Wrec(k) * (Ti(k) / Te(k) - 1.);
				const double ce = -c * nE(k) * nI(k) * (gamma - 1) * mE; // factor in front of implicit term (nT)_e^{m+1}
				const double ci = +c * pow(nE(k), 2) * (gamma - 1) * mI; // factor in front of implicit term (nT)_i^{m+1}
				for (int j = 1; j < 5; j++) {
					if (j < 4) {
						const double valueE = ce * -0.5 / nE(k) * statesE(j,k);
						const double valueI = ci * -0.5 / nI(k) * statesI(j,k);
						triplets.push_back(T(idxEA, idxE+ j, +valueE)); // add to neutral fluid
						triplets.push_back(T(idxEI, idxE + j, -mI * valueE)); // subtract from ion fluid
						triplets.push_back(T(idxEA, idxI + j, +valueI)); // add to neutral fluid
						triplets.push_back(T(idxEI, idxI + j, -mI * valueI)); // subtract from ion fluid
					}
					else {
						triplets.push_back(T(idxEA, idxE + j, +ce)); // add to neutral fluid
						triplets.push_back(T(idxEI, idxE + j, -mI * ce)); // subtract from ion fluid
						triplets.push_back(T(idxEA, idxI + j, +ci)); // add to neutral fluid
						triplets.push_back(T(idxEI, idxI + j, -mI * ci)); // subtract from ion fluid
					}
				}
			}
			{
				// 2 / (1 + 1/mE) * (Tn - Te) * Jion
				const double c = 2. / (1. + 1. / mE) * psi_Jion(k);
				const double ce = -c * nA(k) * (gamma - 1) * mE; // factor in front of implicit term (nT)_e^{m+1}
				const double cn = +c * nE(k) * (gamma - 1) * mA; // factor in front of implicit term (nT)_n^{m+1}
				for (int j = 1; j < 5; j++) {
					if (j < 4) {
						const double valueE = ce * -0.5 / nE(k) * statesE(j, k);
						const double valueN = cn * -0.5 / nA(k) * statesA(j, k);
						triplets.push_back(T(idxEE, idxE + j, +mE * valueE)); // add to electron fluid
						triplets.push_back(T(idxEA, idxE + j, -valueE)); // subtract from neutral fluid
						triplets.push_back(T(idxEE, idxN + j, +mE * valueN)); // add to electron fluid
						triplets.push_back(T(idxEA, idxN + j, -valueN)); // subtract from neutral fluid
					}
					else {
						triplets.push_back(T(idxEE, idxE + j, +mE * ce)); // add to electron fluid
						triplets.push_back(T(idxEA, idxE + j, -ce)); // subtract from neutral fluid
						triplets.push_back(T(idxEE, idxN + j, +mE * cn)); // add to electron fluid
						triplets.push_back(T(idxEA, idxN + j, -cn)); // subtract from neutral fluid
					}
				}
			}
			{
				// 2 / (1 + 1/mE) * (Ti - Te) * Jrec
				const double c = 2. / (1. + 1. / mE) * psi_Jrec(k);
				const double ce = -c * nE(k) * nI(k) * (gamma - 1) * mE; // factor in front of implicit term (nT)_e^{m+1}
				const double ci = +c * pow(nE(k), 2) * (gamma - 1) * mI; // factor in front of implicit term (nT)_i^{m+1}
				for (int j = 1; j < 5; j++) {
					if (j < 4) {
						const double valueE = ce * -0.5 / nE(k) * statesE(j, k);
						const double valueI = ci * -0.5 / nI(k) * statesI(j, k);
						triplets.push_back(T(idxEE, idxE + j, +mE * valueE)); // add to electron fluid
						triplets.push_back(T(idxEI, idxE + j, -mI * valueE)); // subtract from ion fluid
						triplets.push_back(T(idxEE, idxI + j, +mE * valueI)); // add to electron fluid
						triplets.push_back(T(idxEI, idxI + j, -mI * valueI)); // subtract from ion fluid
					}
					else {
						triplets.push_back(T(idxEE, idxE + j, +mE * ce)); // add to electron fluid
						triplets.push_back(T(idxEI, idxE + j, -mI * ce)); // subtract from ion fluid
						triplets.push_back(T(idxEE, idxI + j, +mE * ci)); // add to electron fluid
						triplets.push_back(T(idxEI, idxI + j, -mI * ci)); // subtract from ion fluid
					}
				}
			}
			{
				// mE * (Tn / Te - 1) * Kion * w0 * U0
				const double c = mE * psi_Rion(k);
				const Eigen::Vector3d ve = 1. / (1. + 1. / mE) * 1. / nE(k) * statesE.col(k).segment(1, 3);
				const Eigen::Vector3d vn = 1. / (1. + mE) * 1. / nA(k) * statesA.col(k).segment(1, 3);
				const Eigen::Vector3d v = ve + vn; // vector in dot-product at timestep m

				const Eigen::Vector3d ce = +c * nA(k) * v; // factors in front of implicit term (nu)_e^{m+1}
				const Eigen::Vector3d cn = -c * nE(k) * v; // factors in front of implicit term (nu)_n^{m+1}
				for (int j = 1; j < 4; j++) {
					triplets.push_back(T(idxEA, idxE + j, +ce(j - 1))); // add to neutral fluid
					triplets.push_back(T(idxEA, idxN + j, +cn(j - 1))); // add to neutral fluid
					triplets.push_back(T(idxEI, idxE + j, -mI * ce(j - 1))); // subtract from ion fluid
					triplets.push_back(T(idxEI, idxN + j, -mI * cn(j - 1))); // subtract from ion fluid
				}
			}
			{
				// mE * (Ti/Te - 1) * Krec * w1 * U1
				const double c = mE * psi_Krec(k) * mI * (Ti(k) / Te(k) - 1.);
				const Eigen::Vector3d ve = 2. / (1. + 1./mE) * 1./ nE(k) * statesE.col(k).segment(1,3);
				const Eigen::Vector3d vi = mI / (1. + mE) * 1./nI(k) * statesI.col(k).segment(1,3);
				const Eigen::Vector3d v = ve + vi; // vector in dot-product at timestep m

				const Eigen::Vector3d ce = +c * nE(k) * nI(k) * v; // factor in front of (nu)_e^{m+1}
				const Eigen::Vector3d ci = -c * pow(nE(k), 2) * v; // factor in front of (nu)_i^{m+1}

				for (int j = 1; j < 4; j++) {
					triplets.push_back(T(idxEI, idxE + j, -mI * ce(j-1))); // subtract from ion fluid
					triplets.push_back(T(idxEI, idxI + j, -mI * ci(j-1))); // subtract from ion fluid
					triplets.push_back(T(idxEA, idxE + j, +ce(j-1))); // add to neutral fluid
					triplets.push_back(T(idxEA, idxI + j, +ci(j-1))); // add to neutral fluid
				}
			}
			{
				// mE * Rion * w0 * U0
				const double c = mE * psi_Rion(k);
				const Eigen::Vector3d ve = 1./(1. + 1./mE) * 1./nE(k) * statesE.col(k).segment(1,3); 
				const Eigen::Vector3d vn = 1./(1. + mE) * 1./nA(k) * statesA.col(k).segment(1,3);
				const Eigen::Vector3d v = ve + vn; // vector in dot-product at timestep m

				const Eigen::Vector3d ce = +c * nA(k) * v; // factor in dot-product with (nu)_e^{m+1}
				const Eigen::Vector3d cn = -c * nE(k) * v; // factor in dot-product with (nu)_n^{m+1}
				for (int j = 1; j < 4; j++) {
					triplets.push_back(T(idxEE, idxE + j, -mE * ce(j-1))); // subtract from electron fluid
					triplets.push_back(T(idxEE, idxN + j, -mE * cn(j-1))); // subtract from electron fluid
					triplets.push_back(T(idxEA, idxE + j, +ce(j-1))); // add to neutral fluid
					triplets.push_back(T(idxEA, idxN + j, +cn(j-1))); // add to neutral fluid
				}
			}			
			{
				// mE * Rrec * w1 * U1
				const double c = mE * mI * psi_Rrec(k);
				const Eigen::Vector3d ve = 2./(1.+1./mE) * 1./nE(k) * statesE.col(k).segment(1,3); 
				const Eigen::Vector3d vi = mI / (1 + mE) * 1./nI(k) * statesI.col(k).segment(1,3);
				const Eigen::Vector3d v = ve + vi; // vector in dot-product at timestep m

				const Eigen::Vector3d ce = +c * nE(k) * nI(k) * v; // factors in dot-product with (nu)_e^{m+1}
				const Eigen::Vector3d ci = -c * pow(nE(k),2) * v; // factors in dot-product of (nu)_i^{m+1}
				for (int j = 1; j < 4; j++) {					
					triplets.push_back(T(idxEE, idxE + j, -mE * ce(j - 1))); // subtract from electron fluid
					triplets.push_back(T(idxEE, idxI + j, -mE * ci(j - 1))); // subtract from electron fluid
					triplets.push_back(T(idxEI, idxE + j, +mI * ce(j - 1))); // add to ion fluid
					triplets.push_back(T(idxEI, idxI + j, +mI * ci(j - 1))); // add to ion fluid
				}
			}
		}
	}

	// Create Jacobian matrix from triplets
	const int mSize = 5 * n * nFluids;
	Eigen::SparseMatrix<double> M(mSize, mSize);
	M.setFromTriplets(triplets.begin(), triplets.end());
	M.makeCompressed();
	return M;
}

/**
* This is a reference implementation for inelastic source terms of the Euler equations with explicit discretization. 
*/
Eigen::MatrixXd AppmSolver::getInelasticSourcesExplicit()
{
	const bool showOutput = true;
	const int nFluids = getNFluids();
	const int n = dualMesh.getNumberFluidCells();
	Eigen::MatrixXd src(5 * nFluids, n);
	src.setZero();

	const int nCollisions = inelasticCollisions.size();
	for (int collIdx = 0; collIdx < nCollisions; collIdx++) {
		const InelasticCollision * collision = inelasticCollisions[collIdx];

		// fluid indices
		const int fidxA = collision->getAtomFluidx();
		const int fidxE = collision->getElectronFluidx();
		const int fidxI = collision->getIonFluidx();

		// mass ratios
		const double mA = getSpecies(fidxA).getMassRatio();
		const double mE = getSpecies(fidxE).getMassRatio();
		const double mI = getSpecies(fidxI).getMassRatio();

		// sum of particle masses in collision
		assert(mA == mE + mI);

		/* Ionization energy */
		const double E_ion = collision->getIonizationEnergyScaled();
		//const double E_ion = 0;

		// Fluid states of species (A = neutral atoms, E = electrons, I = ions)
		const Eigen::MatrixXd statesA = getStates(fidxA, n);
		const Eigen::MatrixXd statesE = getStates(fidxE, n);
		const Eigen::MatrixXd statesI = getStates(fidxI, n);

		// Number densities
		const Eigen::VectorXd nA = statesA.row(0);
		const Eigen::VectorXd nE = statesE.row(0);
		const Eigen::VectorXd nI = statesI.row(0);

		// Temperature
		const Eigen::VectorXd TaVec = Physics::getTemperature(statesA, mA);
		const Eigen::VectorXd TeVec = Physics::getTemperature(statesE, mE);
		const Eigen::VectorXd TiVec = Physics::getTemperature(statesI, mI);
		assert(TeVec.allFinite());
		assert((TeVec.array() > 0).all());

		// Velocity vector
		const Eigen::MatrixXd uAmat = Physics::getVelocity(statesA);
		const Eigen::MatrixXd uEmat = Physics::getVelocity(statesE);
		const Eigen::MatrixXd uImat = Physics::getVelocity(statesI);
		assert(uAmat.cols() == n);
		assert(uEmat.cols() == n);
		assert(uImat.cols() == n);

		// Relative velocities in collisions
		const Eigen::Matrix3Xd w0mat = uEmat - uAmat; 
		const Eigen::Matrix3Xd w1mat = mI * (uEmat - uImat);

		// Get ratio of kinetic energy to thermal energy
		Eigen::VectorXd lambdaIon = Eigen::VectorXd::Zero(TeVec.size());
		Eigen::VectorXd lambdaRec = Eigen::VectorXd::Zero(TeVec.size());
		for (int i = 0; i < n; i++) {
			const Eigen::VectorXd w0 = w0mat.col(i);
			const Eigen::VectorXd w1 = w1mat.col(i);
			const double Te = TeVec(i);
			lambdaIon(i) = 0.5 * mE / Te * w0.squaredNorm();
			lambdaRec(i) = 0.5 * mE * Te * w1.squaredNorm();
		}

		// Thermal velocity of electrons
		const Eigen::VectorXd vthE = ((8. / M_PI) * 1. / mE * TeVec).array().sqrt();

		// Ratio of ionization energy to electron thermal energy
		const Eigen::VectorXd xStar = E_ion * TeVec.array().inverse();

		// Coefficients in Le & Cambier (2016), in scale-free numbers
		const Eigen::VectorXd Gion = collision->getGion(nE, nA, vthE, TeVec, lambdaIon);
		const Eigen::VectorXd Grec = collision->getGrec(nI, nE, vthE, xStar, mE, TeVec, lambdaRec);
		const Eigen::VectorXd R0ion = collision->getR0ion(nE, nA, vthE, TeVec, lambdaIon);
		const Eigen::VectorXd R1rec = collision->getR1rec(nE, nI, vthE, xStar, TeVec, lambdaRec);
		const Eigen::VectorXd R2rec = collision->getR2rec(nE, nI, vthE, xStar, TeVec, lambdaRec);
		const Eigen::VectorXd J00ion = collision->getJ00ion(nE, nA, vthE, TeVec, lambdaIon);
		const Eigen::VectorXd J11rec = collision->getJ11rec(nE, nI, vthE, xStar, TeVec, lambdaRec);
		const Eigen::VectorXd J22rec = collision->getJ22rec(nE, nI, vthE, xStar, TeVec, lambdaRec);
		const Eigen::VectorXd J12rec = collision->getJ12rec(nE, nI, vthE, xStar, TeVec, lambdaRec);
		assert(Gion.allFinite());
		assert(Grec.allFinite());
		assert(R0ion.allFinite());
		assert(R1rec.allFinite());
		assert(R2rec.allFinite());
		assert(J00ion.allFinite());
		assert(J11rec.allFinite());
		assert(J22rec.allFinite());
		assert(J12rec.allFinite());

		if (showOutput) {
			int i = 0;
			std::cout << "i = " << i << std::endl;
			std::cout << "nE: " << nE(i) << std::endl;
			std::cout << "nA: " << nA(i) << std::endl;
			std::cout << "vthE: " << vthE(i) << std::endl;
			std::cout << "lambdaIon: " << lambdaIon(i) << std::endl;
			std::cout << "Te: " << TeVec(i) << std::endl;
			std::cout << "Ti: " << TiVec(i) << std::endl;
			std::cout << "Ta: " << TaVec(i) << std::endl;
			std::cout << "Gion: " << Gion(i) << std::endl;
			std::cout << "Grec: " << Grec(i) << std::endl;
		}

		for (int i = 0; i < n; i++) {
			const double Ta = TaVec(i);
			const double Te = TeVec(i);
			const double Ti = TiVec(i);

			const Eigen::VectorXd stateA = statesA.col(i);
			const Eigen::VectorXd stateE = statesE.col(i);
			const Eigen::VectorXd stateI = statesI.col(i);

			// species velocity vector
			const Eigen::Vector3d uA = uAmat.col(i);
			const Eigen::Vector3d uE = uEmat.col(i);
			const Eigen::Vector3d uI = uImat.col(i);

			Eigen::VectorXd srcA = Eigen::VectorXd::Zero(5);
			Eigen::VectorXd srcE = Eigen::VectorXd::Zero(5);
			Eigen::VectorXd srcI = Eigen::VectorXd::Zero(5);


			// Bulk velocity for ionization and recombination
			const Eigen::Vector3d U0 = 1. / (1. + 1. / mE) * uE + 1. / (1. + mE) * uA;
			const Eigen::Vector3d U1 = 2. / (1. + 1. / mE) * uE + mI / (1. + mE) * uI;

			// Relative collision velocities 
			const Eigen::Vector3d w0 = w0mat.col(i);
			const Eigen::Vector3d w1 = w1mat.col(i); 

			// total energy for ionization and recombination in center of mass frame
			const double etot_i_com = 0.5 * (1 + mE) * U0.squaredNorm() + 3. / 2. * Ta;
			const double etot_r_com = 0.5 * (1 + mE) * U1.squaredNorm() + 3. / 2. * Ti; 

			// factors in momentum and energy balance
			const double R_ion = R0ion(i);
			const double K_ion = Gion(i) - R0ion(i);
			const double W_ion = J00ion(i) - 2 * lambdaIon(i) * R0ion(i) + lambdaIon(i) * Gion(i);
			const double J_ion = J00ion(i) - lambdaIon(i) * R0ion(i);
			const double R_rec = R1rec(i) + R2rec(i);
			const double K_rec = 2 * Grec(i) - R1rec(i) - R2rec(i);
			const double W_rec = J11rec(i) + J22rec(i) + 2 * J12rec(i) + 4 * lambdaRec(i) * (Grec(i) - R1rec(i) - R2rec(i));
			const double J_rec = J11rec(i) + J22rec(i) + 2 * J12rec(i) - 2 * lambdaRec(i) * (R1rec(i) + R2rec(i));


			// Number density sources, ionization and recombination
			const double GammaNet = Gion(i) - Grec(i);
			srcA(0) -= GammaNet;
			srcE(0) += GammaNet;
			srcI(0) += GammaNet;

			// Momentum sources
			srcA.segment(1, 3) += -(1 + mE) * (Gion(i) * U0 - Grec(i) * U1) 
				+ mE * (-(Ta - Te) / Te * K_ion * w0 + (Ti - Te)/Te * K_rec * w1 + R_ion * w0);

			srcE.segment(1, 3) += -mE * (R_ion * w0 + R_rec * w1);

			srcI.segment(1, 3) += +(1 + mE) * (Gion(i) * U0 - Grec(i) * U1)
				+ mE * ((Ta - Te) / Te * K_ion * w0 - (Ti - Te) / Te * K_rec * w1 + R_rec * w1);

			// Total energy sources
			srcA(4) += -(Gion(i) * etot_i_com - Grec(i) * etot_r_com) 
				+ 1. / (1. + 1./mE) * (-pow(Ta - Te,2) / Te * W_ion + pow(Ti - Te,2) / Te * W_rec - 2 * (Ta - Te) * J_ion) 
				+ mE * ((Ta - Te)/Te * K_ion * w0.dot(U0) + (Ti - Te)/Te * K_rec * w1.dot(U1) + R_ion * w0.dot(U0));

			srcE(4) += (Grec(i) - Gion(i)) * E_ion
				+ 2. / (1. + 1. / mE) * ((Ta - Te) * J_ion + (Ti - Te) * J_rec)
				- mE * (R_ion * w0.dot(U0) + R_rec * w1.dot(U1));

			srcI(4) += Gion(i) * etot_i_com - Grec(i) * etot_r_com
				+ 1. / (1. + 1./mE) * (pow(Ta - Te, 2) / Te * W_ion - pow(Ti - Te, 2) / Te * W_rec - 2 * (Ti - Te) * J_rec)
				+ mE * (-(Ta - Te) / Te * K_ion * w0.dot(U0) - (Ti - Te) / Te * K_rec * w1.dot(U1) + R_rec * w1.dot(U1));
				
			if (i == 0 && showOutput) {
				std::cout << "Gion, Grec: " << Gion(i) << ", " << Grec(i) << std::endl;
				std::cout << "srcA: " << srcA.transpose() << std::endl;
				std::cout << "srcE: " << srcE.transpose() << std::endl;
				std::cout << "srcI: " << srcI.transpose() << std::endl;

				Eigen::Matrix3Xd data(3, 5);
				data.row(0) = srcA.transpose();
				data.row(1) = srcE.transpose();
				data.row(2) = srcI.transpose();
				//std::ofstream("output.dat") << std::scientific << std::setprecision(20) << data << std::endl;
			}			
			
			// Scaling of momentum and energy sources because they are divided by mass fraction
			srcE.segment(1, 4) *= 1. / mE;
			srcI.segment(1, 4) *= 1. / mI;

			// Update data structure with local species sources
			Eigen::VectorXd localSrc(5 * nFluids);
			localSrc.setZero();
			localSrc.segment(5 * fidxA, 5) += srcA;
			localSrc.segment(5 * fidxE, 5) += srcE;
			localSrc.segment(5 * fidxI, 5) += srcI;
			src.col(i) += localSrc;
		}
	}
	return src;
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
	//M += dt * Q.transpose() * Msigma_inner * Q; // Test: leave this term out. To be activated again
	M.makeCompressed();
	//Eigen::sparseMatrixToFile(M, "M.dat");

	double dt_ratio = dt / dt_previous;

	// Data vector for right hand side
	Eigen::VectorXd rhs(x.size());
	rhs.setZero();
	rhs = M1 * (1 + dt_ratio) * maxwellState;
	rhs -= dt_ratio * M1 * maxwellStatePrevious;
	rhs -= dt * Q.transpose() * (J_h_aux - J_h).topRows(nEdges); 
	//rhs -= dt * Q.transpose() * (J_h_aux - J_h).topRows(nEdges); // Test: leave this term out


	// The system of equations has fixed and free values; 
	// - fixed values: electric potential at terminals (Dirichlet boundary condition)
	// - free  values: electric potential at non-terminal vertices, and electric voltages at non-boundary edges
	int nDirichlet = nPrimalTerminalVertices;
	if (solverParams.getMaxwellCurrentDefined()) {
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
	switch (solverParams.getMaxwellSolverType()) {
	case MaxwellSolverType::CG:
		// remark: Matrix Mf is numerically not SPD. Is it symmetric non-negative definite? 
		xf = solveMaxwell_CG(Mf, rhsFree); 
		break;

	case MaxwellSolverType::PardisoLU:
		xf = solveMaxwell_PardisoLU(Mf, rhsFree);
		break;

	case MaxwellSolverType::BiCGStab:
		xf = solveMaxwell_BiCGStab(Mf, rhsFree);
		break;


	default:
		exit(-1);
	}
	
	// Test: apply guard against truncation errors by adding and subtracting a large number, 
	//       so that small values (of order 1e-10 smaller than largest number) goes to zero 
	//       due to finite precision arithmetic
	//const double x_maxCoeff = xf.cwiseAbs().maxCoeff();
	//const double scale = 1e6;
	//xf.array() += x_maxCoeff * scale;
	//xf.array() -= x_maxCoeff * scale;

	// assemble new state vector
	x.topRows(nFree) = xf;
	x.bottomRows(nDirichlet) = xd;
	//std::ofstream("x.dat") << x << std::endl;

	// Update state vectors for discrete states
	E_h = Q * x;

	const double Emax = E_h.cwiseAbs().maxCoeff();
	const double scale = 1e9;
	E_h.array() += scale * Emax;
	E_h.array() -= scale * Emax;

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
	std::cout << "Primal mesh volume: " << primalMesh.getMeshVolume() << std::endl;
	std::cout << "Dual mesh volume:   " << dualMesh.getMeshVolume() << std::endl;

	std::cout << "Primal mesh" << std::endl;
	std::cout << primalMesh.getMeshInfo() << std::endl;

	std::cout << "Dual mesh" << std::endl;
	std::cout << dualMesh.getMeshInfo() << std::endl;
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
	std::cout << "Write XDMF output file with " << nTimesteps << " timesteps" << std::endl;
	
	std::string gridPrimalEdges;
	std::string gridPrimalFaces;
	std::string gridDualEdges;
	std::string gridDualFaces;

	assert(timeStamps.size() == outputIterations.size());

	std::ofstream file(filename);
	file << "<?xml version = \"1.0\" ?>" << std::endl;
	file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
	file << "<Xdmf Version=\"3.0\" xmlns:xi=\"[http://www.w3.org/2001/XInclude]\">" << std::endl;
	file << "<Domain>" << std::endl;
	file << "<Grid Name=\"Time Grid\" GridType=\"Collection\" CollectionType=\"Temporal\">" << std::endl;
	for (int i = 0; i < nTimesteps; i++) {
		const double time = this->timeStamps[i];
		const int iteration = outputIterations[i];

		file << "<Grid Name=\"Grid of Grids\" GridType=\"Tree\">" << std::endl;
		file << "<Time Value=\"" << time << "\" />" << std::endl;
		file << xdmf_GridPrimalEdges(iteration) << std::endl;
		file << xdmf_GridPrimalFaces(iteration) << std::endl;
		file << xdmf_GridDualEdges(iteration) << std::endl;
		file << xdmf_GridDualFaces(iteration) << std::endl;
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
		const int iteration = this->outputIterations[i];
		file << "<Time Value=\"" << time << "\" />" << std::endl;
		file << xdmf_GridDualCells(iteration) << std::endl;
	}
	file << "</Grid>" << std::endl;
	file << "</Domain>" << std::endl;
	file << "</Xdmf>" << std::endl;
}

/**
* Get formatted string of output filename.
*
* @param iteration  iteration number 
* @return formatted string of output filename
*/
std::string getOutputFilename(const int iteration) {
	return (std::stringstream() << "appm-" << std::setfill('0') << std::setw(5) << iteration << ".h5").str();
}

/**
* Write data to output file.
*/
void AppmSolver::writeOutput(const int iteration, const double time)
{
	std::cout << "Write output at iteration " << iteration << ", time = " << time << std::endl;
	outputIterations.push_back(iteration);
	timeStamps.push_back(time);


	const std::string filename = getOutputFilename(iteration); 

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
		const std::string fluidName = getSpecies(fluidIdx).getName();
		const std::string fluidTag = (std::stringstream() << "/" << fluidName << "-").str();
		const std::string pressureTag = fluidTag + "pressure";
		const std::string velocityTag = fluidTag + "velocity";
		const std::string densityTag = fluidTag + "density";
		const std::string numberDensityTag = fluidTag + "numberDensity";
		Eigen::VectorXd numberDensity = fluidStates.row(5 * fluidIdx);
		Eigen::VectorXd density(nCells);
		Eigen::MatrixXd velocity(3, nCells);
		Eigen::VectorXd pressure(nCells);
		Eigen::VectorXd temperature(nCells);
		Eigen::VectorXd sumMassFluxes = sumOfFaceFluxes.row(5 * fluidIdx);
		Eigen::MatrixXd sumMomentumFluxes = sumOfFaceFluxes.block(5 * fluidIdx + 1, 0, 3, nCells);
		Eigen::VectorXd sumEnergyFluxes = sumOfFaceFluxes.row(5*fluidIdx + 4);
		Eigen::Matrix3Xd frictionSource = frictionForceSourceTerm.block(3 * fluidIdx, 0, 3, nCells);
		Eigen::VectorXd speciesFrictionEnergySource = frictionEnergySourceTerm.row(fluidIdx);
		Eigen::Matrix3Xd speciesDiffusionVelocity = diffusionVelocity.block(3 * fluidIdx, 0, 3, nCells);
		Eigen::VectorXd fluidMassFluxImplicitTerm = massFluxImplicitTerm.row(fluidIdx);
		Eigen::Matrix3Xd speciesMomentumSource = fluidSources.block(5 * fluidIdx + 1, 0, 3, nCells);
		Eigen::VectorXd speciesEnergySource = fluidSources.row(5 * fluidIdx + 4);

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
		temperature.setConstant(std::nan(""));
		qN.setConstant(std::nan("")); 
		qU.setConstant(std::nan("")); 
		qE.setConstant(std::nan("")); 

		const double epsilon2 = getSpecies(fluidIdx).getMassRatio();
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
			//try {
			Physics::state2primitive(epsilon2, state, n, p, u);
			//}
			//catch (std::exception & e) {
			//	std::cout << "Exception: " << e.what() << std::endl;
			//	std::cout << "cell idx: " << i << std::endl;
			//	std::cout << "state:    " << state.transpose() << std::endl;

			//	Eigen::VectorXd sumFaceFluxes = sumOfFaceFluxes.col(i).segment(5 * fluidIdx, 5);
			//	std::cout << "sumOfFaceFluxes: " << sumFaceFluxes.transpose() << std::endl;
			//	createStopFile(1); // write a positive value into stop file, this indicates that the iteration loop should terminate.
			//	//assert(false);
			//}
			if (p < 0 || n < 0) {
				std::cout << "Invalid state " << std::endl;
				std::cout << "fluid: " << getSpecies(fluidIdx).getName() << std::endl;
				std::cout << "cell idx: " << i << std::endl;
				createStopFile(1);
				break;
			}
			density(i) = epsilon2 * n;
			velocity.col(i) = u;
			pressure(i) = p;
			const double T = p / n;
			temperature(i) = T;
		}

		writer.writeData(density, densityTag);
		writer.writeData(numberDensity, numberDensityTag);
		writer.writeData(pressure, pressureTag);
		writer.writeData(velocity, velocityTag);
		writer.writeData(temperature, (std::stringstream() << fluidTag << "temperature").str());
		writer.writeData(sumMassFluxes, (std::stringstream() << fluidTag << "sumMassFlux").str());
		writer.writeData(sumMomentumFluxes, (std::stringstream() << fluidTag << "sumMomentumFlux").str());
		writer.writeData(sumEnergyFluxes, (std::stringstream() << fluidTag << "sumEnergyFlux").str());
		writer.writeData(frictionSource, (std::stringstream() << fluidTag << "frictionSource").str());
		writer.writeData(speciesDiffusionVelocity, (std::stringstream() << fluidTag << "diffusionVelocity").str());
		writer.writeData(speciesFrictionEnergySource, (std::stringstream() << fluidTag << "frictionEnergySource").str());
		writer.writeData(fluidMassFluxImplicitTerm, (std::stringstream() << fluidTag << "massFluxImplicitTerm").str());
		writer.writeData(speciesEnergySource, (std::stringstream() << fluidTag << "energySource").str());
		writer.writeData(speciesMomentumSource, (std::stringstream() << fluidTag << "momentumSource").str());

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
	writer.writeData(bulkVelocity, "/bulkVelocity");
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
		const Eigen::Vector3d fn = face->getNormal().normalized();
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
		const Eigen::Vector3d fn = face->getNormal().normalized();
		currentDensity.col(i) = J_h(i) / fA * fn;
	}
	writer.writeData(currentDensity, "/CurrentDensity");


	assert(dualMesh.getNumberOfFaces() == J_h_aux.size());
	assert(J_h_aux.allFinite());
	Eigen::Matrix3Xd J_h_aux_vector(3, J_h_aux.size());
	for (int i = 0; i < J_h_aux.size(); i++) {
		const Face * face = dualMesh.getFace(i);
		const double fA = face->getArea();
		const Eigen::Vector3d fn = face->getNormal().normalized();
		J_h_aux_vector.col(i) = J_h_aux(i) / fA * fn;
	}
	writer.writeData(J_h_aux_vector, "/j_h_aux_vector");
	writer.writeData(J_h_aux, "/J_h_aux");
	writer.writeData(J_h, "/J_h");

	Jcc = getCurrentDensityAtCellCenter();
	writer.writeData(Jcc, "/Jcc");
	writer.writeData(Jaux_cc, "/Jaux_cc");

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

//void AppmSolver::readParameters(const std::string & filename)
//{
//	appmParams = AppmParameters();
//
//	std::ifstream file(filename);
//	if (!file.is_open()) {
//		std::cout << "File not opened: " << filename;
//		exit(-1);
//	}
//
//	std::string line;
//	const char delim = ':';
//
//	while (std::getline(file, line)) {
//		if (line.size() == 0) { continue; } // skip empty lines
//		int pos = line.find(delim);
//		std::string tag = line.substr(0, pos);
//
//		if (tag == "maxIterations") {
//			std::istringstream(line.substr(pos + 1)) >> appmParams.maxIterations;
//		}
//		if (tag == "maxTime") {
//			std::istringstream(line.substr(pos + 1)) >> appmParams.maxTime;
//		}
//		if (tag == "timestepSize") {
//			std::istringstream(line.substr(pos + 1)) >> appmParams.timestepSize;
//		}
//		if (tag == "isFluidEnabled") {
//			std::istringstream(line.substr(pos + 1)) >> appmParams.isFluidEnabled;
//		}
//		if (tag == "isMaxwellEnabled") {
//			std::istringstream(line.substr(pos + 1)) >> appmParams.isMaxwellEnabled;
//		}
//		if (tag == "lambdaSquare") {
//			std::istringstream(line.substr(pos + 1)) >> this->lambdaSquare;
//		}
//		if (tag == "isMassFluxSchemeImplicit") {
//			std::istringstream(line.substr(pos + 1)) >> appmParams.isMassFluxSchemeImplicit;
//		}
//		if (tag == "initType") {
//			std::istringstream(line.substr(pos + 1)) >> initType;
//		}
//		if (tag == "isElectricLorentzForceActive") {
//			std::istringstream(line.substr(pos + 1)) >> appmParams.isLorentzForceElectricEnabled;
//		}
//		if (tag == "isMagneticLorentzForceActive") {
//			std::istringstream(line.substr(pos + 1)) >> appmParams.isLorentzForceMagneticEnabled;
//		}
//		if (tag == "isShowDataWriterOutput") {
//			std::istringstream(line.substr(pos + 1)) >> isShowDataWriterOutput;
//		}
//		if (tag == "isEulerMaxwellCouplingEnabled") {
//			std::istringstream(line.substr(pos + 1)) >> appmParams.isEulerMaxwellCouplingEnabled;
//		}
//		if (tag == "isMaxwellCurrentDefined") {
//			std::istringstream(line.substr(pos + 1)) >> appmParams.isMaxwellCurrentDefined;
//		}
//		if (tag == "isFrictionActive") {
//			std::istringstream(line.substr(pos + 1)) >> appmParams.isFrictionActive;
//		}
//		if (tag == "MaxwellSolverType") {
//			bool isValid = false;
//			std::string value;
//			std::istringstream(line.substr(pos + 1)) >> value;
//			std::cout << "Maxwell solver type: value = " << value << std::endl;
//			if (value.compare("CG") == 0) {
//				appmParams.maxwellSolverType = MaxwellSolverType::CG;
//				isValid = true;
//			}
//			if (value.compare("PardisoLU") == 0) {
//				appmParams.maxwellSolverType = MaxwellSolverType::PardisoLU;
//				isValid = true;
//			}
//			if (value.compare("BiCGStab") == 0) {
//				appmParams.maxwellSolverType = MaxwellSolverType::BiCGStab;
//				isValid = true;
//			}
//			assert(isValid);
//		}
//	}
//
//	std::cout << std::endl;
//	std::cout << printSolverParameters() << std::endl;
//}

const std::string AppmSolver::xdmf_GridPrimalEdges(const int iteration) const
{
	const std::string meshFilename = primalMesh.getMeshDataFilename();
	std::stringstream ss;
	ss << "<Grid Name=\"Primal Edges\">" << std::endl;
	ss << "<Topology TopologyType=\"Polyline\""
		<< " NumberOfElements=\"" << primalMesh.getNumberOfEdges() << "\""
		<< " NodesPerElement=\"2\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << 2 * primalMesh.getNumberOfEdges() << "\" DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << meshFilename << ":/edge2vertex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Topology>" << std::endl;

	ss << "<Geometry GeometryType=\"XYZ\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << primalMesh.getNumberOfVertices() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << meshFilename << ":/vertexPos" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Geometry>" << std::endl;

	ss << "<Attribute Name=\"Edge index\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << primalMesh.getNumberOfEdges() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << meshFilename << ":/edgeIdx" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	if (isWriteEfield) {
		const std::string dataFilename = getOutputFilename(iteration);
		ss << "<Attribute Name=\"Electric field\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << primalMesh.getNumberOfEdges() << " 3\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << dataFilename << ":/E" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;
	}

	ss << "</Grid>";
	return ss.str();
}
const std::string AppmSolver::xdmf_GridPrimalFaces(const int iteration) const
{
	const std::string meshFilename = primalMesh.getMeshDataFilename();
	H5Reader h5reader;
	h5reader = H5Reader(meshFilename);
	const int nElements = h5reader.readDataSize("/face2vertex");
	assert(nElements > 0);

	std::stringstream ss;
	ss << "<Grid Name=\"Primal Faces\">" << std::endl;
	ss << "<Topology TopologyType=\"Mixed\""
		<< " NumberOfElements=\"" << primalMesh.getNumberOfFaces() << "\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << nElements << "\" DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << meshFilename << ":/face2vertex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Topology>" << std::endl;

	ss << "<Geometry GeometryType=\"XYZ\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << primalMesh.getNumberOfVertices() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << meshFilename << ":/vertexPos" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Geometry>" << std::endl;

	ss << "<Attribute Name=\"Face index\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << primalMesh.getNumberOfFaces() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << meshFilename << ":/faceIndex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"Vertex type\" AttributeType=\"Scalar\" Center=\"Node\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << primalMesh.getNumberOfVertices() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << meshFilename << ":/vertexType" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	const std::string dataFilename = getOutputFilename(iteration);
	if (isWriteBfield) {
		ss << "<Attribute Name=\"Magnetic Flux\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << primalMesh.getNumberOfFaces() << " 3\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << dataFilename << ":/B" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;
	}

	ss << "<Attribute Name=\"Magnetic Flux Vertex Centered\" AttributeType=\"Vector\" Center=\"Node\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << primalMesh.getNumberOfVertices() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << dataFilename << ":/Bcc" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"Magnetic Flux Vertex Centered Mag\" AttributeType=\"Scalar\" Center=\"Node\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfCells() << "\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << dataFilename << ":/Bcc_mag" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "</Grid>";


	return ss.str();
}

const std::string AppmSolver::xdmf_GridDualEdges(const int iteration) const
{
	const std::string meshFilename = dualMesh.getMeshDataFilename();
	std::stringstream ss;
	ss << "<Grid Name=\"Dual Edges\">" << std::endl;
	ss << "<Topology TopologyType=\"Polyline\""
		<< " NumberOfElements=\"" << dualMesh.getNumberOfEdges() << "\"" 
		<< " NodesPerElement=\"2\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << 2 * dualMesh.getNumberOfEdges() << "\" DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << meshFilename << ":/edge2vertex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Topology>" << std::endl;

	ss << "<Geometry GeometryType=\"XYZ\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfVertices() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << meshFilename << ":/vertexPos" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Geometry>" << std::endl;

	ss << "<Attribute Name=\"Edge index\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfEdges() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << meshFilename << ":/edgeIdx" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	const std::string dataFilename = getOutputFilename(iteration);
	if (isWriteHfield) {
		ss << "<Attribute Name=\"Magnetic field\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfEdges() << " 3\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << dataFilename << ":/H" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;
	}
	ss << "</Grid>" << std::endl;
	return ss.str();
}

const std::string AppmSolver::xdmf_GridDualFaces(const int iteration) const
{
	const std::string meshFilename = dualMesh.getMeshDataFilename();
	H5Reader h5reader;
	h5reader = H5Reader(meshFilename);
	const int nElements = h5reader.readDataSize("/face2vertex");
	assert(nElements > 0);

	std::stringstream ss;
	ss << "<Grid Name=\"Dual Faces\">" << std::endl;
	ss << "<Topology TopologyType=\"Mixed\""
		<< " NumberOfElements=\"" << dualMesh.getNumberOfFaces() << "\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << nElements << "\" DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << meshFilename << ":/face2vertex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Topology>" << std::endl;

	ss << "<Geometry GeometryType=\"XYZ\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfVertices() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << meshFilename << ":/vertexPos" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Geometry>" << std::endl;

	ss << "<Attribute Name=\"Face index\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << meshFilename << ":/faceIndex" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;
	
	ss << "<Attribute Name=\"Face Fluid Type\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << "\""
		<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
	ss << meshFilename << ":/faceFluidType" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"" << "Face Normal" << "\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << meshFilename << ":/faceNormal" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	const std::string dataFilename = getOutputFilename(iteration);

	for (int fidx = 0; fidx < getNFluids(); fidx++) {
		const std::string speciesName = getSpecies(fidx).getName();
		std::string attributeName = (std::stringstream() << speciesName <<" Face Type").str();
		std::string dataName = (std::stringstream() << speciesName <<"-faceType").str();
		ss << "<Attribute Name=\"" << attributeName << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << "\""
			<< " DataType=\"Int\" Precision=\"4\" Format=\"HDF\">" << std::endl;
		ss << dataFilename << ":/" << dataName << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;

		attributeName = (std::stringstream() << speciesName << " Face Mass Flux").str();
		ss << "<Attribute Name=\"" << attributeName << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << "\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << dataFilename << ":/" << speciesName << "-massFlux" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;

		attributeName = (std::stringstream() << speciesName << " Face Momentum Flux").str();
		ss << "<Attribute Name=\"" << attributeName << "\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << " 3\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << dataFilename << ":/" << speciesName << "-momentumFlux" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;

		attributeName = (std::stringstream() << speciesName << " Face Energy Flux").str();
		ss << "<Attribute Name=\"" << attributeName << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << "\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << dataFilename << ":/" << speciesName << "-energyFlux" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;

		attributeName = (std::stringstream() << speciesName << " Implicit Mass Flux Term").str();
		ss << "<Attribute Name=\"" << attributeName << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << "\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << dataFilename << ":/" << speciesName << "-massFluxImplicitTerm" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;
	}

	ss << "<Attribute Name=\"Current density\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << dataFilename << ":/CurrentDensity" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;
	

	ss << "<Attribute Name=\"Current density aux\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << dataFilename << ":/j_h_aux_vector" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"J h aux\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << "\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << dataFilename << ":/J_h_aux" << std::endl;
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

	const std::string datafilename = getOutputFilename(iteration); 

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

	ss << "<Attribute Name=\"Current Density Aux Cell Centered\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfCells() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << datafilename << ":/Jaux_cc" << std::endl;
	ss << "</DataItem>" << std::endl;
	ss << "</Attribute>" << std::endl;

	ss << "<Attribute Name=\"" << "Bulk Velocity" << "\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
	ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfCells() << " 3\""
		<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
	ss << datafilename << ":/" << "bulkVelocity" << std::endl;
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
		const std::string speciesName = getSpecies(k).getName();
		const std::string fluidName = (std::stringstream() << speciesName).str();

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

		ss << "<Attribute Name=\"" << fluidName << " temperature" << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << nCells << "\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << datafilename << ":/" << fluidName << "-temperature" << std::endl;
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

		ss << "<Attribute Name=\"" << fluidName << " Diffusion Velocity" << "\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfCells() << " 3\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << datafilename << ":/" << fluidName << "-diffusionVelocity" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;

		ss << "<Attribute Name=\"" << fluidName << " Friction Force Source Term" << "\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfCells() << " 3\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << datafilename << ":/" << fluidName << "-frictionSource" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;

		ss << "<Attribute Name=\"" << fluidName << " Friction Energy Source Term" << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfCells() << "\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << datafilename << ":/" << fluidName << "-frictionEnergySource" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;

		ss << "<Attribute Name=\"" << fluidName << " Momentum Source" << "\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfCells() << " 3\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << datafilename << ":/" << fluidName << "-momentumSource" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;

		ss << "<Attribute Name=\"" << fluidName << " Energy Source Term" << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
		ss << "<DataItem Dimensions=\"" << dualMesh.getNumberOfCells() << "\""
			<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
		ss << datafilename << ":/" << fluidName << "-energySource" << std::endl;
		ss << "</DataItem>" << std::endl;
		ss << "</Attribute>" << std::endl;

		if (isStateWrittenToOutput) {
			ss << "<Attribute Name=\"" << fluidName << " state N" << "\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
			ss << "<DataItem Dimensions=\"" << nCells << "\""
				<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
			ss << datafilename << ":/" << fluidName << "-stateN" << std::endl;
			ss << "</DataItem>" << std::endl;
			ss << "</Attribute>" << std::endl;


			ss << "<Attribute Name=\"" << fluidName << " state U\" AttributeType=\"Vector\" Center=\"Cell\">" << std::endl;
			ss << "<DataItem Dimensions=\"" << nCells << " 3\""
				<< " DataType=\"Float\" Precision=\"8\" Format=\"HDF\">" << std::endl;
			ss << datafilename << ":/" << fluidName << "-stateU" << std::endl;
			ss << "</DataItem>" << std::endl;
			ss << "</Attribute>" << std::endl;

			ss << "<Attribute Name=\"" << fluidName << " state E\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
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

const double AppmSolver::getNumericalSchemeFactor(const Face * face, const int fluidx) const
{
	assert(face != nullptr);
	double numSchemeFactor = 1.0;
	const Face::Type faceType = getFaceTypeOfFluid(face, fluidx);
	if (faceType == Face::Type::INTERIOR) {
		numSchemeFactor = 0.5;
	}
	return numSchemeFactor;
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
	//Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> maxwellSolver;
	//maxwellSolver.preconditioner().setDroptol(0.1);
	//maxwellSolver.preconditioner().setFillfactor(0.1);
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> maxwellSolver;
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
	//xf = maxwellSolver.solve(rhs);
	Eigen::VectorXd guess = maxwellState.topRows(rhs.size());
	xf = maxwellSolver.solveWithGuess(rhs, guess);
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
	//Eigen::SparseMatrix<double> Mf_t = Mf.transpose();
	//Eigen::SparseMatrix<double> A = 0.5 * (Mf + Mf_t);
	maxwellSolver.compute(Mf);
	auto timer_afterCompute = std::chrono::high_resolution_clock::now();

	auto delta_afterCompute = std::chrono::duration<double>(timer_afterCompute - timer_start);
	std::cout << "Solver time for compute step: " << delta_afterCompute.count() << std::endl;

	if (maxwellSolver.info() != Eigen::Success) {
		std::cout << "Maxwell solver setup failed" << std::endl;
	}

	auto timer_startSolve = std::chrono::high_resolution_clock::now();
	//xf = maxwellSolver.solve(rhs);
	Eigen::VectorXd guess = maxwellState.topRows(rhs.size());
	xf = maxwellSolver.solveWithGuess(rhs, guess);
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

	if (solverParams.getMaxwellCurrentDefined()) {
		// Define a uniform current density for those faces being normal to z-axis
		const double t0 = 2;
		const double tscale = 0.5;
		const double currentDensity = 0.5 * (1 + tanh((time + dt - t0) / tscale)); // current density at implicit timestep t + dt
		const double radius = 0.5;
		const double tol = 2 * std::numeric_limits<double>::epsilon();

		for (int i = 0; i < nEdges; i++) {
			const Face * dualFace = dualMesh.getFace(i);
			const Eigen::Vector3d fn = dualFace->getNormal().normalized();
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
			const Eigen::Vector3d nhat = face->getNormal().normalized();

			orientation(i) = Lhat.dot(nhat); // this value should be equal to 1.000
		}
		orientation.array() -= 1;
		const double tol = 4 * std::numeric_limits<double>::epsilon();
		assert((orientation.cwiseAbs().array() <= tol).all());
	}

	const int nFluids = getNFluids();
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;
	std::vector<T> geomTriplets;

	if (solverParams.getFluidEnabled()) {
		std::cout << "isEulerMaxwellCouplingEnabled: " << solverParams.getEulerMaxwellCouplingEnabled() << std::endl;
		if (solverParams.getEulerMaxwellCouplingEnabled()) {

			// accumulator to guard against truncation errors
			Eigen::VectorXd c = Eigen::VectorXd::Zero(nDualFaces);

			assert(nEdges <= faceFluxes.cols());
			// Loop over dual faces (where current density lives). 
			// Do not consider dual boundary faces (those coincident with domain boundary), 
			// this leaves us with inner faces of the dual mesh 
			for (int i = 0; i < nEdges; i++) {
				//std::cout << "i = " << i << std::endl;
				const Face * dualFace = dualMesh.getFace(i);
				const Eigen::Vector3d ni = dualFace->getNormal().normalized();
				const double Ai = dualFace->getArea();
				const Face::Type faceType = dualFace->getType();
				const double numSchemeFactor = (faceType == Face::Type::INTERIOR) ? 0.5 : 1;

				for (int fluidIdx = 0; fluidIdx < nFluids; fluidIdx++) {
					const int q = getSpecies(fluidIdx).getCharge();
					const double massRatio = getSpecies(fluidIdx).getMassRatio();
					const Face::Type faceFluidType = getFaceTypeOfFluid(dualFace, fluidIdx);

					// Term due to explicit mass flux. Note that face fluxes are already multiplied with face area Ai
					const double massFlux = faceFluxes(5 * fluidIdx + 0, i);
					Jaux(i) += q * massFlux;

					// Add extra terms due to implicit formulation
					if (solverParams.getMassfluxSchemeImplicit()) {
						auto adjacientCells = dualFace->getCellList();
						for (auto cell : adjacientCells) {
							// non-fluid cells
							if (cell->getType() != Cell::Type::FLUID) {
								//if (dualFace->getType() == Face::Type::DEFAULT) {
								//	const double elCondSolid = 1e-3;
								//	const double Li = primalMesh.getEdge(i)->getLength();
								//	const double value = elCondSolid * Ai / Li;
								//	triplets.push_back(T(i, i, elCondSolid));
								//}
							}
							else { // is fluid cell
								assert(cell->getType() == Cell::Type::FLUID);
								const double nk = fluidStates(5 * fluidIdx + 0, cell->getIndex()); // number density in cell k
								const double Vk = cell->getVolume(); // volume of cell k

								// Extra term by explicit Fluid sources (e.g., magnetic Lorentz force, or friction force)
								if (solverParams.getLorentzForceEnabled() || solverParams.getFrictionActive()) {
									const Eigen::Vector3d fluidMomentumSource = fluidSources.col(cell->getIndex()).segment(5 * fluidIdx + 1, 3);
									Jaux(i) += numSchemeFactor * q * Ai * dt * fluidMomentumSource.dot(ni);
								}

								auto cellFaces = cell->getFaceList();
								for (auto face : cellFaces) {
									const double Aj = face->getArea();
									const Eigen::Vector3d nj = face->getNormal().normalized();
									const int s_kj = cell->getOrientation(face);
									const int j = face->getIndex();
									double nj_dot_ni = nj.dot(ni); // angle between face i and face j
									if (abs(nj_dot_ni) < 1e-10) {
										nj_dot_ni = 0; // Truncate small angles
									}

									// Extra term by explicit momentum flux
									if (true) {
										const Eigen::Vector3d fj = faceFluxes.col(j).segment(5 * fluidIdx + 1, 3); // (directional) momentum flux at face j, times face area
										const double temp = numSchemeFactor * q * Ai * -dt * 1 / Vk * s_kj * fj.dot(ni);
										c(i) += abs(temp); // accumulator to guard against truncation error
										Jaux(i) += temp;
									}

									// Extra term by implicit electric Lorentz force (electric field)
									const double geomFactor = nj_dot_ni * (Ai * Aj) / Vk; // geometric factor
									if (solverParams.getLorentzForceEnabled()) {
										const double value = numSchemeFactor * 0.5 * dt * pow(q, 2) / massRatio * nk / Vk * Ai * Aj * nj_dot_ni;
										triplets.push_back(T(i, j, value));
									}
									if (fluidIdx == 0) {
										geomTriplets.push_back(T(i, j, geomFactor));
									}
								}
							} // end if cellType != Fluid
						}
					} // end if isMassFluxSchemeImplicit
				}
			}

			// guard against truncation errors
			//std::ofstream("c.dat")    << std::scientific << std::setprecision(16) << c << std::endl;
			//std::ofstream("jaux.dat") << std::scientific << std::setprecision(16) << Jaux << std::endl;
			if (true) {
				const double scale = 1000 * this->fluxTruncationErrorGuardScale;
				Jaux += c * scale;
				Jaux -= c * scale;
			}

			//Eigen::SparseMatrix<double> Msigma_geom(nDualFaces, nEdges);
			//Msigma_geom.setFromTriplets(geomTriplets.begin(), geomTriplets.end());
			//Msigma_geom.makeCompressed();
			//std::cout << "Msigma_geom triplets size: " << geomTriplets.size() << std::endl;
			//std::cout << "Msigma_geom nnz: " << Msigma_geom.nonZeros() << std::endl;
			//Eigen::sparseMatrixToFile(Msigma_geom, "Msigma_geom.dat");
		} // end if isEulerMaxwellCouplingEnabled
	}
	else
	{
		for (int i = 0; i < nEdges; i++) {
			const double elCond = 1;
			const double dL = primalMesh.getEdge(i)->getLength();
			const double dA = dualMesh.getFace(i)->getArea();
			const double value = elCond * dA / dL;
			triplets.push_back(T(i, i, value));
		}
	}
	Msigma.setFromTriplets(triplets.begin(), triplets.end());
	Msigma.makeCompressed();

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
	if (solverParams.getMaxwellCurrentDefined()) {
		nDirichlet = nPrimalTerminalVertices;
	}

	Eigen::VectorXd xd = Eigen::VectorXd::Zero(nDirichlet);

	if (solverParams.getMaxwellCurrentDefined()) {
		xd.topRows(nPrimalTerminalVertices / 2).array() = 0;
	}
	else {
		// Set voltage values at terminal A (stressed electrode)
		auto finalVoltage = 0.25;
		auto t0 = 0.3;
		auto t1 = 0.5;
		auto tscale = 0.05;
		auto timeFunction0 = 0.5 * (1 + tanh((time - t0) / tscale));
		auto timeFunction1 = 0.5 * (1 + tanh(-(time - t1) / tscale));
		auto voltage = finalVoltage * timeFunction0 * timeFunction1;
		voltage = 1e0 * (time > 0);
		xd.topRows(nPrimalTerminalVertices / 2).array() = voltage;
	}
	// Set voltage condition to zero at terminal B (grounded electrode)
	xd.bottomRows(nPrimalTerminalVertices / 2).array() = 0;

	return xd;
}

const Species & AppmSolver::getSpecies(const int idx) const
{
	assert(idx >= 0);
	assert(idx < speciesList.size());
	return speciesList[idx];
}

const int AppmSolver::getSpeciesIndex(const std::string & tag)
{
	int result = -1;
	for (int i = 0; i < speciesList.size(); i++) {
		Species & species = speciesList[i];
		if (species.getSymbol() == tag) {
			result = i;
		}
	}
	if (result < 0) {
		std::cout << std::endl << std::endl;
		std::cout << "Species index not found; tag = " << tag << std::endl;
		std::cout << "Defined species: " << std::endl;
		for (auto species : speciesList) {
			std::cout << species.getSymbol() << std::endl;
		}
	}
	assert(result >= 0);
	return result;
}

const double AppmSolver::getCollisionFrequency(const int alpha, const int beta, const int cellIdx)
{
	auto nFluids = getNFluids();
	assert(alpha >= 0 && alpha < nFluids);
	assert(beta >= 0 && beta < nFluids);
	assert(cellIdx >= 0);

	const Eigen::VectorXd stateA = fluidStates.col(cellIdx).segment(5 * alpha, 5);
	const Eigen::VectorXd stateB = fluidStates.col(cellIdx).segment(5 * beta, 5);
	const double n_b = stateB(0);

	auto m_a = getSpecies(alpha).getMassRatio();
	auto m_b = getSpecies(beta).getMassRatio();
	auto m_ab = m_a * m_b / (m_a + m_b); // reduced mass
	auto T_a = Physics::getTemperature(stateA, m_a); // local temperature of species A 
	auto T_b = Physics::getTemperature(stateB, m_b); // local temperature of species B
	auto T_ab = (m_a * T_b + m_b * T_a) / (m_a + m_b); // reduced temperature
	auto u_therm_ab = std::sqrt(8 / M_PI * T_ab / m_ab); // thermal velocity associated to reduced temperature
	
	double nu_ab = 0;
	// TODO
	//if () {
	//	nu_ab = n_b;
	//}
	//else {
	//	auto Q_ab = 1; // average momentum transfer collision cross-section
	//	nu_ab = 4. / 3. * n_b * u_therm_ab * Q_ab; // collision frequency
	//}

	return nu_ab;
}

/** 
* @return Get reduced mass of the interaction between fluid A and B.
*/
const double AppmSolver::getReducedMass(const int alpha, const int beta) const
{
	auto nFluids = getNFluids();
	assert(alpha >= 0 && alpha < nFluids);
	assert(beta >= 0 && beta < nFluids);
	auto m_a = getSpecies(alpha).getMassRatio();
	auto m_b = getSpecies(beta).getMassRatio();
	auto m_ab = m_a * m_b / (m_a + m_b);
	return m_ab;
}

/**
* Get fluid states.
*
* @param fidx    index of fluid
* @param nCols    number of columns (i.e., cells with index 0 to nCols-1)
* @return fluid states of fluid with index fidx.
*/
const Eigen::MatrixXd AppmSolver::getStates(const int fidx, const int nCols) const
{
	const int nFluids = getNFluids();
	assert(fidx >= 0);
	assert(fidx < nFluids);
	return fluidStates.block(5*fidx, 0, 5, nCols);
}

const int AppmSolver::getLinearIndexInJacobian(const int fluidIdx, const int cellidx) const
{
	const int nFluids = getNFluids();
	const int nFluidCells = dualMesh.getNumberFluidCells();
	assert(fluidIdx >= 0);
	assert(fluidIdx < nFluids);
	assert(cellidx >= 0);
	assert(cellidx < nFluidCells);

	int idx = 5 * (nFluids * cellidx + fluidIdx);
	return idx;
}



std::ostream & operator<<(std::ostream & os, const AppmSolver::MaxwellSolverType & obj)
{
	switch (obj) {
		case AppmSolver::MaxwellSolverType::CG:
			os << "CG";
			break;

		case AppmSolver::MaxwellSolverType::PardisoLU:
			os << "PardisoLU";
			break;

		case AppmSolver::MaxwellSolverType::BiCGStab:
			os << "BiCGStab";
			break;

		default:
			os << "unknown";
	}
	return os;
}

AppmSolver::SolverParameters::SolverParameters()
{
}

AppmSolver::SolverParameters::~SolverParameters()
{
}

void AppmSolver::SolverParameters::setMaxIterations(const int n)
{
	assert(n >= 0);
	this->maxIterations = n;
}

const int AppmSolver::SolverParameters::getMaxIterations() const
{
	return this->maxIterations;
}

void AppmSolver::SolverParameters::setMaxTime(const double tmax)
{
	assert(tmax >= 0);
	this->maxTime = tmax;
}

const double AppmSolver::SolverParameters::getMaxTime() const
{
	return this->maxTime;
}

const bool AppmSolver::SolverParameters::getMaxwellEnabled() const
{
	return this->isMaxwellEnabled;
}

const double AppmSolver::SolverParameters::getMaxTimestepSize() const
{
	return this->timestepSizeMax;
}

void AppmSolver::SolverParameters::setFluidEnabled(const bool b)
{
	this->isFluidEnabled = b;
}

const bool AppmSolver::SolverParameters::getFluidEnabled() const
{
	return this->isFluidEnabled;
}

void AppmSolver::SolverParameters::setLorentzForceEnabled(const bool b)
{
	this->isLorentzForceEnabled = b;
}

const bool AppmSolver::SolverParameters::getLorentzForceEnabled() const
{
	return this->isLorentzForceEnabled;
}

void AppmSolver::SolverParameters::setMassfluxSchemeImplicit(const bool b)
{
	this->isMassFluxSchemeImplicit = b;
}

const bool AppmSolver::SolverParameters::getMassfluxSchemeImplicit() const
{
	return this->isMassFluxSchemeImplicit;
}

const bool AppmSolver::SolverParameters::getFrictionActive() const
{
	return false;
}

const bool AppmSolver::SolverParameters::getMaxwellCurrentDefined() const
{
	return false;
}

const AppmSolver::MaxwellSolverType AppmSolver::SolverParameters::getMaxwellSolverType() const
{
	return AppmSolver::MaxwellSolverType::BiCGStab;
}

const bool AppmSolver::SolverParameters::getEulerMaxwellCouplingEnabled() const
{
	return false;
}

void AppmSolver::SolverParameters::setOutputFrequency(const int n)
{
	assert(n >= 0);
	this->outputFrequency = n;
}

const int AppmSolver::SolverParameters::getOutputFrequency() const
{
	assert(this->outputFrequency >= 0);
	return this->outputFrequency;
}

void AppmSolver::SolverParameters::setFluidInitType(const std::string & s)
{
	assert(s.size() > 0);
	const std::string trimmed = trim(s);
	this->fluidInitType = FluidInitType::DEFAULT;
	if (trimmed == "SHOCKTUBE") {
		this->fluidInitType = FluidInitType::SHOCKTUBE;
	}
	if (trimmed == "UNIFORM") {
		this->fluidInitType = FluidInitType::UNIFORM;
	}
	if (trimmed == "TEST_FRICTION") {
		this->fluidInitType = FluidInitType::TEST_FRICTION;
	}
	if (trimmed == "TEST_FRICTION_TEMPERATURE") {
		this->fluidInitType = FluidInitType::TEST_FRICTION_TEMPERATURE;
	}
	if (trimmed == "TEST_FRICTION_NUMBERDENSITY") {
		this->fluidInitType = FluidInitType::TEST_FRICTION_NUMBERDENSITY;
	}
	if (trimmed == "TEST_FRICTION_ELECTRONS_NONZERO_VELOCITY") {
		this->fluidInitType = FluidInitType::TEST_FRICTION_ELECTRONS_NONZERO_VELOCITY;
	}
	if (trimmed == "INIT_FILE") {
		this->fluidInitType = FluidInitType::INIT_FILE;
	}
	if (this->fluidInitType == FluidInitType::DEFAULT) {
		std::cout << "Warning: it may be that the FluidInitType is not correctly read." << std::endl;
		std::cout << "Value read from imput file: " << trimmed << std::endl;
		std::cout << "Value in solver parameters: " << this->fluidInitType << std::endl;
	}
}

const AppmSolver::FluidInitType AppmSolver::SolverParameters::getFluidInitType() const
{
	return this->fluidInitType;
}

void AppmSolver::SolverParameters::setInitEfield(const Eigen::Vector3d & efield)
{
	assert(efield.allFinite());
	this->initEfield = efield;
}

const Eigen::Vector3d AppmSolver::SolverParameters::getInitEfield() const
{
	return this->initEfield;
}

void AppmSolver::SolverParameters::setApParameter(const double lambdaSquare)
{
	assert(lambdaSquare > 0);
	assert(lambdaSquare <= 1.0);
	this->lambdaSq = lambdaSquare;
}

const double AppmSolver::SolverParameters::getApParameter() const
{
	return this->lambdaSq;
}

void AppmSolver::SolverParameters::setEulerSourcesImplicit(const bool b)
{
	this->isEulerSourcesImplicit = b;
}

const bool AppmSolver::SolverParameters::getEulerSourcesImplicit() const
{
	return isEulerSourcesImplicit;
}

const double AppmSolver::SolverParameters::getTimestepSizeFactor() const
{
	return this->timestepSizeFactor;
}

void AppmSolver::SolverParameters::setTimestepSizeFactor(const double dt_factor)
{
	assert(dt_factor > 0);
	this->timestepSizeFactor = dt_factor;
}

std::ostream & operator<<(std::ostream & os, const AppmSolver::SolverParameters & obj)
{
	os << "APPM Solver Parameters:" << std::endl;
	os << "=======================" << std::endl;
	os << "maxIterations: " << obj.maxIterations << std::endl;
	os << "maxTime: " << obj.maxTime << std::endl;
	os << "timestepSizeFactor: " << obj.timestepSizeFactor << std::endl;
	os << "timestepSizeMax: " << obj.timestepSizeMax << std::endl;
	os << "outputFrequency: " << obj.outputFrequency << std::endl;
	os << "isEulerMaxwellCouplingEnabled: " << obj.isEulerMaxwellCouplingEnabled << std::endl;
	os << "isFluidEnabled: " << obj.isFluidEnabled << std::endl;
	os << "fluidInitState: " << obj.fluidInitType << std::endl;
	os << "initEfield: " << obj.initEfield.transpose() << std::endl;
	os << "isMassfluxSchemeImplicit: " << obj.isMassFluxSchemeImplicit << std::endl;
	os << "isFrictionEnabled: " << obj.isFrictionEnabled << std::endl;
	os << "isMaxwellEnabled: " << obj.isMaxwellEnabled << std::endl;
	os << "maxwellSolverType: " << obj.maxwellSolverType << std::endl;
	os << "AP parameter: " << obj.lambdaSq << std::endl;
	os << "isEulerSourcesImplicit: " << obj.isEulerSourcesImplicit << std::endl;
	os << "=======================" << std::endl;
	return os;
}

std::ostream & operator<<(std::ostream & os, const AppmSolver::FluidInitType & obj)
{
	switch (obj) {
	case AppmSolver::FluidInitType::DEFAULT:
		os << "DEFAULT";
		break;

	case AppmSolver::FluidInitType::UNIFORM:
		os << "UNIFORM";
		break;

	case AppmSolver::FluidInitType::SHOCKTUBE:
		os << "SHOCKTUBE";
		break;

	case AppmSolver::FluidInitType::TEST_FRICTION:
		os << "TEST_FRICTION";
		break;

	case AppmSolver::FluidInitType::TEST_FRICTION_TEMPERATURE:
		os << "TEST_FRICTION_TEMPERATURE";
		break;

	case AppmSolver::FluidInitType::TEST_FRICTION_NUMBERDENSITY:
		os << "TEST_FRICTION_NUMBERDENSITY";
		break;

	case AppmSolver::FluidInitType::TEST_FRICTION_ELECTRONS_NONZERO_VELOCITY:
		os << "TEST_FRICTION_ELECTRONS_NONZERO_VELOCITY";
		break;

	case AppmSolver::FluidInitType::INIT_FILE:
		os << "INIT_FILE";
		break;

	default:
		os << "unknown";
		assert(false);
	}
	return os;
}

std::ostream & operator<<(std::ostream & os, const AppmSolver::InitDataStruct & obj)
{
	std::cout << "fluid name:\t" << obj.fluidName << std::endl;
	std::cout << "n:\t" << obj.n << std::endl;
	std::cout << "T:\t" << obj.T << std::endl;
	std::cout << "uz:\t" << obj.uz << std::endl;
	return os;
}