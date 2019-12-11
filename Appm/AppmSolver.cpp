#include "AppmSolver.h"



AppmSolver::AppmSolver()
{
}


AppmSolver::~AppmSolver()
{
}

void AppmSolver::run()
{
	double dt = 1;
	double time = 0;
	int iteration = 0;

	init_meshes();  // Initialize primal and dual meshes
	init_maxwell(dt); // Initialize Maxwell equations and states
	init_fluid(); 	// Initialize fluid states

	writeOutput(iteration, time);

	// Time integration loop
	//double dT = 0.05;
	const int maxIteration = 10;
	const double maxTime = 10;

	while (iteration < maxIteration && time < maxTime) {
		std::cout << "Iteration " << iteration << ",\t time = " << time << std::endl;
		// Fluid equations
		// dt = update_fluid();
		std::cout << "dt = " << dt << std::endl;
		
		// Maxwell equations
		update_maxwell(dt, time);


		//assert(bvec.array().isFinite().all());
		//assert(dvec.array().isFinite().all());
		//assert(evec.array().isFinite().all());
		//assert(hvec.array().isFinite().all());

		iteration++;
		time += dt;
		writeOutput(iteration, time);
	}
	std::cout << "Final time:      " << time << std::endl;
	std::cout << "Final iteration: " << iteration << std::endl;

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
}

void AppmSolver::init_meshes()
{
	std::cout << "Init primal mesh" << std::endl;
	primalMesh = PrimalMesh();
	primalMesh.init();
	primalMesh.writeToFile();
	primalMesh.writeXdmf();
	primalMesh.check();

	std::cout << "Init dual mesh" << std::endl;
	dualMesh = DualMesh();
	dualMesh.init_dualMesh(primalMesh);
	dualMesh.writeToFile();
	dualMesh.writeXdmf();
}

void AppmSolver::init_fluid()
{
	fluidStates = Eigen::MatrixXd::Zero(5, dualMesh.getNumberOfCells());
	fluidFluxes = Eigen::MatrixXd::Zero(5, dualMesh.getNumberOfFaces());

	const double rho_L = 1.0;
	const double p_L = 1.0;
	const double rho_R = 0.125;
	const double p_R = 0.1;
	const Eigen::Vector3d uZero = Eigen::Vector3d::Zero();

	FluidPrimitiveState primitiveL;
	primitiveL.rho = rho_L;
	primitiveL.p = p_L;
	primitiveL.u = uZero;

	FluidPrimitiveState primitiveR;
	primitiveR.rho = rho_R;
	primitiveR.p = p_R;
	primitiveR.u = uZero;

	const FluidState qL(primitiveL);
	const FluidState qR(primitiveR);

	const Eigen::VectorXd qL_vec = qL.getVecState();
	const Eigen::VectorXd qR_vec = qR.getVecState();

	std::cout << "qL: " << qL_vec.transpose() << std::endl;
	std::cout << "qR: " << qR_vec.transpose() << std::endl;
	std::cout << std::endl;

	for (int i = 0; i < dualMesh.getNumberOfCells(); i++) {
		const Cell * cell = dualMesh.getCell(i);
		const double idxC = cell->getIndex();
		const Eigen::Vector3d cellCenter = cell->getCenter();

		fluidStates.col(i) = (cellCenter(2) < 0.5) ? qL_vec : qR_vec;
	}

	std::cout << "Test for fluid flux: " << std::endl;
	double dt_loc = 0;
	Eigen::VectorXd fluidFlux = Numerics::fluidFlux_rusanov(qL_vec, qR_vec, Eigen::Vector3d(0, 0, 1), 1, dt_loc);
	std::cout << "Flux: " << fluidFlux.transpose() << std::endl;
	std::cout << "dt_loc: " << dt_loc << std::endl;
}

void AppmSolver::init_maxwell(const double dt)
{
	// Define data vectors
	B_h = Eigen::VectorXd::Zero(primalMesh.getNumberOfFaces());
	E_h = Eigen::VectorXd::Zero(primalMesh.getNumberOfEdges());
	J_h = Eigen::VectorXd::Zero(dualMesh.getNumberOfFaces());


	const Eigen::VectorXi vertexType = primalMesh.getVertexTypes();
	const Eigen::VectorXi edgeType = primalMesh.getEdgeTypes();
	const Eigen::VectorXi faceType = primalMesh.getFaceTypes();

	// Number of primal vertices
	const int nPv = primalMesh.getNumberOfVertices();
	// Number of primal vertices on boundary
	const int nPvb = (vertexType.array() == static_cast<int>(Vertex::Type::Boundary)).count()
		+ (vertexType.array() == static_cast<int>(Vertex::Type::Terminal)).count();
	// Number of primal vertices in interior
	const int nPvi = (vertexType.array() == static_cast<int>(Vertex::Type::Inner)).count();
	// Number of primal vertrices at terminals
	const int nPvt = (vertexType.array() == static_cast<int>(Vertex::Type::Terminal)).count();
	assert((nPvb + nPvi) == nPv);

	// Number of primal edges
	const int nPe = primalMesh.getNumberOfEdges();

	// Number of primal edges in interior
	const int nPei = (edgeType.array() == static_cast<int>(Edge::Type::Interior)).count()
		+ (edgeType.array() == static_cast<int>(Edge::Type::InteriorToBoundary)).count();

	// Number of primal faces
	const int nPf = primalMesh.getNumberOfFaces();
	const int nPfi = (faceType.array() == 0).count();

	const int nDof = nPei + nPvb;
	x_m = Eigen::VectorXd::Zero(nDof);
	x_mm1 = Eigen::VectorXd::Zero(nDof);

	// Setup of operator Q = [id, G_b]
	Q = setupOperatorQ().cast<double>();

	Eigen::SparseMatrix<double> Meps;
	Eigen::SparseMatrix<double> Mnu;
	Meps = setupOperatorMeps();
	Mnu = setupOperatorMnu().topLeftCorner(nPfi, nPfi);

	// Curl operator on primal mesh
	Eigen::SparseMatrix<int> C;
	C = primalMesh.get_f2eMap();
	this->C = C.cast<double>();

	// Curl operator for inner-faces and inner-edges
	Eigen::SparseMatrix<int> Cii = C.topLeftCorner(nPfi, nPei);

	// Setup of operator P = [Cii 0]
	Eigen::SparseMatrix<double> P;
	P = Eigen::SparseMatrix<double>(nPfi, nPei + nPvb);
	P.leftCols(Cii.cols()) = Cii.cast<double>();

	Eigen::SparseMatrix<double> B;
	A = Q.transpose() * Meps * Q;
	B = P.transpose() * Mnu * P;
	assert(A.rows() == B.rows());
	assert(A.cols() == B.cols());

	// System of equations
	assert(dt > 0);
	assert(isfinite(dt));
	this->M = A + pow(dt, 2) * B;
	const int nFreeIdx = M.rows() - nPvt;
	assert(nFreeIdx > 0);
	this->M_d = M.rightCols(nPvt);
	M_f = M.topLeftCorner(nFreeIdx, nFreeIdx);

	maxwellSolver.compute(M_f);
	if (maxwellSolver.info() != Eigen::Success) {
		std::cout << "Solver initialization failed" << std::endl;
		exit(-1);
	}
}

const double AppmSolver::update_fluid()
{
	const int nDualFaces = dualMesh.getNumberOfFaces();
	Eigen::VectorXd dt_local(nDualFaces);
	Eigen::MatrixXd dq = Eigen::MatrixXd::Zero(fluidStates.rows(), fluidStates.cols());

	// Compute fluid fluxes at faces
	for (int i = 0; i < nDualFaces; i++) {
		const Face * face = dualMesh.getFace(i);
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
			const Eigen::VectorXd flux = Numerics::fluidFlux_rusanov(qL, qR, faceNormal, dx, dt_loc);
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
				leftCell  = faceCells[0];
				rightCell = faceCells[1];
			}
			else {
				leftCell  = faceCells[1];
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
			const Eigen::VectorXd faceFlux = Numerics::fluidFlux_rusanov(qL, qR, faceNormal, dx, dt_loc);
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

	for (int i = 0; i < dualMesh.getNumberOfCells(); i++) {
		dq.col(i) *= 1./(dualMesh.getCell(i)->getVolume());
	}

	fluidStates -= dt * dq;

	return dt;
}

void AppmSolver::update_maxwell(const double dt, const double time)
{
	const int nDof = x_m.size();

	// Vector of degrees of freedom (dof)
	Eigen::VectorXd x(nDof);
	Eigen::VectorXd rhs = Eigen::VectorXd::Zero(x.size());

	rhs = A * (2 * x_m - x_mm1) + pow(dt, 2) * rhs;

	const int nPvt = M_d.cols();
	const int nFreeIdx = M.rows() - nPvt;

	// Electric potential at terminal vertices
	Eigen::VectorXd phi_t(nPvt);
	double phi1 = 0;
	double phi2 = 0.5 * (1 + tanh((time - 5)));
	phi_t.topRows(nPvt / 2).array() = phi1;
	phi_t.bottomRows(nPvt / 2).array() = phi2;

	Eigen::VectorXd x_d;
	x_d = phi_t;
	rhs -= M_d * x_d;

	Eigen::VectorXd rhs_f = rhs.topRows(nFreeIdx);
	Eigen::VectorXd x_f(nFreeIdx);
	x_f = maxwellSolver.solve(rhs_f);
	if (maxwellSolver.info() != Eigen::Success) {
		std::cout << "Solver solution failed" << std::endl;
		exit(-1);
	}
	x.topRows(x_f.size()) = x_f;
	x.bottomRows(x_d.size()) = x_d;

	// Electric field (from electrostatic boundary conditions)
	E_h = Q * x;

	// Update of magnetic field
	B_h = B_h - dt * C * E_h;

	// Update state vectors
	x_mm1 = x_m;
	x_m = x;
}


void AppmSolver::writeXdmf()
{
	const int nTimesteps = timeStamps.size();
	const std::string filename = "appm.xdmf";
	XdmfRoot root;
	XdmfDomain domain;
	XdmfGrid timeGrid(
		XdmfGrid::Tags(
			"Time Grid", 
			XdmfGrid::GridType::Collection, 
			XdmfGrid::CollectionType::Temporal
		)
	);
	for (int iteration = 0; iteration < nTimesteps; iteration++) {
		const std::string datafilename = (std::stringstream() << "appm-" << iteration << ".h5").str();
		XdmfGrid gridOfGrids(XdmfGrid::Tags("Grid of Grids", XdmfGrid::GridType::Tree));

		const double timeValue = this->timeStamps[iteration];
		gridOfGrids.addChild(XdmfTime(timeValue));

		// Primal edge grid
		XdmfGrid primalEdgeGrid = getOutputPrimalEdgeGrid(iteration, timeValue, datafilename);
		gridOfGrids.addChild(primalEdgeGrid);

		// Primal surface grid
		XdmfGrid primalSurfaceGrid = getOutputPrimalSurfaceGrid(iteration, timeValue, datafilename);
		gridOfGrids.addChild(primalSurfaceGrid);

		// Dual edge grid
		XdmfGrid dualEdgeGrid = getOutputDualEdgeGrid(iteration, timeValue, datafilename);
		gridOfGrids.addChild(dualEdgeGrid);

		// Dual surface grid
		XdmfGrid dualSurfaceGrid = getOutputDualSurfaceGrid(iteration, timeValue, datafilename);
		gridOfGrids.addChild(dualSurfaceGrid);

		//// Dual volume grid
		//XdmfGrid dualVolumeGrid;

		timeGrid.addChild(gridOfGrids);
	}
	domain.addChild(timeGrid);
	root.addChild(domain);
	std::ofstream file(filename);
	file << root << std::endl;
}

void AppmSolver::writeXdmfDualVolume()
{
	const int nTimesteps = timeStamps.size();
	const std::string filename = "appm-volume.xdmf";
	XdmfRoot root;
	XdmfDomain domain;
	XdmfGrid timeGrid(
		XdmfGrid::Tags(
			"Time Grid",
			XdmfGrid::GridType::Collection,
			XdmfGrid::CollectionType::Temporal
		)
	);
	for (int iteration = 0; iteration < nTimesteps; iteration++) {
		const std::string datafilename = (std::stringstream() << "appm-" << iteration << ".h5").str();

		const double timeValue = this->timeStamps[iteration];
		timeGrid.addChild(XdmfTime(timeValue));

		// Dual volume grid
		XdmfGrid dualVolumeGrid;
		dualVolumeGrid = getOutputDualVolumeGrid(iteration, timeValue, datafilename);

		timeGrid.addChild(dualVolumeGrid);
	}
	domain.addChild(timeGrid);
	root.addChild(domain);
	std::ofstream file(filename);
	file << root << std::endl;
}


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
	Eigen::VectorXd q_density = fluidStates.row(0);
	Eigen::MatrixXd q_momentum = fluidStates.middleRows(1,3);
	Eigen::VectorXd q_energy = fluidStates.row(4);
	h5writer.writeData(q_density,  "/qDensity");
	h5writer.writeData(q_momentum, "/qMomentum");
	h5writer.writeData(q_energy,   "/qEnergy");

	// Fluid states in primitive variables
	Eigen::VectorXd pressure(nDualCells);
	Eigen::VectorXd density(nDualCells);
	Eigen::Matrix3Xd velocity(3, nDualCells);
	for (int i = 0; i < nDualCells; i++) {
		const FluidPrimitiveState primitiveState = FluidState(fluidStates.col(i)).getPrimitiveState();
		pressure(i) = primitiveState.p;
		density(i) = primitiveState.rho;
		velocity.col(i) = primitiveState.u;
	}
	h5writer.writeData(pressure, "/pressure");
	h5writer.writeData(density, "/density");
	h5writer.writeData(velocity, "/velocity");




	// Maxwell states
	h5writer.writeData(B_h, "/bvec");
	h5writer.writeData(E_h, "/evec");
	h5writer.writeData(J_h, "/jvec");

	Eigen::Matrix3Xd B(3, nPrimalFaces);
	Eigen::Matrix3Xd J(3, nDualFaces);
	
	for (int i = 0; i < nPrimalFaces; i++) {
		const Eigen::Vector3d fn = primalMesh.getFace(i)->getNormal();
		const double fA = primalMesh.getFace(i)->getArea();
		B.col(i) = (B_h(i) / fA) * fn;
	}
	for (int i = 0; i < nDualFaces; i++) {
		const Eigen::Vector3d fn = dualMesh.getFace(i)->getNormal();
		const double fA = dualMesh.getFace(i)->getArea();
		J.col(i) = J_h(i) / fA * fn;
	}
	h5writer.writeData(B, "/B");
	h5writer.writeData(J, "/J");
	
	const int nPrimalEdges = primalMesh.getNumberOfEdges();
	Eigen::Matrix3Xd E(3, nPrimalEdges);
	for (int i = 0; i < nPrimalEdges; i++) {
		const Edge * edge = primalMesh.getEdge(i);
		E.col(i) = E_h(i) / edge->getLength() * edge->getDirection();
	}
	h5writer.writeData(E, "/E");

	Eigen::VectorXd timeVec(1);
	timeVec(0) = time;
	h5writer.writeData(timeVec, "/time");
	Eigen::VectorXi iterVec(1);
	iterVec(0) = iteration;
	h5writer.writeData(iterVec, "/iteration");
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

	//// Attribute: Magnetic Field H
	//{
	//	XdmfAttribute attribute(
	//		XdmfAttribute::Tags("Magnetic field", XdmfAttribute::Type::Vector, XdmfAttribute::Center::Cell)
	//	);
	//	attribute.addChild(
	//		XdmfDataItem(XdmfDataItem::Tags(
	//			{ dualMesh.getNumberOfEdges(), 3 },
	//			XdmfDataItem::NumberType::Float,
	//			XdmfDataItem::Format::HDF),
	//			(std::stringstream() << dataFilename << ":/H").str()
	//		));
	//	grid.addChild(attribute);
	//}
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

Eigen::SparseMatrix<int> AppmSolver::setupOperatorQ()
{
	const Eigen::VectorXi vertexTypes = primalMesh.getVertexTypes();
	const Eigen::VectorXi edgeTypes = primalMesh.getEdgeTypes();

	const int boundaryVertexType = static_cast<int>(Vertex::Type::Boundary);
	const int terminalVertexType = static_cast<int>(Vertex::Type::Terminal);
	const int innerVertexType = static_cast<int>(Vertex::Type::Inner);
	
	const int boundaryEdgeType = static_cast<int>(Edge::Type::Boundary);
	const int interiorToBoundaryEdgeType = static_cast<int>(Edge::Type::InteriorToBoundary);
	const int interiorEdgeType = static_cast<int>(Edge::Type::Interior);

	const int nPv = primalMesh.getNumberOfVertices();
	const int nPvb = (vertexTypes.array() == boundaryVertexType).count() + (vertexTypes.array() == terminalVertexType).count();
	const int nPvi = (vertexTypes.array() == innerVertexType).count();
	const int nPe = primalMesh.getNumberOfEdges();
	const int nPei = (edgeTypes.array() == interiorEdgeType).count() + (edgeTypes.array() == interiorToBoundaryEdgeType).count();
	const int nPeb = (edgeTypes.array() == boundaryEdgeType).count();

	assert(nPe == nPei + nPeb);

	Eigen::SparseMatrix<int> X(nPv, nPvb);
	typedef Eigen::Triplet<int> T;
	std::vector<T> triplets;
	for (int i = 0; i < nPvb; i++) {
		triplets.push_back(T(nPvi + i, i, 1));
	}
	X.setFromTriplets(triplets.begin(), triplets.end());
	X.makeCompressed();

	const Eigen::SparseMatrix<int> G = primalMesh.get_e2vMap();
	Eigen::SparseMatrix<int> GX = G * X;

	assert(nPei < nPe);
	Eigen::SparseMatrix<int> id(nPe, nPei);
	//id.setIdentity(); // only for square matrices
	for (int i = 0; i < nPei; i++) {
		id.coeffRef(i, i) = 1;
	}

	// Number of degrees of freedom
	const int nDof = id.cols() + GX.cols();
	assert(nDof == id.cols() + GX.cols());
	Eigen::SparseMatrix<int> Q = Eigen::SparseMatrix<int>(nPe, nDof);
	Q.leftCols(id.cols()) = id;
	Q.rightCols(GX.cols()) = GX;
	assert(Q.rows() == nPe);
	assert(Q.cols() == nDof);
	return Q;
}

Eigen::SparseMatrix<double> AppmSolver::setupOperatorMeps()
{
	const int nPe = primalMesh.getNumberOfEdges();
	Eigen::SparseMatrix<double> Meps(nPe, nPe);
	Meps.setIdentity();
	for (int i = 0; i < nPe; i++) {
		const double dualFaceArea = dualMesh.getFace(i)->getArea();
		const double primalEdgeLength = primalMesh.getEdge(i)->getLength();
		Meps.coeffRef(i, i) = dualFaceArea / primalEdgeLength;
	}
	return Meps;
}

Eigen::SparseMatrix<double> AppmSolver::setupOperatorMnu()
{
	const int nPf = primalMesh.getNumberOfFaces();
	Eigen::SparseMatrix<double> Mnu(nPf, nPf);
	Mnu.setIdentity();
	for (int i = 0; i < nPf; i++) {
		const double dualEdgeLength = dualMesh.getEdge(i)->getLength();
		const double primalFaceArea = primalMesh.getFace(i)->getArea();
		Mnu.coeffRef(i, i) = dualEdgeLength / primalFaceArea;
	}
	return Mnu;
}

