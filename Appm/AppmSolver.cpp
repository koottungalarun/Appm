#include "AppmSolver.h"



AppmSolver::AppmSolver()
{
}


AppmSolver::~AppmSolver()
{
}

void AppmSolver::run()
{
	double dt = 0.1;
	double time = 0;
	int iteration = 0;

	init_meshes();  // Initialize primal and dual meshes
	init_maxwell(); // Initialize Maxwell equations and states
	init_fluid(); 	// Initialize fluid states

	test_raviartThomas();

	writeOutput(iteration, time);

	// Time integration loop
	//double dT = 0.05;
	const int maxIteration = 0;
	const double maxTime = 20;

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

AppmSolver::MeshInfo AppmSolver::setMeshInfo(const Mesh & mesh)
{
	MeshInfo meshInfo;

	// Vertices
	const Eigen::VectorXi vertexTypes = mesh.getVertexTypes();
	const int boundaryVertexType = static_cast<int>(Vertex::Type::Boundary);
	const int terminalVertexType = static_cast<int>(Vertex::Type::Terminal);
	meshInfo.nVertices = mesh.getNumberOfVertices();
	meshInfo.nVerticesTerminal = (vertexTypes.array() == terminalVertexType).count();
	meshInfo.nVerticesBoundary = (vertexTypes.array() == boundaryVertexType).count() + meshInfo.nVerticesTerminal;

	// Edges
	const Eigen::VectorXi edgeTypes = mesh.getEdgeTypes();
	const int interiorEdgeType = static_cast<int>(Edge::Type::Interior);
	const int interiorToBoundaryEdgeType = static_cast<int>(Edge::Type::InteriorToBoundary);
	const int boundaryEdgeType = static_cast<int>(Edge::Type::Boundary);
	meshInfo.nEdges = mesh.getNumberOfEdges();
	meshInfo.nEdgesInner = 
		(edgeTypes.array() == interiorEdgeType).count() + 
		(edgeTypes.array() == interiorToBoundaryEdgeType).count();

	// Faces
	const Eigen::VectorXi faceTypes = mesh.getFaceTypes();
	const int isBoundaryFace = 1;
	meshInfo.nFaces = mesh.getNumberOfFaces();
	meshInfo.nFacesInner = (faceTypes.array() != isBoundaryFace).count();

	// Cells
	meshInfo.nCells = mesh.getNumberOfCells();

	return meshInfo;
}

void AppmSolver::test_raviartThomas()
{
	std::cout << "Test Raviart Thomas Interpolation" << std::endl;
	// Magnetic flux B at primal vertices
	const int nVertices = primalMeshInfo.nVertices;
	Eigen::Matrix3Xd B_vertex(3, nVertices);
	B_vertex.setZero();

	// Set magnetic flux B at primal faces.
	const Eigen::Vector3d zUnit(0, 0, 1);
	const int nFaces = primalMeshInfo.nFaces;
	for (int i = 0; i < nFaces; i++) {
		const Face * face = primalMesh.getFace(i);
		const double fa = face->getArea();
		const Eigen::Vector3d fn = face->getNormal();
		if (fn.cross(zUnit).norm() > 0.99) {
			// azimuthal field
			Eigen::Vector3d fc = face->getCenter();
			fc(2) = 0;
			fc.normalize();
			B_h(i) = fn.dot(fc.cross(zUnit));
		}
		else {
			B_h(i) = 0;
		}
	}

	// Primal cells are prisms.
	// For each cell, 
	// - compute Piola map, 
	// - compute reference coordinates of vertices
	// - evaluate RT basis functions 
	// - map RT basis functions from reference to actual coordinates
	// - multiply by magnetic flux B
	// - write data to file: B at vertices
	const int nCells = primalMeshInfo.nCells;
	for (int cidx = 0; cidx < nCells; cidx++) {
#ifdef _RT_ONECELL
		std::cout << "cidx = " << cidx << std::endl;
#endif
		const Cell* cell = primalMesh.getCell(cidx);
		const std::vector<Face*> cellFaces = cell->getFaceList();
		assert(cellFaces.size() == 5); // prism cells
		assert(cellFaces[3]->getVertexList().size() == 3); // triangle faces are at end of list
		assert(cellFaces[4]->getVertexList().size() == 3);

		const Face * bottomFace = cellFaces[3];
		const Face * topFace = cellFaces[4];
		const std::vector<Vertex*> bottomVertices = bottomFace->getVertexList();
		const std::vector<Vertex*> topVertices = topFace->getVertexList();

		// check if vertices of triangle faces have same ordering
		for (int i = 0; i < 3; i++) {
			Eigen::Vector3d v = topVertices[i]->getPosition() - bottomVertices[i]->getPosition();
			v.normalize();
			assert(v.cross(zUnit).norm() < 10 * std::numeric_limits<double>::epsilon());
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

		std::vector<Vertex*> cellVertices(6);
		for (int i = 0; i < 3; i++) {
			cellVertices[i] = bottomVertices[i];
			cellVertices[i + 3] = topVertices[i];
		}

		const Eigen::Vector3d v0 = cellVertices[0]->getPosition();
		const Eigen::Vector3d v1 = cellVertices[1]->getPosition();
		const Eigen::Vector3d v2 = cellVertices[2]->getPosition();
		const Eigen::Vector3d v3 = cellVertices[5]->getPosition();

		// Piola map: xRef -> x := BK*xRef + bk
		Eigen::Matrix3d BK;
		BK.col(0) = v0 - v2;
		BK.col(1) = v1 - v2;
		BK.col(2) = v3 - v2;
		const double detBK = BK.determinant();
		const Eigen::Vector3d bk = v2;
#ifdef _RT_ONECELL
		std::cout << "BK: " << BK << std::endl;
		std::cout << "det(BK): " << detBK << std::endl;
#endif
		// vertex coordinates in actual space
		Eigen::Matrix3Xd vertexCoords(3, 6);
		for (int i = 0; i < 6; i++) {
			vertexCoords.col(i) = cellVertices[i]->getPosition();
		}
		// points in reference coordinates 
		const int nSamples = 5;
		Eigen::Matrix3Xd refCoords3d(3, 6);
		refCoords3d.setZero();
#ifdef _RT_ONECELL
		refCoords3d = getPrismReferenceCoords(nSamples);
#else
		// Reference coordinates of cell vertices
		for (int i = 0; i < 6; i++) {
			const Eigen::Vector3d pos = vertexCoords.col(i);
			const Eigen::Vector3d refPos = BK.inverse() * (pos - bk);
			refCoords3d.col(i) = refPos;
		}
#endif

		assert((refCoords3d.array() >= -1 * std::numeric_limits<double>::epsilon()).all());
		assert((refCoords3d.array() - 1 <= std::numeric_limits<double>::epsilon()).all());

		for (int i = 0; i < cellFaces.size(); i++) {
			const Face * face = cellFaces[i];
			const Eigen::Vector3d fn = face->getNormal();
			const Eigen::Vector3d fc = face->getCenter();

#ifdef _RT_ONECELL
			std::cout << "Face " << i << ": ";
			std::cout << "fc: " << fc.transpose() << "\t ";
			std::cout << "fn: " << fn.transpose() << "\t "; 
			std::cout << "incid: " << cell->getOrientation(face) << "\t ";
			std::cout << std::endl;
#endif
		}

		Eigen::VectorXd coeff = Eigen::VectorXd::Zero(5);
		Eigen::VectorXd B_h_loc(5);
		for (int i = 0; i < 5; i++) {
			const Face * face = cellFaces[i];
			const int fidx = face->getIndex();
			const int orientation = cell->getOrientation(face);
			B_h_loc(i) = B_h(fidx);
			coeff(i) = orientation * B_h(fidx);
		}
#ifdef _RT_ONECELL
		std::cout << "B_h_loc: " << B_h_loc << std::endl;
		std::cout << "coeff: " << coeff << std::endl;
#endif
		
		// for each vertex
		Eigen::Matrix3Xd Btest(3, refCoords3d.cols());
		Eigen::Matrix3Xd xCoords3d(3, refCoords3d.cols());
		for (int i = 0; i < refCoords3d.cols(); i++) {
			const Eigen::Vector3d refPos = refCoords3d.col(i);
			Eigen::Vector3d Bv;
			Bv.setZero();
			// for each basis function
			for (int idx = 0; idx < 5; idx++) {
				const Eigen::Vector3d rt = Numerics::raviartThomasBasis(idx, refPos);
				Bv += coeff(idx) / detBK * BK * rt;
			}
#ifdef _RT_ONECELL
			Btest.col(i) = Bv;
			xCoords3d.col(i) = BK * refPos + bk;
#else
			const int vIdx = cellVertices[i]->getIndex();
			B_vertex.col(vIdx) += Bv;
#endif
		}
#ifdef _RT_ONECELL
		H5Writer tempwriter("temp.h5");
		tempwriter.writeData(refCoords3d, "/refCoords3d");
		tempwriter.writeData(Btest, "/Btest");
		tempwriter.writeData(xCoords3d, "/xCoords3d");
		break;
#endif
	}

	H5Writer h5writer("rt.h5");
	h5writer.writeData(B_vertex, "/Bvertex");
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

void AppmSolver::init_meshes()
{
	std::cout << "Init primal mesh" << std::endl;
	primalMesh = PrimalMesh();
	primalMesh.init();
	primalMeshInfo = setMeshInfo(primalMesh);
	primalMesh.writeToFile();
	primalMesh.writeXdmf();
	primalMesh.check();


	std::cout << "Init dual mesh" << std::endl;
	dualMesh = DualMesh();
	dualMesh.init_dualMesh(primalMesh);
	dualMeshInfo = setMeshInfo(dualMesh);
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

/** Initialize Maxwell objects that are independent of timestep. */
void AppmSolver::init_maxwell()
{
	// Define data vectors
	B_h = Eigen::VectorXd::Zero(primalMesh.getNumberOfFaces());
	E_h = Eigen::VectorXd::Zero(primalMesh.getNumberOfEdges());
	H_h = Eigen::VectorXd::Zero(dualMesh.getNumberOfEdges());
	J_h = Eigen::VectorXd::Zero(dualMesh.getNumberOfFaces());
}

/** Initialize Maxwell objects. */
void AppmSolver::init_maxwell(const double dt)
{
	init_maxwell();
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
	h5writer.writeData(H_h, "/hvec");
	h5writer.writeData(J_h, "/jvec");

	Eigen::Matrix3Xd B(3, nPrimalFaces);
	Eigen::Matrix3Xd H(3, nDualEdges);
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

	for (int i = 0; i < nDualEdges; i++) {
		const Edge * edge = dualMesh.getEdge(i);
		H.col(i) = H_h(i) / edge->getLength() * edge->getDirection();
	}
	h5writer.writeData(H, "/H");

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

Eigen::VectorXd AppmSolver::electricPotentialTerminals(const double time)
{
	const int n = primalMeshInfo.nVerticesTerminal;
	assert(n % 2 == 0);
	const double t0 = 3;
	const double sigma_t = 1;
	const double phi1 = 10 * 0.5 * (1 + tanh( (time - t0) / sigma_t));
	const double phi2 = 0;
	std::cout << "phi1 = " << phi1 << std::endl;
	Eigen::VectorXd phi_terminals(n);
	phi_terminals.topRows(n / 2).array() = phi1;
	phi_terminals.bottomRows(n / 2).array() = phi2;
	return phi_terminals;
}

Eigen::SparseMatrix<double> AppmSolver::speye(const int rows, const int cols)
{
	assert(rows > 0);
	assert(cols > 0);
	const int n = std::min(rows, cols);

	typedef Eigen::Triplet<double> TripletType;
	std::vector<TripletType> triplets;
	triplets.reserve(n);
	for (int i = 0; i < n; i++) {
		triplets.emplace_back(TripletType(i, i, 1.0));
	}
	Eigen::SparseMatrix<double> M(rows, cols);
	M.setFromTriplets(triplets.begin(), triplets.end());
	return M;
}

Eigen::SparseMatrix<double> AppmSolver::hodgeOperatorPrimalEdgeToDualFace()
{
	const int nEdges = primalMeshInfo.nEdges;
	assert(nEdges > 0);

	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;
	triplets.reserve(nEdges);
	for (int i = 0; i < nEdges; i++) {
		const double fA = dualMesh.getFace(i)->getArea();
		const double eL = primalMesh.getEdge(i)->getLength();
		const double value = fA / eL;
		triplets.emplace_back(T(i, i, value));
	}
	Eigen::SparseMatrix<double> M(nEdges, nEdges);
	M.setFromTriplets(triplets.begin(), triplets.end());
	M.makeCompressed();
	return M;
}

Eigen::SparseMatrix<double> AppmSolver::hodgeOperatorDualEdgeToPrimalFace()
{
	const int nFacesInner = primalMeshInfo.nFacesInner;
	assert(nFacesInner > 0);
	
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;
	triplets.reserve(nFacesInner);
	for (int i = 0; i < nFacesInner; i++) {
		const double fA = primalMesh.getFace(i)->getArea();
		const double eL = dualMesh.getEdge(i)->getLength();
		const double value = fA / eL;
		triplets.emplace_back(T(i, i, value));
	}
	Eigen::SparseMatrix<double> M(nFacesInner, nFacesInner);
	M.setFromTriplets(triplets.begin(), triplets.end());
	M.makeCompressed();
	return M;
}

Eigen::SparseMatrix<double> AppmSolver::hodgeOperatorElectricalConductivity()
{
	const int n = primalMeshInfo.nEdges;
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;
	triplets.reserve(n);
	for (int i = 0; i < n; i++) {
		const double sigma = 1; // TODO <<<<----------------------------------------
		const double fA = dualMesh.getFace(i)->getArea();
		const double eL = primalMesh.getEdge(i)->getLength();
		const double value = fA / eL * sigma;
		triplets.emplace_back(T(i, i, value));
	}
	Eigen::SparseMatrix<double> M(n, n);
	M.setFromTriplets(triplets.begin(), triplets.end());
	M.makeCompressed();
	return M;
}

Eigen::SparseMatrix<double> AppmSolver::inclusionOperatorBoundaryVerticesToAllVertices()
{
	const int nVertices = primalMeshInfo.nVertices;
	const int nVerticesBoundary = primalMeshInfo.nVerticesBoundary;
	assert(nVerticesBoundary > 0);
	assert(nVertices > nVerticesBoundary);

	typedef Eigen::Triplet<double> TripletType;
	std::vector<TripletType> triplets;
	triplets.reserve(nVerticesBoundary);

	for (int i = 0; i < nVerticesBoundary; i++) {
		const int offset = nVertices - nVerticesBoundary;
		triplets.emplace_back(TripletType(offset + i, i, 1.0));
	}
	Eigen::SparseMatrix<double> X(nVertices, nVerticesBoundary);
	X.setFromTriplets(triplets.begin(), triplets.end());
	X.makeCompressed();
	return X;
}
