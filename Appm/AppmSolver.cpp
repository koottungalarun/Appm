#include "AppmSolver.h"



AppmSolver::AppmSolver()
{
}


AppmSolver::~AppmSolver()
{
}

void AppmSolver::run()
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

	// Define data vectors
	bvec = Eigen::VectorXd::Zero(primalMesh.getNumberOfFaces());
	dvec = Eigen::VectorXd::Zero(dualMesh.getNumberOfFaces());
	evec = Eigen::VectorXd::Zero(primalMesh.getNumberOfEdges());
	hvec = Eigen::VectorXd::Zero(dualMesh.getNumberOfEdges());
	jvec = Eigen::VectorXd::Zero(dualMesh.getNumberOfFaces());

	// TODO: fluid state vector

	const int nDualFaces = dualMesh.getNumberOfFaces();
	for (int i = 0; i < nDualFaces; i++) {
		const Face * face = dualMesh.getFace(i);

		const Eigen::Vector3d fc = face->getCenter();
		const Eigen::Vector2d fc_2d = fc.segment(0, 2);

		if (fc_2d.norm() < 0.3) {
			const Eigen::Vector3d fn = face->getNormal();
			jvec(i) = 1 * face->getArea() * fn.dot(Eigen::Vector3d(0, 0, 1));
		}
	}
	double time = 0;
	int iteration = 0;

	for (int i = 0; i < primalMesh.getNumberOfEdges(); i++) {
		const Edge * edge = primalMesh.getEdge(i);
		const Eigen::Vector3d edgeDir = edge->getDirection();
		evec(i) = edgeDir.dot(Eigen::Vector3d(0, 0, 1));
	}
	for (int i = 0; i < primalMesh.getNumberOfFaces(); i++) {
		bvec(i) = primalMesh.getFace(i)->getNormal().dot(Eigen::Vector3d(0, 0, 1));
	}
	for (int i = 0; i < dualMesh.getNumberOfFaces(); i++) {
		dvec(i) = dualMesh.getFace(i)->getNormal().dot(Eigen::Vector3d(0, 0, 1));
	}
	for (int i = 0; i < dualMesh.getNumberOfEdges(); i++) {
		const Edge * edge = dualMesh.getEdge(i);
		const Eigen::Vector3d edgeDir = edge->getDirection();
		hvec(i) = edgeDir.dot(Eigen::Vector3d(0, 0, 1));
	}


	writeOutput(iteration, time);

	const Eigen::SparseMatrix<double> Cprimal = primalMesh.get_f2eMap().cast<double>();
	const Eigen::SparseMatrix<double> Cdual = dualMesh.get_f2eMap().cast<double>();
	assert(Cprimal.rows() == bvec.size());
	assert(Cprimal.cols() == evec.size());
	assert(Cdual.rows() == dvec.size());
	assert(Cdual.cols() == hvec.size());

	const int nPrimalFaces = primalMesh.getNumberOfFaces();
	const int nPrimalEdges = primalMesh.getNumberOfEdges();
	typedef Eigen::Triplet<double> T;
	std::vector<T> triplets;
	for (int i = 0; i < nPrimalFaces; i++) {
		double value = dualMesh.getEdge(i)->getLength() / primalMesh.getFace(i)->getArea();
		triplets.push_back(T(i, i, value));
	}
	Eigen::SparseMatrix<double> hodgeOp_muInv(nPrimalFaces, nPrimalFaces);
	hodgeOp_muInv.setFromTriplets(triplets.begin(), triplets.end());
	hodgeOp_muInv.makeCompressed();

	triplets.resize(0);
	for (int i = 0; i < nPrimalEdges; i++) {
		double value = primalMesh.getEdge(i)->getLength() / dualMesh.getFace(i)->getArea();
		triplets.push_back(T(i, i, value));
	}
	Eigen::SparseMatrix<double> hodgeOp_epsInv(nPrimalEdges, nPrimalEdges);
	hodgeOp_epsInv.setFromTriplets(triplets.begin(), triplets.end());
	hodgeOp_epsInv.makeCompressed();

	// Time integration loop
	double dT = 0.05;
	const int maxIteration = 0;
	const double maxTime = 10;
	while (iteration < maxIteration && time < maxTime) {
		//bvec = bvec - dT * Cprimal * evec;
		//hvec.topRows(hodgeOp_muInv.rows()) = hodgeOp_muInv * bvec;
		//dvec = dvec + dT * Cdual * hvec - dT * jvec;
		//evec = hodgeOp_epsInv * dvec.topRows(hodgeOp_epsInv.cols());

		assert(bvec.array().isFinite().all());
		assert(dvec.array().isFinite().all());
		assert(evec.array().isFinite().all());
		assert(hvec.array().isFinite().all());

		iteration++;
		time += dT;
		writeOutput(iteration, time);
	}
	std::cout << "Final time:      " << time << std::endl;
	std::cout << "Final iteration: " << iteration << std::endl;

	writeXdmf();
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


void AppmSolver::writeOutput(const int iteration, const double time)
{
	timeStamps.push_back(time);
	std::cout << "Write output at iteration " << iteration << ", time = " << time << std::endl;
	const std::string filename = (std::stringstream() << "appm-" << iteration << ".h5").str();
	H5Writer h5writer(filename);
	h5writer.writeData(bvec, "/bvec");
	h5writer.writeData(dvec, "/dvec");
	h5writer.writeData(evec, "/evec");
	h5writer.writeData(hvec, "/hvec");
	h5writer.writeData(jvec, "/jvec");

	const int nPrimalFaces = primalMesh.getNumberOfFaces();
	const int nDualEdges = dualMesh.getNumberOfEdges();
	const int nDualFaces = dualMesh.getNumberOfFaces();
	Eigen::Matrix3Xd B(3, nPrimalFaces);
	Eigen::Matrix3Xd D(3, nDualFaces);
	Eigen::Matrix3Xd J(3, nDualFaces);
	Eigen::Matrix3Xd H(3, nDualEdges);

	for (int i = 0; i < nPrimalFaces; i++) {
		const Eigen::Vector3d fn = primalMesh.getFace(i)->getNormal();
		B.col(i) = bvec(i) * fn;
	}
	for (int i = 0; i < nDualEdges; i++) {
		const Edge * edge = dualMesh.getEdge(i);
		H.col(i) = hvec(i) / edge->getLength() * edge->getDirection();
	}
	for (int i = 0; i < nDualFaces; i++) {
		const Eigen::Vector3d fn = dualMesh.getFace(i)->getNormal();
		D.col(i) = dvec(i) * fn;
		J.col(i) = jvec(i) * fn;
	}
	h5writer.writeData(B, "/B");
	h5writer.writeData(D, "/D");
	h5writer.writeData(J, "/J");
	h5writer.writeData(H, "/H");

	const int nPrimalEdges = primalMesh.getNumberOfEdges();
	Eigen::Matrix3Xd E(3, nPrimalEdges);
	for (int i = 0; i < nPrimalEdges; i++) {
		const Edge * edge = primalMesh.getEdge(i);
		E.col(i) = evec(i) / edge->getLength() * edge->getDirection();
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

	// Attribute: Displacement Field D
	{
		XdmfAttribute attribute(
			XdmfAttribute::Tags("Displacement Field", XdmfAttribute::Type::Vector, XdmfAttribute::Center::Cell)
		);
		attribute.addChild(
			XdmfDataItem(XdmfDataItem::Tags(
				{ dualMesh.getNumberOfFaces(), 3 },
				XdmfDataItem::NumberType::Float,
				XdmfDataItem::Format::HDF),
				(std::stringstream() << dataFilename << ":/D").str()
			));
		grid.addChild(attribute);
	}

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

