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
	return;

	// Define data vectors
	bvec = Eigen::VectorXd::Zero(primalMesh.getNumberOfFaces());
	dvec = Eigen::VectorXd::Zero(dualMesh.getNumberOfFaces());
	evec = Eigen::VectorXd::Zero(primalMesh.getNumberOfEdges());
	hvec = Eigen::VectorXd::Zero(dualMesh.getNumberOfEdges());
	jvec = Eigen::VectorXd::Zero(dualMesh.getNumberOfFaces());

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

	//for (int i = 0; i < primalMesh.getNumberOfEdges(); i++) {
	//	const Edge * edge = primalMesh.getEdge(i);
	//	const Eigen::Vector3d edgeDir = edge->getDirection();
	//	evec(i) = edgeDir.dot(Eigen::Vector3d(0, 0, 1));
	//}

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
	const int maxIteration = 1000;
	const double maxTime = 10;
	while (iteration < maxIteration && time < maxTime) {
		bvec = bvec - dT * Cprimal * evec;
		hvec.topRows(hodgeOp_muInv.rows()) = hodgeOp_muInv * bvec;
		dvec = dvec + dT * Cdual * hvec - dT * jvec;
		evec = hodgeOp_epsInv * dvec.topRows(hodgeOp_epsInv.cols());

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

	//writeXdmfPrimalMesh();
	//writeXdmfDualMesh();
	const int nTimesteps = timeStamps.size();

	const std::string filename = "appm-surface.xdmf";
	std::ofstream file(filename);
	file << "<?xml version=\"1.0\" ?>" << std::endl;
	file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;

	//XmlElement root("<Xdmf Version=\"3.0\" xmlns:xi=\"[http://www.w3.org/2001/XInclude]\">", "</Xdmf>");
	//XmlElement * domain = new XmlElement("<Domain>", "</Domain>");
	//root.addChild(domain);

	//XmlElement * timeGrid = new XmlElement("<Grid Name=\"Time Grid\" GridType=\"Collection\" CollectionType=\"Temporal\">", "</Grid>");
	//domain->addChild(timeGrid);

	//std::stringstream body;
	//std::stringstream headTag;

	//XmlElement * topology;
	//XmlElement * dataItem;
	//XmlElement * geometry;
	//XmlElement * attribute;

	//H5Reader h5reader;
	//h5reader = H5Reader("primal-mesh.h5");
	//const int nElemPrimalTopo = h5reader.readDataSize("/face2vertex");
	//assert(nElemPrimalTopo > 0);
	//
	//h5reader = H5Reader("dual-mesh.h5");
	//const int nElemDualTopo = h5reader.readDataSize("/face2vertex");
	//assert(nElemDualTopo > 0);



	//for (int i = 0; i < nTimesteps; i++) {
	//	
	//	const std::string datafilename = (std::stringstream() << "appm-" << i << ".h5").str();
	//	XmlElement * gridList = new XmlElement("<Grid Name=\"Two Grids\" GridType=\"Tree\">", "</Grid>");
	//	
	//	XmlElement * timeTag = new XmlElement((std::stringstream() << "<Time Value=\"" << timeStamps[i] << "\"/>").str());
	//	gridList->addChild(timeTag);
	//
	//	// Primal surface grid
	//	XmlElement * primalGrid = new XmlElement("<Grid Name=\"Primal\" GridType=\"Uniform\">", "</Grid>");
	//	topology = new XmlElement((std::stringstream() << "<Topology TopologyType=\"Mixed\" NumberOfElements=\"" << primalMesh.getNumberOfFaces() << "\">").str(), "</Topology>");
	//	body = std::stringstream();
	//	body << primalMesh.getPrefix() << "-mesh.h5:/face2vertex";
	//	headTag = std::stringstream();
	//	headTag << "<DataItem Dimensions=\"" << nElemPrimalTopo << "\" DataType=\"Int\" Format=\"HDF\">";
	//	dataItem = new XmlElement(headTag.str(), "</DataItem>", body.str());
	//	topology->addChild(dataItem);
	//	primalGrid->addChild(topology);

	//	geometry = new XmlElement("<Geometry GeometryType=\"XYZ\">", "</Geometry>");
	//	headTag = std::stringstream();
	//	headTag << "<DataItem Dimensions=\"" << primalMesh.getNumberOfVertices() << " 3\" DataType=\"Float\" Precision=\"8\" Format=\"HDF\">";
	//	body = std::stringstream() << primalMesh.getPrefix() << "-mesh.h5:/vertexPos";
	//	dataItem = new XmlElement(headTag.str(), "</DataItem>", body.str());
	//	geometry->addChild(dataItem);
	//	primalGrid->addChild(geometry);

	//	// Attribute: B-field
	//	attribute = new XmlElement("<Attribute Name=\"Magnetic Flux\" AttributeType=\"Vector\" Center=\"Cell\">", "</Attribute>");
	//	body = std::stringstream() << datafilename << ":" << "/B";
	//	headTag = std::stringstream() << "<DataItem Dimensions=\"" << primalMesh.getNumberOfFaces() << " 3\" DataType=\"Float\" Precision=\"8\" Format=\"HDF\">";
	//	dataItem = new XmlElement(headTag.str(), "</DataItem>", body.str());
	//	attribute->addChild(dataItem);
	//	primalGrid->addChild(attribute);

	//	// Dual surface grid
	//	XmlElement * dualGrid = new XmlElement("<Grid Name=\"Dual\" GridType=\"Uniform\">", "</Grid>");
	//	topology = new XmlElement((std::stringstream() << "<Topology TopologyType=\"Mixed\" NumberOfElements=\"" << dualMesh.getNumberOfFaces() << "\">").str(), "</Topology>");
	//	body = std::stringstream();
	//	body << dualMesh.getPrefix() << "-mesh.h5:/face2vertex";
	//	headTag = std::stringstream();
	//	headTag << "<DataItem Dimensions=\"" << nElemDualTopo << "\" DataType=\"Int\" Format=\"HDF\">";
	//	dataItem = new XmlElement(headTag.str(), "</DataItem>", body.str());
	//	topology->addChild(dataItem);
	//	dualGrid->addChild(topology);

	//	geometry = new XmlElement("<Geometry GeometryType=\"XYZ\">", "</Geometry>");
	//	headTag = std::stringstream();
	//	headTag << "<DataItem Dimensions=\"" << dualMesh.getNumberOfVertices() << " 3\" DataType=\"Float\" Precision=\"8\" Format=\"HDF\">";
	//	body = std::stringstream() << dualMesh.getPrefix() << "-mesh.h5:/vertexPos";
	//	dataItem = new XmlElement(headTag.str(), "</DataItem>", body.str());
	//	geometry->addChild(dataItem);
	//	dualGrid->addChild(geometry);

	//	// Attribute: D-field
	//	attribute = new XmlElement("<Attribute Name=\"Displacement Field\" AttributeType=\"Vector\" Center=\"Cell\">", "</Attribute>");
	//	body = std::stringstream() << datafilename << ":" << "/D";
	//	headTag = std::stringstream() << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << " 3\" DataType=\"Float\" Precision=\"8\" Format=\"HDF\">";
	//	dataItem = new XmlElement(headTag.str(), "</DataItem>", body.str());
	//	attribute->addChild(dataItem);
	//	dualGrid->addChild(attribute);

	//	// Attribute: J-field
	//	attribute = new XmlElement("<Attribute Name=\"Current\" AttributeType=\"Vector\" Center=\"Cell\">", "</Attribute>");
	//	body = std::stringstream() << datafilename << ":" << "/J";
	//	headTag = std::stringstream() << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << " 3\" DataType=\"Float\" Precision=\"8\" Format=\"HDF\">";
	//	dataItem = new XmlElement(headTag.str(), "</DataItem>", body.str());
	//	attribute->addChild(dataItem);
	//	dualGrid->addChild(attribute);

	//	// Primal edge grid
	//	XmlElement * primalEdgeGrid = new XmlElement("<Grid Name=\"PrimalEdges\" GridType=\"Uniform\">", "</Grid>");
	//	topology = new XmlElement((std::stringstream() << "<Topology TopologyType=\"Polyline\" NumberOfElements=\"" << primalMesh.getNumberOfEdges() << "\" NodesPerElement=\"2\">").str(), "</Topology>");
	//	body = std::stringstream() << primalMesh.getPrefix() << "-mesh.h5" << ":/edge2vertex";
	//	dataItem = new XmlElement((std::stringstream() << "<DataItem Dimensions=\"" << 2*primalMesh.getNumberOfEdges() << "\" DataType=\"Int\" Precision=\"4\" Format=\"HDF\">").str(), "</DataItem>", body.str());
	//	topology->addChild(dataItem);
	//	primalEdgeGrid->addChild(topology);

	//	geometry = new XmlElement("<Geometry GeometryType=\"XYZ\">", "</Geometry>");
	//	body = std::stringstream() << primalMesh.getPrefix() << "-mesh.h5" << ":/vertexPos";
	//	dataItem = new XmlElement((std::stringstream() << "<DataItem Dimensions=\"" << primalMesh.getNumberOfVertices() << " 3\" DataType=\"Float\" Precision=\"8\" Format=\"HDF\">").str(), "</DataItem>", body.str());
	//	geometry->addChild(dataItem);
	//	primalEdgeGrid->addChild(geometry);

	//	// Attribute: edge index
	//	//attribute = new XmlElement("<Attribute Name=\"EdgeIndex\" AttributeType=\"Scalar\" Center=\"Cell\">", "</Attribute>");
	//	//body = std::stringstream() << primalMesh.getPrefix() << "-mesh.h5:" << "/edgeIdx";
	//	//headTag = std::stringstream() << "<DataItem Dimensions=\"" << primalMesh.getNumberOfEdges() << "\" DataType=\"Int\" Precision=\"4\" Format=\"HDF\">";
	//	//dataItem = new XmlElement(headTag.str(), "</DataItem>", body.str());
	//	//attribute->addChild(dataItem);
	//	//primalEdgeGrid->addChild(attribute);

	//	// Attribute: E-field
	//	attribute = new XmlElement("<Attribute Name=\"Electric Field\" AttributeType=\"Vector\" Center=\"Cell\">", "</Attribute>");
	//	body = std::stringstream() << datafilename << ":" << "/E";
	//	dataItem = new XmlElement((std::stringstream() << "<DataItem Dimensions=\"" << primalMesh.getNumberOfEdges() << " 3\" DataType=\"Float\" Precision=\"8\" Format=\"HDF\">").str(), "</DataItem>", body.str());
	//	attribute->addChild(dataItem);
	//	primalEdgeGrid->addChild(attribute);

	//	attribute = new XmlElement("<Attribute Name=\"evec\" AttributeType=\"Scalar\" Center=\"Cell\">", "</Attribute>");
	//	body = std::stringstream() << datafilename << ":" << "/evec";
	//	dataItem = new XmlElement((std::stringstream() << "<DataItem Dimensions=\"" << primalMesh.getNumberOfEdges() << "\" DataType=\"Float\" Precision=\"8\" Format=\"HDF\">").str(), "</DataItem>", body.str());
	//	attribute->addChild(dataItem);
	//	primalEdgeGrid->addChild(attribute);


	//	// Add grids to gridList
	//	gridList->addChild(primalGrid);
	//	gridList->addChild(dualGrid);
	//	gridList->addChild(primalEdgeGrid);
	//	timeGrid->addChild(gridList);
	//}


	//file << root << std::endl;
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
	const int nDualFaces = dualMesh.getNumberOfFaces();
	Eigen::Matrix3Xd B(3, nPrimalFaces);
	Eigen::Matrix3Xd D(3, nDualFaces);
	Eigen::Matrix3Xd J(3, nDualFaces);

	for (int i = 0; i < nPrimalFaces; i++) {
		const Eigen::Vector3d fn = primalMesh.getFace(i)->getNormal();
		B.col(i) = bvec(i) * fn;
	}
	for (int i = 0; i < nDualFaces; i++) {
		const Eigen::Vector3d fn = dualMesh.getFace(i)->getNormal();
		D.col(i) = dvec(i) * fn;
		J.col(i) = jvec(i) * fn;
	}
	h5writer.writeData(B, "/B");
	h5writer.writeData(D, "/D");
	h5writer.writeData(J, "/J");

	const int nPrimalEdges = primalMesh.getNumberOfEdges();
	Eigen::Matrix3Xd E(3, nPrimalEdges);
	E.setZero();
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

