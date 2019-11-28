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

	const int nDualFaces = dualMesh.getNumberOfFaces();
	for (int i = 0; i < nDualFaces; i++) {
		const Face * face = dualMesh.getFace(i);

		const Eigen::Vector3d fc = face->getCenter();
		const Eigen::Vector2d fc_2d = fc.segment(0, 2);

		if (fc_2d.norm() < 0.3) {
			const Eigen::Vector3d fn = face->getNormal();
			jvec(i) = 1 * fn.dot(Eigen::Vector3d(0, 0, 1));
		}
	}
	double time = 0;
	int iteration = 0;
	writeOutput(iteration, time);


	double timestepSize = 1;
	const int maxIteration = 1;
	const double maxTime = 10;
	while (iteration < maxIteration && time < maxTime) {
		iteration++;
		time += timestepSize;
		writeOutput(iteration, time);
	}

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

	XmlElement root("<Xdmf Version=\"3.0\" xmlns:xi=\"[http://www.w3.org/2001/XInclude]\">", "</Xdmf>");
	XmlElement * domain = new XmlElement("<Domain>", "</Domain>");
	root.addChild(domain);

	XmlElement * timeGrid = new XmlElement("<Grid Name=\"Time Grid\" GridType=\"Collection\" CollectionType=\"Temporal\">", "</Grid>");
	domain->addChild(timeGrid);

	std::stringstream body;
	std::stringstream headTag;

	XmlElement * topology;
	XmlElement * dataItem;
	XmlElement * geometry;
	XmlElement * attribute;

	H5Reader h5reader;
	h5reader = H5Reader("primal-mesh.h5");
	const int nElemPrimalTopo = h5reader.readDataSize("/face2vertex");
	assert(nElemPrimalTopo > 0);
	
	h5reader = H5Reader("dual-mesh.h5");
	const int nElemDualTopo = h5reader.readDataSize("/face2vertex");
	assert(nElemDualTopo > 0);



	for (int i = 0; i < nTimesteps; i++) {
		
		const std::string datafilename = (std::stringstream() << "appm-" << i << ".h5").str();
		std::cout << "loop: " << i << ", " << datafilename << std::endl;

		XmlElement * twoGrids = new XmlElement("<Grid Name=\"Two Grids\" GridType=\"Tree\">", "</Grid>");
		
		XmlElement * timeTag = new XmlElement((std::stringstream() << "<Time Value=\"" << timeStamps[i] << "\"/>").str());
		twoGrids->addChild(timeTag);
	
		XmlElement * primalGrid = new XmlElement("<Grid Name=\"Primal\" GridType=\"Uniform\">", "</Grid>");
		topology = new XmlElement((std::stringstream() << "<Topology TopologyType=\"Mixed\" NumberOfElements=\"" << primalMesh.getNumberOfFaces() << "\">").str(), "</Topology>");
		body = std::stringstream();
		body << primalMesh.getPrefix() << "-mesh.h5:/face2vertex";
		headTag = std::stringstream();
		headTag << "<DataItem Dimensions=\"" << nElemPrimalTopo << "\" DataType=\"Int\" Format=\"HDF\">";
		dataItem = new XmlElement(headTag.str(), "</DataItem>", body.str());
		topology->addChild(dataItem);
		primalGrid->addChild(topology);

		geometry = new XmlElement("<Geometry GeometryType=\"XYZ\">", "</Geometry>");
		headTag = std::stringstream();
		headTag << "<DataItem Dimensions=\"" << primalMesh.getNumberOfVertices() << " 3\" DataType=\"Float\" Precision=\"8\" Format=\"HDF\">";
		body = std::stringstream() << primalMesh.getPrefix() << "-mesh.h5:/vertexPos";
		dataItem = new XmlElement(headTag.str(), "</DataItem>", body.str());
		geometry->addChild(dataItem);
		primalGrid->addChild(geometry);

		// Attribute: B-field
		attribute = new XmlElement("<Attribute Name=\"Magnetic Flux\" AttributeType=\"Vector\" Center=\"Cell\">", "</Attribute>");
		body = std::stringstream() << datafilename << ":" << "/B";
		headTag = std::stringstream() << "<DataItem Dimensions=\"" << primalMesh.getNumberOfFaces() << " 3\" DataType=\"Float\" Precision=\"8\" Format=\"HDF\">";
		dataItem = new XmlElement(headTag.str(), "</DataItem>", body.str());
		attribute->addChild(dataItem);
		primalGrid->addChild(attribute);

		// Dual grid
		XmlElement * dualGrid = new XmlElement("<Grid Name=\"Dual\" GridType=\"Uniform\">", "</Grid>");
		topology = new XmlElement((std::stringstream() << "<Topology TopologyType=\"Mixed\" NumberOfElements=\"" << dualMesh.getNumberOfFaces() << "\">").str(), "</Topology>");
		body = std::stringstream();
		body << dualMesh.getPrefix() << "-mesh.h5:/face2vertex";
		headTag = std::stringstream();
		headTag << "<DataItem Dimensions=\"" << nElemDualTopo << "\" DataType=\"Int\" Format=\"HDF\">";
		dataItem = new XmlElement(headTag.str(), "</DataItem>", body.str());
		topology->addChild(dataItem);
		dualGrid->addChild(topology);

		geometry = new XmlElement("<Geometry GeometryType=\"XYZ\">", "</Geometry>");
		headTag = std::stringstream();
		headTag << "<DataItem Dimensions=\"" << dualMesh.getNumberOfVertices() << " 3\" DataType=\"Float\" Precision=\"8\" Format=\"HDF\">";
		body = std::stringstream() << dualMesh.getPrefix() << "-mesh.h5:/vertexPos";
		dataItem = new XmlElement(headTag.str(), "</DataItem>", body.str());
		geometry->addChild(dataItem);
		dualGrid->addChild(geometry);

		// Attribute: D-field
		attribute = new XmlElement("<Attribute Name=\"Displacement Field\" AttributeType=\"Vector\" Center=\"Cell\">", "</Attribute>");
		body = std::stringstream() << datafilename << ":" << "/D";
		headTag = std::stringstream() << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << " 3\" DataType=\"Float\" Precision=\"8\" Format=\"HDF\">";
		dataItem = new XmlElement(headTag.str(), "</DataItem>", body.str());
		attribute->addChild(dataItem);
		dualGrid->addChild(attribute);

		// Attribute: J-field
		attribute = new XmlElement("<Attribute Name=\"Current\" AttributeType=\"Vector\" Center=\"Cell\">", "</Attribute>");
		body = std::stringstream() << datafilename << ":" << "/J";
		headTag = std::stringstream() << "<DataItem Dimensions=\"" << dualMesh.getNumberOfFaces() << " 3\" DataType=\"Float\" Precision=\"8\" Format=\"HDF\">";
		dataItem = new XmlElement(headTag.str(), "</DataItem>", body.str());
		attribute->addChild(dataItem);
		dualGrid->addChild(attribute);

		twoGrids->addChild(primalGrid);
		twoGrids->addChild(dualGrid);
		timeGrid->addChild(twoGrids);
	}


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


	Eigen::VectorXd timeVec(1);
	timeVec(0) = time;
	h5writer.writeData(timeVec, "/time");
	Eigen::VectorXi iterVec(1);
	iterVec(0) = iteration;
	h5writer.writeData(iterVec, "/iteration");

}

