#include "Mesh.h"



Mesh::Mesh()
{
}

Mesh::Mesh(const std::string & meshPrefix)
{
	this->meshPrefix = meshPrefix;
}


Mesh::~Mesh()
{
}

void Mesh::writeToFile()
{
	std::cout << std::endl;
	std::cout << "Mesh:     " << this->meshPrefix << std::endl;
	std::cout << "Vertices: " << vertexList.size() << std::endl;
	std::cout << "Edges:    " << edgeList.size() << std::endl;
	std::cout << "Faces:    " << faceList.size() << std::endl;
	std::ofstream file;

	const int nVertices = vertexList.size();
	const int nEdges = edgeList.size();

	// Write vertices to file
	Eigen::Matrix3Xd positions(3, nVertices);
	for (int i = 0; i < nVertices; i++) {
		positions.col(i) = vertexList[i]->getPosition();
	}
	file = std::ofstream(this->meshPrefix + "-vertices.dat");
	file << positions << std::endl;

	// Create incidence maps and write them to file
	createIncidenceMaps();
	Eigen::sparseMatrixToFile(edge2vertexMap, this->meshPrefix + "-e2v.dat");
	Eigen::sparseMatrixToFile(face2edgeMap, this->meshPrefix + "-f2e.dat");
	Eigen::sparseMatrixToFile(cell2faceMap, this->meshPrefix + "-c2f.dat");



	//for (auto v : vertexList) {
	//	std::cout << v->getIndex() << ": " << v->getPosition().transpose() << std::endl;
	//}
	//for (auto e : edgeList) {
	//	std::cout << *e << std::endl;
	//}
	//for (auto f : faceList) {
	//	std::cout << *f << std::endl;
	//}
}

void Mesh::writeXdmf()
{
	writeXdmf_surface();
	writeXdmf_volume();
}

Vertex * Mesh::addVertex(const Eigen::Vector3d & position)
{
	const int index = vertexList.size();
	Vertex * vertex = new Vertex(position, index);
	vertexList.push_back(vertex);
	return vertex;
}

Edge * Mesh::addEdge(Vertex * A, Vertex * B)
{
	assert(A != nullptr);
	assert(B != nullptr);
	if (!isConnected(A, B)) {
		Edge * edge = new Edge(A, B);
		edge->setIndex(edgeList.size());
		edgeList.push_back(edge);
		return edge;
	}
	else {
		return getEdge(A, B);
	}
}

Face * Mesh::addFace(const std::vector<Edge*> & faceEdges)
{
	assert(faceEdges.size() >= 3);
	if (!isConnected(faceEdges)) {
		Face * face = new Face(faceEdges);
		face->setIndex(faceList.size());
		faceList.push_back(face);
		return face;
	}
	else {
		Face * face = getFace(faceEdges);
	}
	return nullptr;
}

Vertex * Mesh::getVertex(const int index)
{
	assert(index >= 0);
	assert(index < vertexList.size());
	Vertex * vertex = vertexList[index];
	assert(vertex != nullptr);
	assert(vertex->getIndex() == index);
	return vertex;
}

Face * Mesh::getFace(const std::vector<Edge*>& faceEdges)
{
	assert(this->isConnected(faceEdges));
	assert(faceEdges[0]->hasConnectedFace(faceEdges));
	return faceEdges[0]->getConnectedFace(faceEdges);
}

void Mesh::createIncidenceMaps()
{
	create_vertexCoordinates();
	create_edge2vertex_map();
	create_face2edge_map();
	create_cell2face_map();
}

bool Mesh::isConnected(const Vertex * A, const Vertex * B) const
{
	return (A->isAdjacientTo(B) && B->isAdjacientTo(A));
}

bool Mesh::isConnected(const std::vector<Edge*> & faceEdges) 
{
	bool isAlreadyDefinedFace = true;
	for (auto edge : faceEdges) {
		isAlreadyDefinedFace &= edge->hasConnectedFace(faceEdges);
		if (!isAlreadyDefinedFace) {
			break;
		}
	}
	return isAlreadyDefinedFace;
}

Edge * Mesh::getEdge(Vertex * A, Vertex * B) const
{
	assert(A != nullptr);
	assert(B != nullptr);
	if (isConnected(A, B)) {
		return A->getAdjacientEdge(B);
	}
	else {
		return nullptr;
	}
}

void Mesh::create_vertexCoordinates()
{
	const int nVertices = vertexList.size();
	vertexCoordinates = Eigen::Matrix3Xd(3, nVertices);
	for (int i = 0; i < nVertices; i++) {
		vertexCoordinates.col(i) = vertexList[i]->getPosition();
	}
}

void Mesh::create_edge2vertex_map()
{
	const int nEdges = edgeList.size();
	const int nVertices = vertexList.size();
	this->edge2vertexMap = Eigen::SparseMatrix<int>(nEdges, nVertices);
	typedef Eigen::Triplet<int> T;
	std::vector<T> triplets;
	for (int i = 0; i < nEdges; i++) {
		auto e = edgeList[i];
		int idxE = e->getIndex();
		const Vertex * A = e->getVertexA();
		const Vertex * B = e->getVertexB();
		const int idxA = A->getIndex();
		const int idxB = B->getIndex();
		triplets.push_back(T(idxE, idxA, -1));
		triplets.push_back(T(idxE, idxB, +1));
	}
	edge2vertexMap.setFromTriplets(triplets.begin(), triplets.end());
	edge2vertexMap.makeCompressed();
}

void Mesh::create_face2edge_map()
{
	const int nFaces = faceList.size();
	const int nEdges = edgeList.size();
	this->face2edgeMap = Eigen::SparseMatrix<int>(nFaces, nEdges);
	typedef Eigen::Triplet<int> T;
	std::vector<T> triplets;
	for (int i = 0; i < nFaces; i++) {
		auto face = faceList[i];
		const int idxF = face->getIndex();
		const std::vector<Edge*> faceEdges = face->getEdgeList();
		for (auto edge : faceEdges) {
			auto idxE = edge->getIndex();
			const int orientation = face->getOrientation(edge);
			triplets.push_back(T(idxF, idxE, orientation));
		}
	}
	face2edgeMap.setFromTriplets(triplets.begin(), triplets.end());
	face2edgeMap.makeCompressed();
}

void Mesh::create_cell2face_map()
{
}

void Mesh::writeXdmf_surface()
{
	std::stringstream ss;
	std::stringstream body;
	XmlElement * dataItem;

	std::string filename = "appm-" + this->meshPrefix + "-surface.xdmf";
	std::cout << "Write XDMF file: " << filename << std::endl;

	XmlElement root("<Xdmf Version=\"3.0\" xmlns:xi=\"[http://www.w3.org/2001/XInclude]\">", "</Xdmf>");
	XmlElement * domain = new XmlElement("<Domain>", "</Domain>");
	root.addChild(domain);

	XmlElement * grid = new XmlElement("<Grid Name=\"Appm faces\">", "</Grid>");
	domain->addChild(grid);

	ss << "<Topology TopologyType=\"Mixed\">";
	XmlElement * topology = new XmlElement(ss.str(), "</Topology>");
	grid->addChild(topology);

	const int nFaces = faceList.size();
	Eigen::MatrixXd face2vertex(nFaces, 4);
	for (int i = 0; i < nFaces; i++) {
		Eigen::VectorXd temp(4);
		temp(0) = 4;
		Face * face = faceList[i];
		std::vector<Vertex*> faceVertices = face->getVertexList();
		const int nFaceVertices = faceVertices.size();
		assert(nFaceVertices == 3);
		for (int j = 0; j < faceVertices.size(); j++) {
			temp(j + 1) = faceVertices[j]->getIndex();
		}
		face2vertex.row(i) = temp;
	}
	body << face2vertex;
	
	ss = std::stringstream();
	ss << "<DataItem Dimensions=\"" << face2vertex.array().size() << "\" DataType=\"Int\" Format=\"XML\">";

	dataItem = new XmlElement(ss.str(), "</DataItem>", body.str());
	topology->addChild(dataItem);

	XmlElement * geometry = new XmlElement("<Geometry GeometryType=\"XYZ\">", "</Geometry>");
	grid->addChild(geometry);
	ss = std::stringstream();
	ss << "<DataItem Dimensions=\"" << vertexCoordinates.cols() << " " << 3 << "\" DataType=\"Float\" Precision=\"8\" Format=\"XML\">";
	body = std::stringstream();
	body << vertexCoordinates.transpose();
	dataItem = new XmlElement(ss.str(), "</DataItem>", body.str());
	geometry->addChild(dataItem);


	std::ofstream file(filename);
	file << "<?xml version=\"1.0\" ?>" << std::endl;
	file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
	file << root << std::endl;
}

void Mesh::writeXdmf_volume()
{
}
