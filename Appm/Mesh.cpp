#include "Mesh.h"



Mesh::Mesh()
{
	//std::cout << "Call to Mesh()" << std::endl;
}

Mesh::Mesh(const std::string & meshPrefix)
{
	//std::cout << "Call to Mesh(string)" << std::endl;
	this->meshPrefix = meshPrefix;
}


Mesh::~Mesh()
{
	//std::cout << "Call to ~Mesh()" << std::endl;
	if (vertexList.size() > 0) {
		for (auto v : vertexList) {
			delete v;
		}
		vertexList.resize(0);
	}
	if (edgeList.size() > 0) {
		for (auto e : edgeList) {
			delete e;
		}
		edgeList.resize(0);
	}
	if (faceList.size() > 0) {
		for (auto f : faceList) {
			delete f;
		}
		faceList.resize(0);
	}
	//if (cellList.size() > 0) {
	//	for (auto c : cellList) {
	//		delete c;
	//	}
	//	cellList.resize(0);
	//}
}

void Mesh::writeToFile()
{
	std::cout << std::endl;
	std::cout << "Mesh:     " << this->meshPrefix << std::endl;
	std::cout << "Vertices: " << vertexList.size() << std::endl;
	std::cout << "Edges:    " << edgeList.size() << std::endl;
	std::cout << "Faces:    " << faceList.size() << std::endl;
	std::cout << "Cells:    " << cellList.size() << std::endl;

	std::ofstream file;

	const int nVertices = vertexList.size();
	const int nEdges = edgeList.size();

	// Write vertices to file
	Eigen::Matrix3Xd positions(3, nVertices);
	for (int i = 0; i < nVertices; i++) {
		positions.col(i) = vertexList[i]->getPosition();
	}
	file = std::ofstream(this->meshPrefix + "-vertices.dat");
	file << positions.transpose() << std::endl;

	file = std::ofstream(this->meshPrefix + "-coords.dat");
	file << vertexCoordinates.transpose() << std::endl;

	// Create incidence maps and write them to file
	createIncidenceMaps();
	Eigen::sparseMatrixToFile(edge2vertexMap, this->meshPrefix + "-e2v.dat");
	Eigen::sparseMatrixToFile(face2edgeMap, this->meshPrefix + "-f2e.dat");
	Eigen::sparseMatrixToFile(cell2faceMap, this->meshPrefix + "-c2f.dat");

	const int nFaces = faceList.size();
	Eigen::MatrixXi f2v(3, nFaces);
	for (int j = 0; j < nFaces; j++) {
		const Face * face = faceList[j];
		std::vector<Vertex*> faceVertices = face->getVertexList();
		const int nFaceVertices = faceVertices.size();
		if (nFaceVertices > f2v.rows()) {
			f2v.conservativeResize(nFaceVertices, f2v.cols());
		}
		for (int k = 0; k < nFaceVertices; k++) {
			f2v(k, j) = faceVertices[k]->getIndex();
		}
	}
	file = std::ofstream(this->meshPrefix + "-f2v.dat");
	file << f2v.transpose() << std::endl;

	file = std::ofstream(this->meshPrefix + "-faceCenter.dat");
	Eigen::Matrix3Xd fc(3, nFaces);
	for (int i = 0; i < nFaces; i++) {
		fc.col(i) = getFace(i)->getCenter();
	}
	file << fc.transpose() << std::endl;

	file = std::ofstream(this->meshPrefix + "-faceNormal.dat");
	Eigen::Matrix3Xd fn(3, nFaces);
	for (int i = 0; i < nFaces; i++) {
		fn.col(i) = getFace(i)->getNormal();
	}
	file << fn.transpose() << std::endl;

	file = std::ofstream(this->meshPrefix + "-faceBoundary.dat");
	Eigen::VectorXi faceBoundary(nFaces);
	for (int i = 0; i < nFaces; i++) {
		faceBoundary(i) = getFace(i)->isBoundary();
	}
	file << faceBoundary << std::endl;

	file = std::ofstream(this->meshPrefix + "-cellCenter.dat");
	const int nCells = cellList.size();
	Eigen::Matrix3Xd cellCenters(3, nCells);
	for (int i = 0; i < nCells; i++) {
		cellCenters.col(i) = getCell(i)->getCenter();
	}
	file << cellCenters.transpose() << std::endl;

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
		return face;
	}
	return nullptr;
}
//Face * Mesh::addFace(const std::vector<Vertex*> & faceVertices)
//{
//	assert(faceVertices.size() >= 3);
//	std::vector<Edge*> faceEdges;
//	const int nVertices = faceVertices.size();
//	for (int i = 0; i < nVertices; i++) {
//		Vertex * A = faceVertices[i];
//		Vertex * B = faceVertices[(i + 1) % nVertices];
//		Edge * edge = getEdge(A, B);
//		assert(edge != nullptr);
//		faceEdges.push_back(edge);
//	}
//	if (!isConnected(faceEdges)) {
//		Face * face = new Face(faceVertices);
//		face->setIndex(faceList.size());
//		faceList.push_back(face);
//		return face;
//	}
//	else {
//		Face * face = getFace(faceEdges);
//		return face;
//	}
//	return nullptr;
//}



TriPrism * Mesh::addTriPrism(const std::vector<Face*>& sideFaces, Face * bottomFace, Face * topFace)
{
	for (auto f : sideFaces) {
		assert(f != nullptr);
	}
	assert(bottomFace != nullptr);
	assert(topFace != nullptr);

	std::vector<Face*> cellFaces;
	for (auto f : sideFaces) {
		cellFaces.push_back(f);
	}
	cellFaces.push_back(bottomFace);
	cellFaces.push_back(topFace);

	if (!isConnected(cellFaces)) {
		const int idx = cellList.size();
		TriPrism * prism = new TriPrism(cellFaces);
		prism->setIndex(idx);
		cellList.push_back(prism);
		return prism;
	}
	else {
		Cell * cell = getCell(cellFaces);
		assert(cell != nullptr);
		return static_cast<TriPrism*>(cell);
	}
	return nullptr;
}

Cell * Mesh::addCell(const std::vector<Face*>& cellFaces)
{
	assert(cellFaces.size() >= 4);
	if (!isConnected(cellFaces)) {
		const int idx = cellList.size();
		Cell * cell = new Cell(cellFaces);
		cell->setIndex(idx);
		cellList.push_back(cell);
		return cell;
	}
	else {
		Cell * cell = getCell(cellFaces);
		assert(cell != nullptr);
		return cell;
	}
	return nullptr;
}

Cell * Mesh::getCell(const std::vector<Face*>& cellFaces)
{
	assert(cellFaces.size() >= 2);
	const Face * face0 = cellFaces[0];
	const Face * face1 = cellFaces[1];
	assert(face0 != face1);
	assert(face0->hasCommonCell(face1));
	Cell * cell = face0->getCommonCell(face1);
	for (auto face : cellFaces) {
		assert(cell->hasFace(face));
	}
	return cell;
}

Vertex * Mesh::getVertex(const int index) const
{
	assert(index >= 0);
	assert(index < vertexList.size());
	Vertex * vertex = vertexList[index];
	assert(vertex != nullptr);
	assert(vertex->getIndex() == index);
	return vertex;
}

Edge * Mesh::getEdge(const int index) const
{
	assert(index >= 0);
	assert(index < edgeList.size());
	return edgeList[index];
}

Face * Mesh::getFace(const int index) const
{
	assert(index >= 0);
	assert(index < faceList.size());
	return faceList[index];
}

Face * Mesh::getFace(const std::vector<Edge*>& faceEdges)
{
	assert(this->isConnected(faceEdges));
	assert(faceEdges[0]->hasConnectedFace(faceEdges));
	return faceEdges[0]->getConnectedFace(faceEdges);
}

Cell * Mesh::getCell(const int index) const
{
	assert(index >= 0);
	assert(index < cellList.size());
	return cellList[index];
}

void Mesh::createIncidenceMaps()
{
	create_vertexCoordinates();
	create_edge2vertex_map();
	create_face2edge_map();
	create_cell2face_map();
}

const int Mesh::getNumberOfVertices() const
{
	return vertexList.size();
}

const int Mesh::getNumberOfEdges() const
{
	return edgeList.size();
}

const int Mesh::getNumberOfFaces() const
{
	return faceList.size();
}

const std::vector<Cell*> Mesh::getCells() const
{
	return cellList;
}

const std::vector<Face*> Mesh::getFaces() const
{
	return faceList;
}

void Mesh::check() const
{
	//for (auto cell : cellList) {
	//	std::cout << "Cell idx " << cell->getIndex() << std::endl;
	//	const std::vector<Face*> cellFaces = cell->getFaceList();
	//	for (auto face : cellFaces) {
	//		std::cout << "Face " << face->getIndex() << "; ";
	//		std::cout << "orientation = " << cell->getOrientation(face);
	//		std::cout << std::endl;
	//	}
	//}

	if (cellList.size() > 0) {
		// face is either a boundary face, or it has two adjacient cells; orientation +/-1.
		for (auto face : faceList) {
			std::vector<Cell*> faceCells = face->getCellList();
			if (face->isBoundary()) {
				assert(faceCells.size() == 1);
			}
			else {
				assert(faceCells.size() == 2);
				const Cell * cell0 = faceCells[0];
				const Cell * cell1 = faceCells[1];
				std::vector<int> orientations = { cell0->getOrientation(face), cell1->getOrientation(face) };
				int sumOrientations = 0;
				for (auto value : orientations) {
					sumOrientations += value;
				}
				if (sumOrientations != 0) {
					std::cout << "Face idx " << face->getIndex() << " has wrong set of orientation: "  << orientations[0] << ", " << orientations[1] << std::endl;
					std::cout << "face normal:  " << face->getNormal().transpose() << std::endl;
					std::cout << "face center:  " << face->getCenter().transpose() << std::endl;
					std::cout << "cell0 center: " << cell0->getCenter().transpose() << std::endl;
					std::cout << "cell1 center: " << cell1->getCenter().transpose() << std::endl;
					std::cout << "(c0 - fc):    " << (cell0->getCenter() - face->getCenter()).transpose() << std::endl;
					std::cout << "(c1 - fc):    " << (cell1->getCenter() - face->getCenter()).transpose() << std::endl;
					std::cout << "(c0 - fc).fn:    " << (cell0->getCenter() - face->getCenter()).dot(face->getNormal()) << std::endl;
					std::cout << "(c1 - fc).fn:    " << (cell1->getCenter() - face->getCenter()).dot(face->getNormal()) << std::endl;

				}

				assert(orientations[0] + orientations[1] == 0);
			}
		}
	}

}

const Eigen::SparseMatrix<int>& Mesh::get_f2eMap() const
{
	return face2edgeMap;
}

const Eigen::SparseMatrix<int>& Mesh::get_e2vMap() const
{
	return edge2vertexMap;
}

std::vector<Edge*> Mesh::makeContinuousLoop(std::vector<Edge*> edges)
{
	std::vector<Edge*> loop;
	std::vector<Edge*>::iterator it;
	Edge * edge = edges[0];
	Vertex * V = edge->getVertexB();
	const int nEdges = edges.size();

	for (int i = 0; i < nEdges; i++) {
		// add edge to new list
		loop.push_back(edge);
		// remove edge from old list
		for (it = edges.begin(); it != edges.end(); it++) {
			if (*it == edge) { break; }
		}
		edges.erase(it);
		// update vertex that is searched 
		V = edge->getOppositeVertex(V);
		// find edge that has this vertex
		for (it = edges.begin(); it != edges.end(); it++) {
			edge = *it;
			if (edge->hasVertex(V)) {
				break;
			}
		}
	}
	assert(edges.size() == 0);
	assert(loop.size() == nEdges);
	return loop;
}

bool Mesh::isConnected(const Vertex * A, const Vertex * B) const
{
	return (A->isAdjacientTo(B) && B->isAdjacientTo(A));
}

bool Mesh::isConnected(const std::vector<Edge*> & faceEdges) const
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

bool Mesh::isConnected(const std::vector<Face*> & cellFaces) const
{
	bool isAlreadyDefined = true;
	for (auto face : cellFaces) {
		isAlreadyDefined &= face->hasCommonCell(face);
	}
	return isAlreadyDefined;
}

Edge * Mesh::getEdge(Vertex * A, Vertex * B)
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
	const int nFaces = faceList.size();
	const int nCells = cellList.size();
	this->cell2faceMap = Eigen::SparseMatrix<int>(nCells, nFaces);
	typedef Eigen::Triplet<int> T;
	std::vector<T> triplets;
	for (int i = 0; i < nCells; i++) {
		auto cell = cellList[i];
		const int idxC = cell->getIndex();
		const std::vector<Face*> cellFaces = cell->getFaceList();
		for (auto face : cellFaces) {
			auto idxF = face->getIndex();
			const int orientation = cell->getOrientation(face);
			triplets.push_back(T(idxC, idxF, orientation));
		}
	}
	cell2faceMap.setFromTriplets(triplets.begin(), triplets.end());
	cell2faceMap.makeCompressed();

	// Check cell-to-face map
	// Each face has either one adjacient cell (boundary face), or 
	// it has two adjacient cells (non-boundary face = inner face). 
	// Moreover, because of positive/negative orientations, 
	// a colwise sum of the matrix elements yields +/-1 (boundary) or 0 (non-boundary).
	//Eigen::VectorXi colwiseSum = cell2faceMap.transpose() * Eigen::VectorXi::Ones(nCells);
	//Eigen::VectorXi isBoundaryFace(nFaces);
	//for (int i = 0; i < nFaces; i++) {
	//	isBoundaryFace(i) = getFace(i)->isBoundary();
	//}
	//assert(colwiseSum.size() == nFaces);
	//Eigen::VectorXi notBoundaryFace = isBoundaryFace.array() - 1;
	//assert((colwiseSum.cwiseProduct(notBoundaryFace).array() == 0).all());
	//Eigen::VectorXi temp = colwiseSum.cwiseAbs().cwiseProduct(isBoundaryFace);
	//assert((temp.array() >= 0).all());
	//assert((temp.array() <= 1).all());
	//assert(temp.sum() == isBoundaryFace.sum());
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
	std::vector<int> f2v;
	for (int i = 0; i < nFaces; i++) {
		const Face * face = faceList[i];
		std::vector<Vertex*> faceVertices = face->getVertexList();
		const int nFaceVertices = faceVertices.size();

		int faceType = 0;
		switch (nFaceVertices) {
		case 3: faceType = 4; break;
		case 4: faceType = 5; break;
		default: faceType = 3; break;
		}
		f2v.push_back(faceType);
		for (int j = 0; j < faceVertices.size(); j++) {
			f2v.push_back(faceVertices[j]->getIndex());
		}
	}
	for (int i = 0; i < f2v.size(); i++) {
		body << f2v[i] << " ";
		if ((i+1) % 10 == 0) {
			body << std::endl;
		}
	}
	ss = std::stringstream();
	ss << "<DataItem Dimensions=\"" << f2v.size() << "\" DataType=\"Int\" Format=\"XML\">";

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
	std::stringstream ss;
	std::stringstream body;
	XmlElement * dataItem;

	std::string filename = "appm-" + this->meshPrefix + "-volume.xdmf";
	std::cout << "Write XDMF file: " << filename << std::endl;

	XmlElement root("<Xdmf Version=\"3.0\" xmlns:xi=\"[http://www.w3.org/2001/XInclude]\">", "</Xdmf>");
	XmlElement * domain = new XmlElement("<Domain>", "</Domain>");
	root.addChild(domain);

	XmlElement * grid = new XmlElement("<Grid Name=\"Appm cells\">", "</Grid>");
	domain->addChild(grid);

	ss << "<Topology TopologyType=\"Mixed\">";
	XmlElement * topology = new XmlElement(ss.str(), "</Topology>");
	grid->addChild(topology);
	
	std::vector<int> data;
	const int nCells = cellList.size();
	for (int i = 0; i < nCells; i++) {
		const Cell * cell = cellList[i];
		data.push_back(16); // polygon type
		const std::vector<Face*> cellFaces = cell->getFaceList();
		const int nCellFaces = cellFaces.size();
		data.push_back(nCellFaces); // number of faces
		for (int j = 0; j < nCellFaces; j++) {
			const Face * face = cellFaces[j];
			const std::vector<Vertex*> faceVertices = face->getVertexList();
			const int nFaceVertices = faceVertices.size();
			data.push_back(nFaceVertices); // number of vertices
			for (auto vertex : faceVertices) {
				data.push_back(vertex->getIndex());
			}
		}
	}
	for (int i = 0; i < data.size(); i++) {
		body << data[i] << " ";
		if ((i + 1) % 10 == 0) {
			body << std::endl;
		}
	}
	ss = std::stringstream();
	ss << "<DataItem Dimensions=\"" << data.size() << "\" DataType=\"Int\" Format=\"XML\">";

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
