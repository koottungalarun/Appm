#include "Face.h"



Face::Face()
{
}

Face::Face(std::vector<Edge*> faceEdges)
{
	this->edgeList = faceEdges;
	for (auto edge : edgeList) {
		edge->setAdjacient(this);
	}

	if (edgeList.size() == 3) {
		center = getCircumCenter();
	}
	else {
		// arithmetic mean of face vertices
		center = Eigen::Vector3d::Zero();
		std::vector<Vertex*> vertexList = getVertexList();
		for (auto v : vertexList) {
			center += v->getPosition();
		}
		center /= vertexList.size();
	}
}


Face::~Face()
{
}

std::vector<Edge*> Face::getEdgeList() const
{
	return this->edgeList;
}

std::vector<Vertex*> Face::getVertexList() const
{
	std::vector<Vertex*> faceVertices;

	Edge * e = edgeList[0];
	Vertex * v = e->getVertexB();

	const int nEdges = edgeList.size();
	for (int i = 0; i < nEdges; i++) {
		faceVertices.push_back(v);
		Vertex * otherVertex = e->getOppositeVertex(v);
		assert(otherVertex != nullptr);
		Edge * otherEdge = nullptr;
		for (auto edge : edgeList) {
			if (edge != e && edge->hasVertex(otherVertex)) {
				otherEdge = edge;
			}
		}
		assert(otherEdge != nullptr);
		e = otherEdge;
		v = otherVertex;
	}

	assert(faceVertices.size() == edgeList.size());
	return faceVertices;
}

bool Face::hasFaceEdges(const std::vector<Edge*> faceEdges) const
{
	std::vector<int>::iterator it;
	std::vector<int> edgeIdx;
	for (auto edge : faceEdges) {
		edgeIdx.push_back(edge->getIndex());
	}
	std::sort(edgeIdx.begin(), edgeIdx.end());
	// remove duplicate edge indices 
	it = std::unique(edgeIdx.begin(), edgeIdx.end());
	edgeIdx.resize( std::distance(edgeIdx.begin(), it) );

	if (faceEdges.size() != edgeList.size()) {
		return false;
	}

	std::vector<int> thisFaceEdgeIdx;
	for (auto edge : edgeList) {
		thisFaceEdgeIdx.push_back(edge->getIndex());
	}
	std::sort(thisFaceEdgeIdx.begin(), thisFaceEdgeIdx.end());
	
	std::vector<int> result(edgeIdx.size());
	it = std::set_intersection(edgeIdx.begin(), edgeIdx.end(), thisFaceEdgeIdx.begin(), thisFaceEdgeIdx.end(), result.begin());
	const int distance = std::distance(result.begin(), it);
	if (distance < edgeList.size()) {
		return false;
	}

	return true;
}

const int Face::getOrientation(const Edge * edge)
{
	return 0;
}

bool Face::isBoundary() const
{
	for (auto edge : edgeList) {
		if (edge->isBoundary()) {
			return true;
		}
	}
	return false;
}

const Eigen::Vector3d Face::getCenter() const
{
	return center;
}

const Eigen::Vector3d Face::getCircumCenter() const
{
	assert(edgeList.size() == 3);
	const std::vector<Vertex*> vertexList = getVertexList();

	// Definition of circumcenter: http://mathworld.wolfram.com/Circumcircle.html
	double z = 0;
	Eigen::MatrixXd D(3, 4);
	for (int i = 0; i < 3; i++) {
		const Eigen::Vector3d pos = vertexList[i]->getPosition();
		D(i, 0) = pow(pos(0), 2) + pow(pos(1), 2);
		D(i, 1) = pos(0);
		D(i, 2) = pos(1);
		D(i, 3) = 1.0;
		z += pos(2);
	}
	z /= 3;

	Eigen::Matrix3d Bx, By;
	Bx.col(0) = D.col(0);
	Bx.col(1) = D.col(2);
	Bx.col(2) = D.col(3);

	By.col(0) = D.col(0);
	By.col(1) = D.col(1);
	By.col(2) = D.col(3);
	
	const double a = D.rightCols(3).determinant();
	assert(abs(a) > 0);
	const double bx = Bx.determinant();
	const double by = By.determinant();
	return Eigen::Vector3d(bx / (2*a), by / (2*a), z);
}

std::ostream & operator<<(std::ostream & os, const Face & obj)
{
	os << "Face " << obj.getIndex() << ": {";
	for (auto edge : obj.edgeList) {
		os << edge->getIndex() << ",";
	}
	os << "}";
	return os;
}
