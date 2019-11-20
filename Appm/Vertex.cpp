#include "Vertex.h"



Vertex::Vertex()
{
}

Vertex::Vertex(const int index)
	: GeometryItem(index)
{
}

Vertex::Vertex(const Eigen::Vector3d & position)
{
	this->position = position;
}

Vertex::Vertex(const Eigen::Vector3d & position, const int index)
	: Vertex(position)
{
	setIndex(index);
}


Vertex::~Vertex()
{
}

const Eigen::Vector3d Vertex::getPosition() const
{
	return this->position;
}

void Vertex::setAdjacient(Edge * edge) 
{
	assert(edge != nullptr);
	adjacientEdges.push_back(edge);
}

bool Vertex::isAdjacientTo(const Vertex * other) const
{
	assert(other != nullptr);
	for (auto edge : adjacientEdges) {
		if (edge->isAdjacient(this, other)) {
			return true;
		}
	}
	return false;
}

Edge * Vertex::getAdjacientEdge(const Vertex * other)
{
	assert(other != nullptr);
	for (auto edge : adjacientEdges) {
		if (edge->isAdjacient(this, other)) {
			return edge;
		}
	}
	return nullptr;
}

std::vector<Edge*> Vertex::getEdges()
{
	return adjacientEdges;
}

std::ostream & operator<<(std::ostream & os, const Vertex & obj)
{
	os << "Vertex idx: " << obj.getIndex() << "; ";
	os << "edges = {";
	for (auto e : obj.adjacientEdges) {
		os << e->getIndex() << ",";
	}
	os << "}";
	return os;
}
