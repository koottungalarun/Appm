#include "PrimalMesh.h"



PrimalMesh::PrimalMesh() :
	PrimalMesh("primal")
{
	init();
}

PrimalMesh::PrimalMesh(const std::string & meshPrefix)
	: Mesh(meshPrefix)
{
}


PrimalMesh::~PrimalMesh()
{
}

void PrimalMesh::init()
{
	Vertex * origin = addVertex(Eigen::Vector3d(0, 0, 0));
	const int corners = 6;
	for (int k = 0; k < corners; k++) {
		const double phi = 2 * M_PI * k / corners;
		const Eigen::Vector3d pos(sin(phi), cos(phi), 0);
		Vertex * v = addVertex(pos);
		addEdge(origin, v);
	}

	for (int k = 0; k < corners; k++) {
		Vertex * A = getVertex(k + 1);
		Vertex * B = getVertex(((k + 1) % corners) + 1);
		auto edgeA = addEdge(origin, A);
		auto edgeB = addEdge(origin, B);
		auto edgeC = addEdge(A, B);
		std::vector<Edge*> faceEdges = { edgeA, edgeB, edgeC };
		Face * face = addFace(faceEdges);
	}

}
