#include "DualMesh.h"



DualMesh::DualMesh()
	: Mesh("dual")
{
}

DualMesh::DualMesh(const std::string & meshPrefix) 
	: Mesh(meshPrefix)
{
}

DualMesh::~DualMesh()
{
}

void DualMesh::init_dualMesh(const PrimalMesh & primal)
{
	// add dual vertices at primal cell centers
	const std::vector<Cell*> primalCells = primal.getCells();
	for (auto cell : primalCells) {
		addVertex(cell->getCenter());
	}

	// add dual edges across primal faces
	const std::vector<Face*> primalFaces = primal.getFaces();
	for (auto face : primalFaces) {
		const std::vector<Cell*> adjCells = face->getCellList();
		if (face->isBoundary()) {
			assert(adjCells.size() == 1);
		}
		else {
			assert(adjCells.size() == 2);
			int idxA = -1; 
			int idxB = -1;
			for (auto cell : adjCells) {
				int orientation = cell->getOrientation(face);
				if (orientation > 0) {
					idxA = cell->getIndex();
				}
				else {
					idxB = cell->getIndex();
				}				
			}
			assert(idxA >= 0);
			assert(idxB >= 0);
			Vertex * A = getVertex(idxA);
			Vertex * B = getVertex(idxB);
			addEdge(A, B);
		}
	}

	// TODO set dual faces
}
