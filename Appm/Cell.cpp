#include "Cell.h"



Cell::Cell()
{
}

Cell::Cell(const std::vector<Face*>& faces)
{
	//std::cout << "Create cell from faces: {";
	//for (auto face : faces) {
	//	std::cout << face->getIndex() << ",";
	//}
	//std::cout << "}" << std::endl;
	//std::cout << "Face vertices: " << std::endl;
	//for (auto face : faces) {
	//	std::cout << "face idx " << face->getIndex() << ": ";
	//	const std::vector<Vertex*> f2v = face->getVertexList();
	//	for (auto v : f2v) {
	//		std::cout << v->getIndex() << ",";
	//	}
	//	std::cout << std::endl;
	//}

	this->faceList = faces;
	for (auto face : faceList) {
		face->setAdjacient(this);
	}

	// determine cell center
	std::vector<Face*> zFaces;
	for (auto face : faceList) {
		const Eigen::Vector3d fn = face->getNormal();
		if (fn.cross(Eigen::Vector3d(0, 0, 1)).norm() < 4 * std::numeric_limits<double>::epsilon()) {
			zFaces.push_back(face);
		}
	}
	assert(zFaces.size() == 2);
	//std::cout << "z-Faces: ";
	//for (auto f : zFaces) {
	//	std::cout << f->getIndex() << ",";
	//}
	//std::cout << std::endl;
	this->center = Eigen::Vector3d::Zero();
	center = 1. / 2. * (zFaces[0]->getCenter() + zFaces[1]->getCenter());

	//Eigen::Matrix3d M;
	//M.col(0) = zFaces[0]->getCenter();
	//M.col(1) = zFaces[1]->getCenter();
	//M.col(2) = this->center;
	//std::cout << "face/cell center: " << std::endl << M << std::endl;
}


Cell::~Cell()
{
}

const double Cell::getVolume() const
{
	return volume;
}

const std::vector<Face*> Cell::getFaceList() const
{
	return faceList;
}

bool Cell::hasFace(const Face * face) const
{
	assert(face != nullptr);
	for (auto f : faceList) {
		if (f == face) {
			return true;
		}
	}
	return false;
}

const int Cell::getOrientation(const Face * face) const
{
	assert(face != nullptr);
	bool isMember = false;
	for (auto f : faceList) {
		if (f == face) {
			isMember |= true;
			break;
		}
	}
	assert(isMember);
	const Eigen::Vector3d a = face->getCenter() - this->center;
	const Eigen::Vector3d fn = face->getNormal();
	const int orientation = (a.dot(fn) > 0) ? 1 : -1;
	return orientation;
}

const Eigen::Vector3d & Cell::getCenter() const
{
	return center;
}
