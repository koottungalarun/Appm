#pragma once
#include "GeometryItem.h"

#include "Face.h"
#include <limits>
class Face;

class Cell :
	public GeometryItem
{
public:

	enum class Type {
		DEFAULT, FLUID, SOLID, GHOST, ELECTRODE
	};

	Cell();
	Cell(const std::vector<Face*> & faces);
	~Cell();

	const double getVolume() const;
	const std::vector<Face*> getFaceList() const;
	bool hasFace(const Face * face) const;
	const int getOrientation(const Face * face) const;

	const Eigen::Vector3d & getCenter() const;
	const Eigen::Matrix3Xd getVertexCoordinates() const;

	void setType(const Type & type);
	const Type getType() const;

	bool validateCellGeometry() const;

private:
	double volume = 0;
	Eigen::Vector3d center;
	std::vector<Face*> faceList;

	Type fluidType = Type::DEFAULT;

	Eigen::Vector3d computeCellCenter() const;

};

