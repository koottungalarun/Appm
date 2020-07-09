#pragma once
#include <ostream>
#include "InterpolationTable.h"

class FrictionModel {

public:
	enum class Type
	{
		OFF, ON, TEST
	};

	friend std::ostream & operator<<(std::ostream & os, const Type & obj);

	FrictionModel();
	~FrictionModel();

	const Eigen::VectorXd getAverageMomentumCrossSection(const Eigen::VectorXd & x);


private:
	InterpolationTable table;

};