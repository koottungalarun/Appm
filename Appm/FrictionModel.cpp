#include "FrictionModel.h"

FrictionModel::FrictionModel()
{
}

FrictionModel::~FrictionModel()
{
}

const Eigen::VectorXd FrictionModel::getAverageMomentumCrossSection(const Eigen::VectorXd & x)
{
	return table.interpolate(x);
}

std::ostream & operator<<(std::ostream & os, const FrictionModel::Type & obj) {
	switch (obj) {
	case FrictionModel::Type::OFF:
		os << "OFF";
		break;
	case FrictionModel::Type::ON:
		os << "ON";
		break;
	case FrictionModel::Type::TEST:
		os << "TEST";
		break;
	}

	return os;
}


