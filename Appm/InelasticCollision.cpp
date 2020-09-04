#include "InelasticCollision.h"

InelasticCollision::InelasticCollision()
{
}

InelasticCollision::InelasticCollision(const std::string & folderPath)
{
	const DataTransform xTrans = DataTransform::INVERSE;
	const DataTransform yTrans = DataTransform::NONE;
	const DataTransform fTrans = DataTransform::LOG;

	std::string filename;
	filename = folderPath + "/I_Gion.csv";
	data_Gion = Interpolation1d(filename, xTrans, fTrans);

	filename = folderPath + "/I_Grec.csv";
	data_Grec = Interpolation2d(filename, xTrans, yTrans, fTrans);

	filename = folderPath + "/I_J00ion.csv";
	data_J00ion = Interpolation2d(filename, xTrans, yTrans, fTrans);

	filename = folderPath + "/I_J11rec.csv";
	data_J11rec = Interpolation2d(filename, xTrans, yTrans, fTrans);

	filename = folderPath + "/I_J22rec.csv";
	data_J22rec = Interpolation2d(filename, xTrans, yTrans, fTrans);

	filename = folderPath + "/I_J12rec.csv";
	data_J12rec = Interpolation2d(filename, xTrans, yTrans, fTrans);

	filename = folderPath + "/I_R0ion.csv";
	data_R0ion = Interpolation2d(filename, xTrans, yTrans, fTrans);

	filename = folderPath + "/I_R1rec.csv";
	data_R1rec = Interpolation2d(filename, xTrans, yTrans, fTrans);

	filename = folderPath + "/I_R2rec.csv";
	data_R2rec = Interpolation2d(filename, xTrans, yTrans, fTrans);
}

InelasticCollision::~InelasticCollision()
{
}

Eigen::VectorXd InelasticCollision::getGionInterpolated(const Eigen::VectorXd & Te)
{
	return Eigen::VectorXd();
}

Eigen::VectorXd InelasticCollision::getR0ionInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda)
{
	return data_R0ion.bicubicInterp(lambda, Te);
}

Eigen::VectorXd InelasticCollision::getJ00ionInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda)
{
	return data_J00ion.bicubicInterp(lambda, Te);
}

Eigen::VectorXd InelasticCollision::getGrecInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda)
{
	return data_Grec.bicubicInterp(lambda, Te);
}

Eigen::VectorXd InelasticCollision::getR1recInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda)
{
	return data_R1rec.bicubicInterp(lambda, Te);
}

Eigen::VectorXd InelasticCollision::getR2recInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda)
{
	return data_R2rec.bicubicInterp(lambda, Te);
}

Eigen::VectorXd InelasticCollision::getJ11recInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda)
{
	return data_J11rec.bicubicInterp(lambda, Te);
}

Eigen::VectorXd InelasticCollision::getJ22recInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda)
{
	return data_J22rec.bicubicInterp(lambda, Te);
}

Eigen::VectorXd InelasticCollision::getJ12recInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda)
{
	return data_J12rec.bicubicInterp(lambda, Te);
}

