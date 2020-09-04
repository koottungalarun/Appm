#include "InelasticCollision.h"

InelasticCollision::InelasticCollision()
	: idxE(0), idxA(0), idxI(0)
{
}

InelasticCollision::InelasticCollision(const std::string & folderPath, const int idxE_, const int idxA_, const int idxI_) 
	: idxE(idxE_), idxA(idxA_), idxI(idxI_)
{
	const DataTransform xTrans = DataTransform::NONE;
	const DataTransform yTrans = DataTransform::INVERSE;
	const DataTransform fTrans = DataTransform::LOG;

	assert(folderPath.back() == '/'); // folder path should end with a path separator

	std::string filename;
	filename = folderPath + "I_Gion.csv";
	data_Gion = Interpolation1d(filename, xTrans, fTrans);

	filename = folderPath + "I_Grec.csv";
	data_Grec = Interpolation2d(filename, xTrans, yTrans, fTrans);

	filename = folderPath + "I_J00ion.csv";
	data_J00ion = Interpolation2d(filename, xTrans, yTrans, fTrans);

	filename = folderPath + "I_J11rec.csv";
	data_J11rec = Interpolation2d(filename, xTrans, yTrans, fTrans);

	filename = folderPath + "I_J22rec.csv";
	data_J22rec = Interpolation2d(filename, xTrans, yTrans, fTrans);

	filename = folderPath + "I_J12rec.csv";
	// if lambda = 0 we have J12rec = 0, transformation with exp() cannot be executed. 
	// Therefore, pull lambda out of integral 
	data_J12rec = Interpolation2d(filename, xTrans, yTrans, fTrans); 

	filename = folderPath + "I_R0ion.csv";
	data_R0ion = Interpolation2d(filename, xTrans, yTrans, fTrans);

	filename = folderPath + "I_R1rec.csv";
	data_R1rec = Interpolation2d(filename, xTrans, yTrans, fTrans);

	filename = folderPath + "I_R2rec.csv";
	data_R2rec = Interpolation2d(filename, xTrans, yTrans, fTrans);
}

InelasticCollision::~InelasticCollision()
{
}

const Eigen::VectorXd InelasticCollision::getGionInterpolated(const Eigen::VectorXd & Te) const
{
	return data_Gion.cubicInterp(Te);
}

const Eigen::VectorXd InelasticCollision::getR0ionInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda) const
{
	return data_R0ion.bicubicInterp(lambda, Te);
}

const Eigen::VectorXd InelasticCollision::getJ00ionInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda) const
{
	return data_J00ion.bicubicInterp(lambda, Te);
}

const Eigen::VectorXd InelasticCollision::getGrecInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda) const
{
	return data_Grec.bicubicInterp(lambda, Te);
}

const Eigen::VectorXd InelasticCollision::getR1recInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda) const
{
	return data_R1rec.bicubicInterp(lambda, Te);
}

const Eigen::VectorXd InelasticCollision::getR2recInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda) const
{
	return data_R2rec.bicubicInterp(lambda, Te);
}

const Eigen::VectorXd InelasticCollision::getJ11recInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda) const
{
	return data_J11rec.bicubicInterp(lambda, Te);
}

const Eigen::VectorXd InelasticCollision::getJ22recInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda) const
{
	return data_J22rec.bicubicInterp(lambda, Te);
}

const Eigen::VectorXd InelasticCollision::getJ12recInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda) const
{
	return data_J12rec.bicubicInterp(lambda, Te);
}

const int InelasticCollision::getElectronFluidx() const
{
	return this->idxE;
}

const int InelasticCollision::getIonFluidx() const
{
	return this->idxI;
}

const int InelasticCollision::getAtomFluidx() const
{
	return this->idxA;
}

