#include "InelasticCollision.h"

InelasticCollision::InelasticCollision()
	: idxE(0), idxA(0), idxI(0)
{
}

InelasticCollision::InelasticCollision(const std::string & folderPath, 
	const int idxE_, const int idxA_, const int idxI_, 
	const ScalingParameters & scales) 
	: idxE(idxE_), idxA(idxA_), idxI(idxI_)
{
	this->sigmaScale = scales.getCrossSectionsScale();
	std::cout << "sigmaScale: " << sigmaScale << std::endl;
	const double sigmaScaleInverse = 1. / sigmaScale;

	// set parameters 
	const double ionizationEnergy_eV = 15.75;
	this->ionizationEnergyScaled = Physics::electronVolt_to_Joule(ionizationEnergy_eV) / scales.getEnergyScale();
	this->nScale = scales.getNumberDensityScale();
	this->g0 = 1;
	this->g1 = 6;

	PhysicsConstants & pc = PhysicsConstants::instance();
	const double pi = pc.pi();
	const double kB = pc.kB();
	const double h = pc.planckConstant();
	const double massScale = scales.getMassScale();
	const double Tscale = scales.getTemperatureScale();
	this->thermalDeBroglieScale = h / sqrt(2 * pi* massScale * kB * Tscale);
	assert(std::isfinite(thermalDeBroglieScale));
	assert(thermalDeBroglieScale > 0);


	const DataTransform xTrans = DataTransform::NONE;
	const DataTransform yTrans = DataTransform::INVERSE;
	const DataTransform fTrans = DataTransform::LOG;

	const double yScale = 1. / Tscale;

	assert(folderPath.back() == '/'); // folder path should end with a path separator

	std::string filename;
	filename = folderPath + "I_Gion.csv";
	data_Gion = Interpolation2d(filename, xTrans, yTrans, fTrans, yScale, sigmaScaleInverse);
	data_Gion.writeData("I_Gion_temp.dat", yScale, sigmaScaleInverse);

	filename = folderPath + "I_Grec.csv";
	data_Grec = Interpolation2d(filename, xTrans, yTrans, fTrans, yScale, sigmaScaleInverse);

	filename = folderPath + "I_J00ion.csv";
	data_J00ion = Interpolation2d(filename, xTrans, yTrans, fTrans, yScale, sigmaScaleInverse);

	filename = folderPath + "I_J11rec.csv";
	data_J11rec = Interpolation2d(filename, xTrans, yTrans, fTrans, yScale, sigmaScaleInverse);

	filename = folderPath + "I_J22rec.csv";
	data_J22rec = Interpolation2d(filename, xTrans, yTrans, fTrans, yScale, sigmaScaleInverse);

	filename = folderPath + "I_J12rec.csv";
	// if lambda = 0 we have J12rec = 0, transformation with exp() cannot be executed. 
	// Therefore, pull lambda out of integral 
	data_J12rec = Interpolation2d(filename, xTrans, yTrans, fTrans, yScale, sigmaScaleInverse);

	filename = folderPath + "I_R0ion.csv";
	data_R0ion = Interpolation2d(filename, xTrans, yTrans, fTrans, yScale, sigmaScaleInverse);

	filename = folderPath + "I_R1rec.csv";
	data_R1rec = Interpolation2d(filename, xTrans, yTrans, fTrans, yScale, sigmaScaleInverse);

	filename = folderPath + "I_R2rec.csv";
	data_R2rec = Interpolation2d(filename, xTrans, yTrans, fTrans, yScale, sigmaScaleInverse);
}

InelasticCollision::~InelasticCollision()
{
}

const Eigen::VectorXd InelasticCollision::getGion(const Eigen::VectorXd & nE, const Eigen::VectorXd & nA, const Eigen::VectorXd & vthE, const Eigen::VectorXd & TeVec, const Eigen::VectorXd & lambdaIon) const
{	
	const Eigen::VectorXd I_Gion = getGionInterpolated(TeVec, lambdaIon);
	const Eigen::VectorXd Gion = nE.array() * nA.array() * vthE.array() * (-1 * lambdaIon.array()).exp() * I_Gion.array();
	assert(Gion.allFinite());
	return Gion;
}

const Eigen::VectorXd InelasticCollision::getGrec(const Eigen::VectorXd & nI, const Eigen::VectorXd & nE, const Eigen::VectorXd & vthE, const Eigen::VectorXd & xStar, const double mE, const Eigen::VectorXd & Te, const Eigen::VectorXd & lambdaVec) const
{
	const double factor0 = nScale * pow(thermalDeBroglieScale, 3) * g0 / (2 * g1);
	const Eigen::VectorXd factor1 = (mE * Te).array().pow(-3./2.) * nI.array() * nE.array().pow(2) * vthE.array();
	const Eigen::VectorXd factor2 = (-2 * lambdaVec).array().exp() * xStar.array().exp();
	const Eigen::VectorXd I_Grec = getGrecInterpolated(Te, lambdaVec);
	const Eigen::VectorXd Grec = factor0 * factor1.array() * factor2.array() * I_Grec.array();
	const bool isGrecFinite = Grec.allFinite();
	if (!isGrecFinite) {
		std::cout << "is factor0 finite: " << std::isfinite(factor0) << std::endl;
		std::cout << "is factor1 finite: " << factor1.allFinite() << std::endl;
		std::cout << "is factor2 finite: " << factor2.allFinite() << std::endl;
		std::cout << "is I_Grec finite: " << I_Grec.allFinite() << std::endl;

		const int rows = nI.size();
		Eigen::MatrixXd data(rows, 4);
		data.col(0) = factor1;
		data.col(1) = factor2;
		data.col(2) = I_Grec;
		data.col(3) = Grec;
		std::ofstream("temp.dat") << data << std::endl;

		std::cout << "Data for cell idx = 0" << std::endl;
		std::cout << "nI: " << nI(0) << std::endl;
		std::cout << "nE: " << nE(0) << std::endl;
		std::cout << "vthE: " << vthE(0) << std::endl;
		std::cout << "xStar: " << xStar(0) << std::endl;
		std::cout << "mE: " << mE << std::endl;
		std::cout << "Te: " << Te(0) << std::endl;
		std::cout << "lambdaVec: " << lambdaVec(0) << std::endl;
		std::cout << "factor0: " << factor0 << std::endl;
		std::cout << "factor1: " << factor1(0) << std::endl;
		std::cout << "factor2: " << factor2(0) << std::endl;
		std::cout << "I_Grec: " << I_Grec(0) << std::endl;
		std::cout << "Grec: " << Grec(0) << std::endl;
	}
	assert(isGrecFinite);
	return Grec;
}

const Eigen::VectorXd InelasticCollision::getR0ion(const Eigen::VectorXd & nE, const Eigen::VectorXd & nA, const Eigen::VectorXd & vthE, const Eigen::VectorXd & TeVec, const Eigen::VectorXd & lambdaIon) const
{
	const Eigen::VectorXd I_R0ion = getR0ionInterpolated(TeVec, lambdaIon);
	const Eigen::VectorXd factor0 = 2./3. * nE.array() * nA.array() * vthE.array();
	const Eigen::VectorXd factor1 = (-1 * lambdaIon).array().exp() * I_R0ion.array();
	const Eigen::VectorXd R0ion = factor0.array() * factor1.array();
	assert(R0ion.allFinite());
	return R0ion;
}

const Eigen::VectorXd InelasticCollision::getJ00ion(const Eigen::VectorXd & nE, const Eigen::VectorXd & nA, const Eigen::VectorXd & vthE, const Eigen::VectorXd & TeVec, const Eigen::VectorXd & lambdaIon) const
{
	const Eigen::VectorXd I_J00ion = getJ00ionInterpolated(TeVec, lambdaIon);
	const Eigen::VectorXd factor0 = 2. / 3. * nE.array() * nA.array() * vthE.array();
	const Eigen::VectorXd factor1 = (-1 * lambdaIon).array().exp() * I_J00ion.array();
	const Eigen::VectorXd J00ion = factor0.array() * factor1.array();
	assert(J00ion.allFinite());
	return J00ion;
}

const Eigen::VectorXd InelasticCollision::getR1rec(const Eigen::VectorXd & nE, const Eigen::VectorXd & nI, const Eigen::VectorXd & vthE, const Eigen::VectorXd & xStar, const Eigen::VectorXd & TeVec, const Eigen::VectorXd & lambdaRec) const
{
	const Eigen::VectorXd I_R1rec = getR1recInterpolated(TeVec, lambdaRec);
	const double factor0 = 2. / 3. * (0.5 * g0) / g1 * nScale * pow(thermalDeBroglieScale, 3);
	const Eigen::VectorXd factor1 = nI.array() * nE.array().pow(2) * vthE.array();
	const Eigen::VectorXd factor2 = (-2 * lambdaRec).array().exp() * xStar.array().exp() * I_R1rec.array();
	const Eigen::VectorXd R1rec = factor0 * factor1.array() * factor2.array();
	assert(R1rec.allFinite());
	return R1rec;
}

const Eigen::VectorXd InelasticCollision::getR2rec(const Eigen::VectorXd & nE, const Eigen::VectorXd & nI, const Eigen::VectorXd & vthE, const Eigen::VectorXd & xStar, const Eigen::VectorXd & TeVec, const Eigen::VectorXd & lambdaRec) const
{
	const Eigen::VectorXd I_R2rec = getR2recInterpolated(TeVec, lambdaRec);
	const double factor0 = 2. / 3. * (0.5 * g0) / g1 * nScale * pow(thermalDeBroglieScale, 3);
	const Eigen::VectorXd factor1 = nI.array() * nE.array().pow(2) * vthE.array();
	const Eigen::VectorXd factor2 = (-2 * lambdaRec).array().exp() * xStar.array().exp() * I_R2rec.array();
	const Eigen::VectorXd R2rec = factor0 * factor1.array() * factor2.array();
	assert(R2rec.allFinite());
	return R2rec;
}

const Eigen::VectorXd InelasticCollision::getJ11rec(const Eigen::VectorXd & nE, const Eigen::VectorXd & nI, const Eigen::VectorXd & vthE, const Eigen::VectorXd & xStar, const Eigen::VectorXd & TeVec, const Eigen::VectorXd & lambdaRec) const
{
	const Eigen::VectorXd I_J11rec = getJ11recInterpolated(TeVec, lambdaRec);
	const double factor0 = (0.5 * g0) / g1 * nScale * pow(thermalDeBroglieScale, 3);
	const Eigen::VectorXd factor1 = nI.array() * nE.array().pow(2) * vthE.array();
	const Eigen::VectorXd factor2 = (-2 * lambdaRec).array().exp() * xStar.array().exp() * I_J11rec.array();
	const Eigen::VectorXd R2rec = factor0 * factor1.array() * factor2.array();
	assert(R2rec.allFinite());
	return R2rec;
}

const Eigen::VectorXd InelasticCollision::getJ22rec(const Eigen::VectorXd & nE, const Eigen::VectorXd & nI, const Eigen::VectorXd & vthE, const Eigen::VectorXd & xStar, const Eigen::VectorXd & TeVec, const Eigen::VectorXd & lambdaRec) const
{
	const Eigen::VectorXd I_J22rec = getJ22recInterpolated(TeVec, lambdaRec);
	const double factor0 = (0.5 * g0) / g1 * nScale * pow(thermalDeBroglieScale, 3);
	const Eigen::VectorXd factor1 = nI.array() * nE.array().pow(2) * vthE.array();
	const Eigen::VectorXd factor2 = (-2 * lambdaRec).array().exp() * xStar.array().exp() * I_J22rec.array();
	const Eigen::VectorXd R2rec = factor0 * factor1.array() * factor2.array();
	assert(R2rec.allFinite());
	return R2rec;
}

const Eigen::VectorXd InelasticCollision::getJ12rec(const Eigen::VectorXd & nE, const Eigen::VectorXd & nI, const Eigen::VectorXd & vthE, const Eigen::VectorXd & xStar, const Eigen::VectorXd & TeVec, const Eigen::VectorXd & lambdaRec) const
{
	const Eigen::VectorXd I_J12rec = getJ12recInterpolated(TeVec, lambdaRec);
	const double factor0 = 4./9. * (0.5 * g0) / g1 * nScale * pow(thermalDeBroglieScale, 3);
	const Eigen::VectorXd factor1 = nI.array() * nE.array().pow(2) * vthE.array();
	const Eigen::VectorXd factor2 = (-2 * lambdaRec).array().exp() * xStar.array().exp() * lambdaRec.array() * I_J12rec.array();
	const Eigen::VectorXd R2rec = factor0 * factor1.array() * factor2.array();
	assert(R2rec.allFinite());
	return R2rec;
}

const Eigen::VectorXd InelasticCollision::getGionInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda) const
{
	return data_Gion.bicubicInterp(lambda, Te);
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

const double InelasticCollision::getIonizationEnergyScaled() const
{
	return this->ionizationEnergyScaled;
}

