#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include "Interpolation1d.h"
#include "Interpolation2d.h"
#include "DataTransform.h"
#include <fstream>
#include <Eigen/Dense>

class InelasticCollision
{
public:
	InelasticCollision();

	/* Read data for inelastic collisions from folder */
	InelasticCollision(const std::string & folderPath);

	~InelasticCollision();

	Eigen::VectorXd getGionInterpolated(const Eigen::VectorXd & Te);
	Eigen::VectorXd getR0ionInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda);
	Eigen::VectorXd getJ00ionInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda);
	Eigen::VectorXd getGrecInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda);
	Eigen::VectorXd getR1recInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda);
	Eigen::VectorXd getR2recInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda);
	Eigen::VectorXd getJ11recInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda);
	Eigen::VectorXd getJ22recInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda);
	Eigen::VectorXd getJ12recInterpolated(const Eigen::VectorXd & Te, const Eigen::VectorXd & lambda);



private:
	Interpolation1d data_Gion;
	Interpolation2d data_R0ion;
	Interpolation2d data_J00ion;

	Interpolation2d data_Grec;
	Interpolation2d data_R1rec;
	Interpolation2d data_R2rec;
	Interpolation2d data_J11rec;
	Interpolation2d data_J22rec;
	Interpolation2d data_J12rec;

};

