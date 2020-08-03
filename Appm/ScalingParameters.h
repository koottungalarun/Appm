#pragma once

#include <cmath>
#include <cassert>
#include "PhysicsConstants.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

class ScalingParameters
{
private:
	/** Temperature scale (K) */
	double T0 = 0;
	/** Number density scale (m^-3) */
	double n0 = 0;
	/** Spatial length scale (m) */
	double x0 = 0;
	/** Time scale (s) */
	double t0 = 0;
	/** Energy scale (J) */
	double energy0 = 0;
	/** E-field scale (V m^-1) */
	double E0 = 0;
	/** pressure scale (Pa) */
	double p0 = 0;

	/* Scaled Debye length squared */
	double lambdaSq = 0;

	/** 
	* Get number of undefined scales in the set of T0, n0, x0, lambdaSq.
	*/
	const int undefinedScales();


public:
	ScalingParameters();
	ScalingParameters(const std::string & filename);
	~ScalingParameters();

	void applyIdealGasLaw();
	void compute();

	void setTemperatureScale(const double T);
	void setNumberDensityScale(const double n);
	void setLengthScale(const double x);
	void setScaledDebyeLength(const double lambdasq);
	void setPressureScale(const double p);

	double getTemperatureScale() const;
	double getLengthScale() const;
	double getNumberDensityScale() const;
	double getScaledDebyeLengthSquared() const;
	double getCrossSectionsScale() const;

	friend std::ostream & operator<<(std::ostream & os, const ScalingParameters & obj);

};

