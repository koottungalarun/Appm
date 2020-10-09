#pragma once

#define _USE_MATH_DEFINES
#define EIGEN_USE_MKL_ALL

#include <cmath>
#include <iostream>
#include <iomanip>
#include <exception>
#include <string>

#include "PrimalMesh.h"
#include "DualMesh.h"
#include "Numerics.h"

#include "Physics.h"
#include "ScalingParameters.h"
#include "Species.h"

#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers> 	
#include <Eigen/PardisoSupport>

#include <chrono>
#include <stdexcept>

#include "ElasticCollision.h"
#include "InelasticCollision.h"

#define _RT_ONECELL 
#undef _RT_ONECELL


class AppmSolver
{

public:

	enum class MaxwellSolverType {
		CG, PardisoLU, BiCGStab
	};
	friend std::ostream & operator<<(std::ostream & os, const AppmSolver::MaxwellSolverType & obj);

	enum class FluidInitType {
		DEFAULT, 
		UNIFORM, 
		SHOCKTUBE, 
		TEST_FRICTION, 
		TEST_FRICTION_TEMPERATURE, 
		TEST_FRICTION_NUMBERDENSITY, 
		TEST_FRICTION_ELECTRONS_NONZERO_VELOCITY,
		INIT_FILE
	};
	friend std::ostream & operator<<(std::ostream & os, const AppmSolver::FluidInitType & obj);

	class SolverParameters {
	public:
		SolverParameters();
		~SolverParameters();

		void setMaxIterations(const int n);
		const int getMaxIterations() const;
		void setMaxTime(const double tmax);
		const double getMaxTime() const;
		void setMaxwellEnabled(const bool b);
		const bool getMaxwellEnabled() const;
		void setMaxTimestepSize(const double dtMax);
		const double getMaxTimestepSize() const;
		void setFluidEnabled(const bool b);
		const bool getFluidEnabled() const;
		void setLorentzForceEnabled(const bool b);
		const bool getLorentzForceEnabled() const;
		void setMassfluxSchemeImplicit(const bool b);
		const bool getMassfluxSchemeImplicit() const;
		void setFrictionActive(const bool b);
		const bool getFrictionActive() const;
		void setMaxwellCurrentDefined(const bool b);
		const bool getMaxwellCurrentDefined() const;
		void setMaxwellSolverType(const AppmSolver::MaxwellSolverType type);
		const AppmSolver::MaxwellSolverType getMaxwellSolverType() const;
		void setEulerMaxwellCouplingEnabled(const bool b);
		const bool getEulerMaxwellCouplingEnabled() const;
		void setOutputFrequency(const int n);
		const int getOutputFrequency() const;
		void setFluidInitType(const std::string & s);
		const AppmSolver::FluidInitType getFluidInitType() const;
		void setInitEfield(const Eigen::Vector3d & efield);
		const Eigen::Vector3d getInitEfield() const;
		void setApParameter(const double lambdaSquare);
		const double getApParameter() const;
		void setEulerSourcesImplicit(const bool b);
		const bool getEulerSourcesImplicit() const;
		const double getTimestepSizeFactor() const;
		void setTimestepSizeFactor(const double dt_factor);

		friend std::ostream & operator<<(std::ostream & os, const AppmSolver::SolverParameters & obj);
	private:
		double lambdaSq = 1;
		int maxIterations = 0;
		double maxTime = 0;
		double timestepSizeFactor = 1;
		double timestepSizeMax = 1;
		int outputFrequency = 1;

		bool isEulerMaxwellCouplingEnabled = false;
	
		bool isFluidEnabled = false;
		bool isMassFluxSchemeImplicit = true;
		bool isFrictionEnabled = false;
		bool isLorentzForceEnabled = false;
		bool isEulerSourcesImplicit = false;

		bool isMaxwellEnabled = false;
		MaxwellSolverType maxwellSolverType = MaxwellSolverType::BiCGStab;		
		AppmSolver::FluidInitType fluidInitType = AppmSolver::FluidInitType::UNIFORM;

		Eigen::Vector3d initEfield = Eigen::Vector3d::Zero();
	};


	AppmSolver();
	//AppmSolver(const PrimalMesh::PrimalMeshParams & primalMeshParams, 
	//	const AppmSolver::SolverParameters & appmParams);
	~AppmSolver();

	void run();

	void setSolverParameters(const AppmSolver::SolverParameters & solverParams);
	void setMeshParameters(const PrimalMesh::PrimalMeshParams & meshParams);
	void setSpecies(const std::vector<Species> & speciesList);
	void setElasticCollisions(const std::vector<std::string> & list);
	void setInelasticCollisions(const std::vector<std::string> & list);
	void setScalingParameters(const std::string & filename);


	// Structure to save init data
	struct InitDataStruct {
		std::string fluidName;
		double n = 0;
		double T = 0;
		double uz = 0;

	};
	friend std::ostream & operator<<(std::ostream & os, const AppmSolver::InitDataStruct & obj);


protected:
	PrimalMesh primalMesh;
	DualMesh dualMesh;

	Eigen::Matrix3Xd B_vertex;

	bool isMaxwellCurrentSource = false;

private:
	PrimalMesh::PrimalMeshParams primalMeshParams;
	SolverParameters solverParams;
	std::vector<Species> speciesList;
	ScalingParameters scalingParameters;

	//struct ParticleParameters {
	//	std::string name = "neutral";
	//	double mass = 1.0;
	//	int electricCharge = 0;
	//};
	//std::vector<ParticleParameters> particleParams;

	// List of elastic collisions
	std::vector<ElasticCollision*> elasticCollisions;
	void setElasticCollisionSourceTerms();

	std::vector<InelasticCollision*> inelasticCollisions;

	//struct AppmParameters {
	//	int maxIterations = 0;
	//	double maxTime = 0;
	//	bool isFluidEnabled = false;
	//	bool isMaxwellEnabled = false;
	//	bool isEulerMaxwellCouplingEnabled = false;
	//	bool isLorentzForceElectricEnabled = false;
	//	bool isLorentzForceMagneticEnabled = false;
	//	bool isMassFluxSchemeImplicit = false;
	//	double timestepSize = 1;
	//	bool isMaxwellCurrentDefined = false;
	//	bool isFrictionActive = true;
	//	MaxwellSolverType maxwellSolverType = MaxwellSolverType::CG;
	//} appmParams;

	std::ofstream timeFile;

	void init();
	std::string getIterationHeader(const int iter, const double time, const double dt) const;

	void debug_checkCellStatus() const;

	/** Get number of fluids. */
	const int getNFluids() const;

	/** */
	const Eigen::VectorXd getState(const int cIdx, const int fluidx) const;

	bool isStateWrittenToOutput = false;


	//double lambdaSquare = 1.0;
	//int initType = 1;
	//bool isElectricLorentzForceActive = false;
	//bool isMagneticLorentzForceActive = false;
	bool isMomentumFluxActive = false;

	bool isWriteEfield = false;
	bool isWriteBfield = false;
	bool isWriteHfield = false;
	
	// State vector for solving Maxwell's equations
	Eigen::VectorXd maxwellState;
	Eigen::VectorXd maxwellStatePrevious;
	Eigen::SparseMatrix<double> Q;
	Eigen::SparseMatrix<double> Meps;
	Eigen::SparseMatrix<double> Mnu;
	Eigen::SparseMatrix<double> C;

	Eigen::VectorXd E_h;
	Eigen::VectorXd B_h;
	Eigen::VectorXd J_h;
	Eigen::VectorXd J_h_previous;
	Eigen::VectorXd J_h_aux;
	Eigen::VectorXd J_h_aux_mm1;
	Eigen::Matrix3Xd Jcc;
	Eigen::Matrix3Xd Jaux_cc;
	Eigen::Matrix3Xd E_cc; // Electric field at cell center

	// M1 = lambda^2 * Q' * Meps * Q in the reformulated Ampere equation. 
	Eigen::SparseMatrix<double> M1;

	// M2 = Cdual * Mnu * C in the reformulated Ampere equation, restricted to inner edges. 
	Eigen::SparseMatrix<double> M2;


	// Fluid data vectors
	Eigen::MatrixXd fluidStates;
	Eigen::MatrixXd faceFluxes;
	Eigen::MatrixXd sumOfFaceFluxes;
	Eigen::MatrixXd fluidSources;
	Eigen::MatrixXd LorentzForce_magnetic;
	Eigen::MatrixXd LorentzForce_electric;
	Eigen::MatrixXd frictionForceSourceTerm;
	Eigen::MatrixXd frictionEnergySourceTerm;
	Eigen::MatrixXd diffusionVelocity;
	Eigen::MatrixXd bulkVelocity;
	Eigen::MatrixXd massFluxImplicitTerm;

	Eigen::MatrixXd faceFluxesImExRusanov;

	//std::string printSolverParameters() const;

	// Isentropic expansion coefficient, aka ratio of heat capacities
	//const double gamma = 1.4; 

	/* scaling factor for effectively removing truncation errors in the low bits */
	const double fluxTruncationErrorGuardScale = 1e4; 


	void init_maxwellStates();
	Eigen::SparseMatrix<double> getBoundaryGradientInnerInclusionOperator();
	Eigen::SparseMatrix<double> getElectricPermittivityOperator();
	Eigen::SparseMatrix<double> getMagneticPermeabilityOperator();

	void init_multiFluid();
	void applyFluidInitializationType();

	void init_SodShockTube(const double zRef);
	void init_Uniformly(const double n, const double p, const Eigen::Vector3d u);
	void init_Explosion(const Eigen::Vector3d refPos, const double radius);
	void init_testcase_frictionTerm();
	void init_ignitionWire();
	void init_fluid_frictonTest();
	void init_fluid_frictionTest_temperature();
	void init_fluid_frictionTest_numberDensity();
	void init_fluid_frictionTest_electrons_nonzero_velocity();

	void init_multiFluid_readFromFile(const std::string & filename);

	void set_Efield_uniform(const Eigen::Vector3d direction);
	void set_Bfield_azimuthal();

	Eigen::SparseMatrix<double> M_perot;
	Eigen::SparseMatrix<double> initPerotInterpolationMatrix();
	const Eigen::Matrix3Xd getEfieldAtCellCenter();
	const Eigen::Matrix3Xd getCurrentDensityAtCellCenter();

	//void get_Msigma_consistent(const double dt, Eigen::SparseMatrix<double> & Msigma, Eigen::VectorXd & jaux);
	
	void setRadiationSource();
	void setFrictionSourceTerms();
	void setMagneticLorentzForceSourceTerms();

	const int getFluidStateLength() const;
	const double getNextFluidTimestepSize() const;

	//const Eigen::VectorXd getFluidState(const int cellIdx, const int fluidIdx) const;
	const Eigen::Vector3d getFluidState(const int cellIdx, const int fluidIdx, const Eigen::Vector3d & faceNormal) const;
	
	const int getOrientation(const Cell * cell, const Face * face) const;
	
	//const Eigen::Vector3d getFluidFluxFromState(const Eigen::Vector3d & q) const;
	//const Eigen::Vector3d getRusanovFluxExplicit(const int faceIdx, const int fluidIdx) const;
	//const Eigen::Vector3d getRusanovFluxImEx(const int faceIdx, const int fluidIdx, const double dt);

	//const double getImplicitExtraTermMomentumFlux(const int cellIdx, const Eigen::Vector3d & faceNormal, const int fluidIdx) const;

	const std::pair<int,int> getAdjacientCellStates(const Face * face, const int fluidIdx, Eigen::Vector3d & qL, Eigen::Vector3d & qR) const;

	//const double getMomentumUpdate(const int k, const Eigen::Vector3d & nvec, const int fluidIdx) const;
	const std::string stopFilename = "stop.txt";
	void createStopFile(const int value);
	bool isStopFileActive();

	void setFluidFaceFluxes();
	Eigen::Vector3d getSpeciesFaceFlux(const Face * face, const int fluidIdx);
	const Eigen::Vector3d getSpeciesFaceFluxAtCathode(const Face * face, const int fluidIdx);
	//void setFluidSourceTerm();
	void setSumOfFaceFluxes();
	void setImplicitMassFluxTerms(const double dt);
	void updateFluidStates(const double dt, const bool isImplicitSources);
	Eigen::SparseMatrix<double> getJacobianEulerSourceElasticCollisions() const;
	Eigen::SparseMatrix<double> getJacobianEulerSourceInelasticCollisions(Eigen::VectorXd & rhs) const;
	Eigen::MatrixXd getInelasticSourcesExplicit();

	void solveMaxwellSystem(const double time, const double dt, const double dt_previous, const Eigen::SparseMatrix<double> & Msigma);
	

	void interpolateMagneticFluxToPrimalVertices();

	// void test_raviartThomas();
	const Eigen::Matrix3Xd getPrismReferenceCoords(const int nSamples);

	// Timestamps and iterations at which output data has been written
	std::vector<double> timeStamps;
	std::vector<int> outputIterations;

	void init_meshes(const PrimalMesh::PrimalMeshParams & primalParams);

	void writeXdmf(const std::string & filename);
	void writeXdmfDualVolume(const std::string & filename);
	
	bool isShowDataWriterOutput = false;
	void writeOutput(const int iteration, const double time);

	void writeFluidStates(H5Writer & writer);
	void writeMaxwellStates(H5Writer & writer);

	XdmfGrid getOutputPrimalEdgeGrid(const int iteration, const double time, const std::string & dataFilename);
	XdmfGrid getOutputPrimalSurfaceGrid(const int iteration, const double time, const std::string & dataFilename);
	
	XdmfGrid getOutputDualEdgeGrid(const int iteration, const double time, const std::string & dataFilename);
	XdmfGrid getOutputDualSurfaceGrid(const int iteration, const double time, const std::string & dataFilename);
	XdmfGrid getOutputDualVolumeGrid(const int iteration, const double time, const std::string & dataFilename);

	void init_RaviartThomasInterpolation();

	std::vector<Eigen::Matrix3d> rt_piolaMatrix;
	std::vector<Eigen::Vector3d> rt_piolaVector;

	//void readParameters(const std::string & filename);

	const std::string xdmf_GridPrimalEdges(const int iteration) const;
	const std::string xdmf_GridPrimalFaces(const int iteration) const;
	const std::string xdmf_GridDualEdges(const int iteration) const;
	const std::string xdmf_GridDualFaces(const int iteration) const;
	const std::string xdmf_GridDualCells(const int iteration) const;

	const std::string fluidXdmfOutput(const std::string & filename) const;

	Eigen::MatrixXi faceTypeFluids;
	const Face::Type getFaceTypeOfFluid(const Face * face, const int fluidIdx) const;

	const double getNumericalSchemeFactor(const Face * face, const int fluidx) const;

	const Eigen::VectorXd solveMaxwell_PardisoLU(Eigen::SparseMatrix<double> & Mf, Eigen::VectorXd & rhs);
	const Eigen::VectorXd solveMaxwell_sparseLU(Eigen::SparseMatrix<double> & Mf, Eigen::VectorXd & rhs);
	const Eigen::VectorXd solveMaxwell_BiCGStab(Eigen::SparseMatrix<double> & Mf, Eigen::VectorXd & rhs);
	const Eigen::VectorXd solveMaxwell_LSCG(Eigen::SparseMatrix<double> & Mf, Eigen::VectorXd & rhs);
	const Eigen::VectorXd solveMaxwell_CG(Eigen::SparseMatrix<double> & Mf, Eigen::VectorXd & rhs);


	Eigen::SparseMatrix<double> get_Msigma_spd(Eigen::VectorXd & Jaux, const double dt, const double time);

	const Eigen::VectorXd testcase_001_FluidSourceTerm(const double time, const Cell * cell, const int fluidIdx) const;
	const Eigen::VectorXd setVoltageBoundaryConditions(const double time) const;

	const Species & getSpecies(const int idx) const;
	const int getSpeciesIndex(const std::string & tag);

	const double getCollisionFrequency(const int alpha, const int beta, const int cellIdx);
	const double getReducedMass(const int alpha, const int beta) const;
	const Eigen::MatrixXd getStates(const int fidx, const int nCols) const;


	const int getLinearIndexInJacobian(const int fluidIdx, const int cellidx) const;

	const std::string message_howToVisualizeData(const std::string & outFilename_volume, const std::string & outFilename_surface) const;

};

/**
* Trim white space of a string.
* @see https://stackoverflow.com/a/17976541
*/
inline std::string trim(const std::string &s)
{
	auto  wsfront = std::find_if_not(s.begin(), s.end(), [](int c) {return std::isspace(c); });
	return std::string(wsfront, std::find_if_not(s.rbegin(), std::string::const_reverse_iterator(wsfront), [](int c) {return std::isspace(c); }).base());
}

