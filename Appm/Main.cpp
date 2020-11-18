#include "Main.h"
#include "Interpolation2d.h"

// This library is not available in Intel C++ Compiler 2018.1.163
// #include <filesystem>
// namespace fs = std::experimental::filesystem;
//namespace fs = std::filesystem;

/**
* This is a test for checking bicubic spline interpolation
*/
//void test_bicubicInterpolation() {
//	Eigen::VectorXd lambdaVec(10);
//	lambdaVec << 0, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 1, 3, 5, 10;
//
//	const std::string filename = "collisions/inelastic/I_R0ion.csv";
//	const DataTransform xTrans = DataTransform::NONE;
//	const DataTransform yTrans = DataTransform::INVERSE;
//	const DataTransform fTrans = DataTransform::LOG;
//
//	const double fScale = 1;
//
//	Interpolation2d I_R0ion(filename, xTrans, yTrans, fTrans, fScale);
//
//	Eigen::VectorXd result;
//	Eigen::VectorXd xSites;
//	Eigen::VectorXd ySites = Eigen::VectorXd::LinSpaced(298,3e2,40e3);
//
//	Eigen::MatrixXd Mout(ySites.size(), lambdaVec.size() + 1);
//	Mout.col(0) = ySites;
//	for (int j = 0; j < lambdaVec.size(); j++) {
//		double lambda = lambdaVec(j);
//		xSites = lambda * Eigen::VectorXd::Ones(ySites.size());
//		result = I_R0ion.bicubicInterp(xSites, ySites);
//		Mout.col(j+1) = result;
//	}
//	std::ofstream("output.dat") << Mout << std::endl;
//}

void test_interpolation2d() {
	const std::string filename = "collisions/inelastic/I_Gion.csv";
	DataTransform xTrans = DataTransform::NONE;
	DataTransform yTrans = DataTransform::INVERSE;
	DataTransform fTrans = DataTransform::LOG;
	const double Tscale = 10e3;
	const double yScale = 1. / Tscale;
	const double sigmaScale = 1.36e-21;
	const double fScale = 1. / sigmaScale;
	Interpolation2d table2d(filename, xTrans, yTrans, fTrans, yScale, fScale);

	const int n = 101;
	Eigen::VectorXd xSites(n), ySites(n);
	ySites = Eigen::VectorXd::LinSpaced(n, 0.3, 4);
	xSites.array() = 0;

	Eigen::VectorXd fSites;
	fSites = table2d.bicubicInterp(xSites, ySites);

	Eigen::MatrixXd data(n, 3);
	data.col(0) = xSites;
	data.col(1) = ySites;
	data.col(2) = fSites;
	std::ofstream file("test_table2d_out.dat");
	file << data << std::endl;
}


void test_Gion() {
	std::string folderPath = "collisions/inelastic/";
	int idxA = -1;
	int idxE = -1;
	int idxI = -1;
	ScalingParameters scales("scales.txt");
	std::cout << scales << std::endl;
	InelasticCollision collision(folderPath, idxE, idxA, idxI, scales);

	Eigen::VectorXd Te = Eigen::VectorXd::LinSpaced(398, 300, 30e3);
	
	const int n = Te.size();
	Eigen::VectorXd nE = Eigen::VectorXd::Ones(n);
	Eigen::VectorXd nA = Eigen::VectorXd::Ones(n);
	Eigen::VectorXd nI = Eigen::VectorXd::Ones(n);
	Eigen::VectorXd vthE = Eigen::VectorXd::Ones(n);
	Eigen::VectorXd xStar = Eigen::VectorXd::Zero(n);
	const double mE = 1e-4;
	Eigen::VectorXd lambdaIon = Eigen::VectorXd::Zero(n);
	Eigen::VectorXd lambdaRec = Eigen::VectorXd::Zero(n);

	Eigen::VectorXd I_Gion = collision.getGionInterpolated(Te, lambdaIon);
	Eigen::VectorXd Gion = collision.getGion(nE, nA, vthE, Te, lambdaIon);
	Eigen::VectorXd Grec = collision.getGrec(nI, nE, vthE, xStar, mE, Te, lambdaRec);

	Eigen::MatrixXd data(n, 4);
	data.col(0) = Te;
	data.col(1) = I_Gion;
	data.col(2) = Gion;
	data.col(3) = Grec;

	std::cout << "\tTe\t" << "I_Gion\t" << "Gion\t" << "Grec\t" << std::endl;
	std::cout << data << std::endl;
	
	std::ofstream file("output.dat");
	file << "\tTe\t" << "I_Gion\t" << "Gion\t" << "Grec\t" << std::endl;
	file << data << std::endl;
}

// void removeFiles(const fs::path folder, const std::string extension) {
	// //const std::string extension = ".h5";
	// std::cout << "Remove files with extension " << extension << std::endl;
	// auto it = fs::directory_iterator(folder);
	// auto counter = 0;
	// for (auto & p : it) {
		// const std::string filename = p.path().string();
		// const std::string ext = filename.substr(filename.size() - extension.size());
		// if (ext == extension) {
			// //std::cout << filename << std::endl;
			// const bool isRemoved = fs::remove(p.path()); // delete file
			// assert(isRemoved);
			// counter++;
		// }
	// }
	// std::cout << "Number of deleted files: " << counter << std::endl;
// }


/**
* Main function.
*/
int main() {
	std::cout << "***********************" << std::endl;
	std::cout << "*    APPM             *" << std::endl;
	std::cout << "***********************" << std::endl;

	// Note: not available with Intel C++ compiler 2018.1.163
	// auto cwd = fs::current_path();
	// std::cout << "Current directory: " << cwd << std::endl;
	// removeFiles(cwd, ".h5");

	const std::string inputFilename = "input.txt";

	//test_Gion();
	//test_interpolation2d();
	//return EXIT_SUCCESS;
	//test_bicubicInterpolation();

	Main main;
	main.setInputFilename(inputFilename);
	main.run();
	std::cout << "TERMINATED" << std::endl;
	return EXIT_SUCCESS;
}


Main::Main()
{	
	std::cout << showVersionInfo();
	std::cout << std::endl;
	std::cout << std::endl;
}


Main::~Main()
{
}

/**
* Run the application. 
*
* Read the input file that specifies the mesh parameter file. 
* Then, read mesh parameters and pass them to the application.
*/
void Main::run()
{
	readInputFile();

	PrimalMesh::PrimalMeshParams primalMeshParams;
	assert(this->meshFilename.size() > 0);
	primalMeshParams = PrimalMesh::PrimalMeshParams(this->meshFilename);

	AppmSolver appm;
	appm.setScalingParameters("scales.txt");
	appm.setSolverParameters(solverParameters);
	appm.setMeshParameters(primalMeshParams);
	appm.setSpecies(speciesList);
	appm.setElasticCollisions(elasticCollisionList);
	appm.setInelasticCollisions(inelasticCollisionList);

	appm.run();	
}

void Main::setInputFilename(const std::string & filename)
{
	if (filename.size() == 0) {
		std::cout << "Input filename is too short (length = " << filename.size() << ")";
		exit(-1);
	}
	this->inputFilename = filename;
}

void Main::readInputFile()
{
	assert(inputFilename.size() > 0);
	std::cout << "Read input file: " << this->inputFilename << std::endl;
	std::ifstream file(this->inputFilename);
	if (!file.is_open()) {
		std::cout << "File could not be opened: " << inputFilename << std::endl;
		exit(-1);
	}
	const std::string speciesPrefix = "species";
	const std::string elasticCollisionPrefix = "elastic-collision";
	const std::string inelasticCollisionPrefix = "inelastic-collision";
	std::vector<std::string> collisionInputList;

	const char delim = ':';
	std::string line;
	while (std::getline(file, line)) {
		// Skip empty lines
		if (line.size() == 0) { continue; }
		// Skip comment lines
		if (line.front() == '#') { continue; }

		const int pos = line.find(delim); // position of delimiter
		if (pos == std::string::npos) {
			std::cout << "delimiter (" << delim << ") not found in text line: " << std::endl;
			std::cout << line << std::endl;
			exit(-1);
		}
		const std::string tag = line.substr(0, pos);    // tag
		const std::string value = line.substr(pos + 1); // value
		assert(value.size() > 0);

		// Apply tag-value pairs
		if (tag == "mesh") {
			std::istringstream(value) >> this->meshFilename;
		}
		if (tag == "maxTime") {
			double maxTime = 0;
			std::istringstream(value) >> maxTime;
			if (maxTime < 0) {
				maxTime = std::numeric_limits<int>::max();
			}
			this->solverParameters.setMaxTime(maxTime);
		}
		if (tag == "maxIterations") {
			int maxIterations = 0;
			std::istringstream(value) >> maxIterations;
			this->solverParameters.setMaxIterations(maxIterations);
		}
		if (tag == "outputFrequency") {
			int outputFrequency = 1;
			std::istringstream(value) >> outputFrequency;
			this->solverParameters.setOutputFrequency(outputFrequency);
		}
		if (tag == "isFluidEnabled") {
			bool b = false;
			std::istringstream(value) >> b;
			solverParameters.setFluidEnabled(b);
		}
		if (tag == "isMassfluxImplicit") {
			bool b = false;
			std::istringstream(value) >> b;
			solverParameters.setMassfluxSchemeImplicit(b);
		}
		if (tag == "isLorentzForceEnabled") {
			bool b = false;
			std::istringstream(value) >> b;
			solverParameters.setLorentzForceEnabled(b);
		}
		if (tag.substr(0,speciesPrefix.size()) == speciesPrefix) {
			speciesList.push_back(Species(value));
		}
		if (tag.substr(0, elasticCollisionPrefix.size()) == elasticCollisionPrefix) {
			this->elasticCollisionList.push_back(trim(value));
		}
		if (tag.substr(0, inelasticCollisionPrefix.size()) == inelasticCollisionPrefix) {
			this->inelasticCollisionList.push_back(trim(value));
		}
		if (tag == "initFluidState") {
			solverParameters.setFluidInitType(value);
		}
		if (tag == "initEfield") {
			double x, y, z;
			std::stringstream ss(value);
			ss >> x >> y >> z;
			solverParameters.setInitEfield(Eigen::Vector3d(x, y, z));
		}
		if (tag == "isEulerSourcesImplicit") {
			bool b = false;
			std::istringstream(value) >> b;
			solverParameters.setEulerSourcesImplicit(b);
		}
		if (tag == "timestepSizeFactor") {
			double dtFactor = 1;
			std::istringstream(value) >> dtFactor;
			solverParameters.setTimestepSizeFactor(dtFactor);
		}
		if (tag == "isMaxwellEnabled") {
			bool b = false;
			std::istringstream(value) >> b;
			solverParameters.setMaxwellEnabled(b);
		}
		if (tag == "voltageBCfilename") {
			solverParameters.setVoltageBCfilename(trim(value));
		}

	}

	std::cout << solverParameters << std::endl;

	std::cout << "Species: " << std::endl;
	std::cout << "======================" << std::endl;
	for (auto species : speciesList) {
		std::cout << species << std::endl;
		std::cout << "========" << std::endl;
	}
	std::cout << "======================" << std::endl;

	std::cout << std::endl;
	std::cout << "Elastic collisions: " << std::endl;
	std::cout << "======================" << std::endl;
	for (auto item : elasticCollisionList) {
		std::cout << item << std::endl;
		std::cout << "========" << std::endl;
	}
	std::cout << "======================" << std::endl;

	std::cout << std::endl;
	std::cout << "Inelastic collisions: " << std::endl;
	std::cout << "======================" << std::endl;
	for (auto item : inelasticCollisionList) {
		std::cout << item << std::endl;
		std::cout << "========" << std::endl;
	}
	std::cout << "======================" << std::endl;

	//exit(-1);
}

/**
* @return a text showing compile time and versions of used libraries.
*/
std::string Main::showVersionInfo()
{
	std::stringstream ss;
	ss << "Compiled on " << __TIMESTAMP__ << std::endl;
	ss << "with libraries: " << std::endl;
	
	// Version of Eigen library
	ss << Eigen::showVersion() << std::endl;

	// Version of MKL library
	const int len = 200;
	char mklVersion[len];
	mkl_get_version_string(mklVersion, len);
	ss << mklVersion << std::endl;

	// Version of HDF5 library
	ss << H5_VERS_INFO;
	return ss.str();
}

