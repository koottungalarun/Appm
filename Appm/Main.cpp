#include "Main.h"

/**
* Main function.
*/
int main() {
	std::cout << "***********************" << std::endl;
	std::cout << "*    APPM             *" << std::endl;
	std::cout << "***********************" << std::endl;

	//int N = 10;
	//Eigen::VectorXd x(N), y(N);
	//const double dx = 1. / (N - 1);
	//for (int i = 0; i < N; i++) {
	//	x(i) = i * dx;
	//}
	//y = x.array().pow(2);

	//InterpolationTable table(x, y);
	//Eigen::VectorXd sites(N);
	//sites = Eigen::VectorXd::LinSpaced(N, 0.0, 1.0);
	//sites(0) = std::nanf("");
	//Eigen::VectorXd result;
	//Eigen::MatrixXd data(sites.rows(), 2);

	//result = table.interpolate(sites);
	//data.col(0) = y;
	//data.col(1) = result;
	//std::cout << data << std::endl;

	//result = table.interpolate(sites);
	//data.col(0) = y;
	//data.col(1) = result;
	//std::cout << data << std::endl;
	//return EXIT_SUCCESS;

	const std::string inputFilename = "input.txt";

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
	//primalMeshParams = PrimalMesh::PrimalMeshParams("primalMeshParams.txt");
	assert(this->meshFilename.size() > 0);
	primalMeshParams = PrimalMesh::PrimalMeshParams(this->meshFilename);

	//AppmSolver appm(primalMeshParams, appmParameters);
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

