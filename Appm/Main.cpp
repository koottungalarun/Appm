#include "Main.h"

int main() {
	std::cout << "***********************" << std::endl;
	std::cout << "*    APPM             *" << std::endl;
	std::cout << "***********************" << std::endl;
	Main main;
	main.run();
	std::cout << "TERMINATED" << std::endl;
}


Main::Main()
{
}


Main::~Main()
{
}

void Main::run()
{
	AppmSolver appm;
	appm.run();
}
