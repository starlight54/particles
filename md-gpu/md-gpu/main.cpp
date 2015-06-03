/*Simple serial implementation of MD with paramterised inputs
 *Michael J. Craig 15/02/15*/

#include "factoryRegister.h"
#include "simulationSpace.h"
#include "simulationIteratorFactory.h"
#include "configManager.h"

int main(int argc, char* argv[])
{
	//Register factory classes
	FactoryRegister::GlobalFactoryRegister();

	ConfigManager::Get()->Initialise();

	Visualiser::Get()->GlutInit(argc, argv);

        SimulationSpace space = SimulationSpace();
        space.InitCube(100);
	space.AutoInitParticles(10000);
	//system("pause");
	space.InitIterator(100, 0.005, 1);
}