/*Simple serial implementation of MD with paramterised inputs
 *Michael J. Craig 15/02/15*/

#include "factoryRegister.h"
#include "simulationSpace.h"
#include "simulationIteratorFactory.h"

int main(int argc, char* argv[])
{
	//Register factory classes
	FactoryRegister::GlobalFactoryRegister();
	Visualiser::Get()->GlutInit(argc, argv);

        SimulationSpace space = SimulationSpace();
        space.InitCube(100);
	space.AutoInitParticles(2);
	space.InitIterator(200, 0.005, 1);
}