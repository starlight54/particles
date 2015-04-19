#include "simulationSpace.h"

SimulationSpace::SimulationSpace()
{

}

SimulationSpace::~SimulationSpace()
{
        free(particles);
}

void SimulationSpace::InitCube(double r)
{
        maxX = maxY = maxZ = r;
}

void SimulationSpace::InitCuboid(double x, double y, double z)
{

}

void SimulationSpace::LoadParticles(ParticleSystem* particles)
{
        this->particles = particles;
}

void SimulationSpace::AutoInitParticles(unsigned long numberParticles)
{
        particles = new ParticleSystem();
        particles->AutoInit(maxX, maxY, maxZ, numberParticles);
}

void SimulationSpace::InitIterator(double time, double deltaT, double sigma)
{
	simulationIterator = SimulationIteratorFactory::Get()->
		Create("MolDynIterator");
	unsigned long numIterations = time / deltaT;
	simulationIterator->Initialise(particles, numIterations, 273, deltaT, 
		2.5, maxX, maxY, maxZ);
	for (int i = 0; i < numIterations; ++i) {
		simulationIterator->Iterate(particles);
	}

	simulationIterator->Print(particles);
}

void SimulationSpace::ExecuteSimulation()
{

}

//1.12246204831