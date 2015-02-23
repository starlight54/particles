#include "simulationSpace.h"
#include "simulationIteratorFactory.h"

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

void SimulationSpace::InitIterator()
{
	simulationIterator = SimulationIteratorFactory::Get()->
		Create("MolDynIterator");
	this->ExecuteSimulation();
}

void SimulationSpace::ExecuteSimulation()
{
	simulationIterator->Initialise(particles, 5000, 273, 0.003, 2.5, maxX, maxY, maxZ);
	simulationIterator->Iterate(particles);
}
