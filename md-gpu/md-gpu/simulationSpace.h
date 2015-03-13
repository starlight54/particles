#ifndef SIMULATION_SPACE_H
#define SIMULATION_SPACE_H

#include "particleSystem.h"
#include "iSimulationIterator.h"
#include "simulationIteratorFactory.h"

class SimulationSpace
{
public:
        SimulationSpace();
        ~SimulationSpace();
        void InitCube(double r);
        void InitCuboid(double x, double y, double z);
        void LoadParticles(ParticleSystem* particles);
        void AutoInitParticles(unsigned long numberParticles);
        void InitIterator();
	void InitIterator(double time, double deltaT, double sigma);
	void ExecuteSimulation();
private:
        double maxX;
        double maxY;
        double maxZ;
        ParticleSystem* particles;
        ISimulationIterator* simulationIterator;
};

#endif //SIMULATION_SPACE_H
