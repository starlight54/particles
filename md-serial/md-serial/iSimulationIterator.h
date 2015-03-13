#ifndef SIMULATION_ITERATOR_H
#define SIMULATION_ITERATOR_H

#include "particleSystem.h"

class ISimulationIterator
{
public:
        ISimulationIterator();
        virtual ~ISimulationIterator();
        virtual void Iterate(ParticleSystem* particles) = 0;
	virtual void Print(ParticleSystem* particles) = 0;
	virtual void Initialise(unsigned long numberIterations) = 0;
	virtual void Initialise(ParticleSystem* particles,
		unsigned long numberIterations, double temperature, double deltaT, 
		double cutoff, double maxX, double maxY, double maxZ) = 0;
protected:
	//Common to all iterators?
	unsigned long numberIterations;
	float simulationTime;
	float deltaT;
	float nIters;
	float sigma;
	float epsilon;

        //Storage variables for outputs
        //Centre of Mass Velocity
        double comVel[3];
        //Kinetic Energy
        double kinEn;
        //Total Energy
        double totEn;
	//Scaled Total Energy
	double scaledTotEn;
	//Instantaneous Temp
	double instantTemp;
	//Energy per particle
	double energyPerParticle;
};

#endif //SIMULATION_ITERATOR_H