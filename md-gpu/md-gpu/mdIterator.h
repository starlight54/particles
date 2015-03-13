#ifndef MD_ITERATOR_H
#define MD_ITERATOR_H

#include "iSimulationIterator.h"
#include "iForceEvaluator.h"
#include "iIntegrationEvaluator.h"
#include "visualiser.h"
#include <stdio.h>
#include <math.h>
#include <fstream>

//Nasty
#include "forceEvaluatorFactory.h"
#include "integrationEvaluatorFactory.h"
#include <cuda_runtime_api.h> //Particularly nasty, solely for IDE intellisense!

class MolDynIterator : public ISimulationIterator
{
public:
        MolDynIterator();
        ~MolDynIterator();
	void Print(ParticleSystem* particles);
	void Iterate(ParticleSystem* particles);
	void Initialise(unsigned long numberIterations)
	{};
	void Initialise(ParticleSystem* particles, unsigned long numberIterations
		, double temperature, double deltaT, double cutoff, double maxX,
		double maxY, double maxZ);
private:
	void UpdatePositions(ParticleSystem* particles);
	void UpdateForces(unsigned long numParticles);
	void UpdateVelocitiesT(unsigned long numParticles);
        void CreateVelocities();
        void CreateInitPrevPos();
	IForceEvaluator* forceEvaluator;
        IIntegrationEvaluator* integrator;

        double* prevPos;
        double* vel;
        double* force;

	double* devicePos;
	double* deviceVel;
	double* deviceForce;

	double maxX;
	double maxY;
	double maxZ;
	double cutoff;
	double xDist;
	double yDist;
	double zDist;
};

#endif //MD_ITERATOR_H