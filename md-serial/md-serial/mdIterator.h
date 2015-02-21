#ifndef MD_ITERATOR_H
#define MD_ITERATOR_H

#include "iSimulationIterator.h"
#include "forceEvaluator.h"
#include "integrationEvaluator.h"
#include <math.h>

class MolDynIterator : public ISimulationIterator
{
public:
        MolDynIterator();
        ~MolDynIterator();
	void Iterate(ParticleSystem* particles);
	void Initialise(unsigned long numberIterations)
	{};
	void Initialise(ParticleSystem* particles, unsigned long numberIterations
		, double temperature, double deltaT, double cutoff, double maxX,
		double maxY, double maxZ);
private:
        void CreateVelocities();
        void CreateInitPrevPos();
        ForceEvaluator* forceEvaluator;
        IntegrationEvaluator* integrator;
        double* prevPos;
        double* vel;
        double* force;
	double maxX;
	double maxY;
	double maxZ;
	double cutoffSquared;
};

#endif //MD_ITERATOR_H