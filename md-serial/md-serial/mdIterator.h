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
	IForceEvaluator* forceEvaluator;
        IIntegrationEvaluator* integrator;
        double* prevPos;
        double* vel;
        double* force;
	double maxX;
	double maxY;
	double maxZ;
	double cutoff;
	double xDist;
	double yDist;
	double zDist;
};

#endif //MD_ITERATOR_H