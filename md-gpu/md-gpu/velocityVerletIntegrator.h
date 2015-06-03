#ifndef VELOCITY_VERLET_INTEGRATOR_H
#define VELOCITY_VERLET_INTEGRATOR_H

#include "iIntegrationEvaluator.h"

class VelocityVerletIntegrator : public IIntegrationEvaluator
{
public:
	VelocityVerletIntegrator();
	~VelocityVerletIntegrator();
	void Evaluate();
	__device__ double Evaluate(double pos, double &vel, double force, double deltaT);
};

#endif