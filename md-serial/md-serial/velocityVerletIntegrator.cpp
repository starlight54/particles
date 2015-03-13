#include "velocityVerletIntegrator.h"

VelocityVerletIntegrator::VelocityVerletIntegrator()
{

}

VelocityVerletIntegrator::~VelocityVerletIntegrator()
{

}

double VelocityVerletIntegrator::Evaluate(double pos, double &vel, double force, 
	double deltaT)
{
	double forceT = force * deltaT * 0.5;
	vel += forceT;
	double newPos = pos + vel * deltaT + deltaT * forceT;

	return newPos;
}

void VelocityVerletIntegrator::Evaluate()
{

}
