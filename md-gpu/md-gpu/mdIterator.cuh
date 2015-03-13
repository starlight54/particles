#include "particleSystem.h"
#include "iIntegrationEvaluator.h"
#include "iForceEvaluator.h"

extern "C" {
	void UpdatePositions(ParticleSystem* particles);
	void UpdateForces(unsigned long numParticles);
	void UpdateVelocitiesT(unsigned long numParticles);
	__global__ void DeviceUpdatePositions(ParticleSystem* particles, IIntegrationEvaluator* integrator,
		double* devicePos, double* deviceVel, double* deviceForce, double deltaT, double maxX,
		double maxY, double maxZ, double &kinEn);
	__global__ void DeviceUpdateForces(unsigned long numParticles, double* deviceForce,
		double* devicePos, double maxX, double maxY, double maxZ, IForceEvaluator* forceEvaluator);
	__global__ void DeviceUpdateVelocitiesT(unsigned long numParticles, double* deviceForce,
		double deltaT, double* deviceVel);
}