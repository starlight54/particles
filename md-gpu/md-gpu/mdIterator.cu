#include "mdIterator.h"
#include "mdIterator.cuh"

void MolDynIterator::UpdatePositions(ParticleSystem* particles)
{
	DeviceUpdatePositions << < particles.numBlocks,
		particles->numThreadsPerBlock >> > (particles, integrator,
		devicePos, deviceVel, deviceForce, deltaT, maxX, maxY, maxZ, kinEn);
}

void MolDynIterator::UpdateForces(unsigned long numParticles)
{
	DeviceUpdateForces << < particles.numBlocks,
		particles->numThreadsPerBlock >> > (particles->numParticles, deviceForce,
		devicePos, maxX, maxY, maxZ, forceEvaluator);
}

void MolDynIterator::UpdateVelocitiesT(unsigned long numParticles)
{
	UpdateVelocitiesT << < particles.numBlocks,
		particles->numThreadsPerBlock >> > (particles->numParticles, deviceForce,
		deltaT, deviceVel);
}

__global__ void DeviceUpdatePositions(ParticleSystem* particles, IIntegrationEvaluator* integrator, 
	double* devicePos, double* deviceVel, double* deviceForce, double deltaT, double maxX,
	double maxY, double maxZ, double &kinEn)
{
	int c = 0;

	for (int i = 0; i < 3 * particles->numParticles; i++) {

		if (c == 3) {
			c = 0;
		}

		double newPos = integrator->Evaluate(devicePos[i], deviceVel[i],
			deviceForce[i], deltaT);
		/*
		double tempXrX = particles->pos[i * 3 + 0] - prevPos[i * 3 + 0];
		double tempYrY = particles->pos[i * 3 + 1] - prevPos[i * 3 + 1];
		double tempZrZ = particles->pos[i * 3 + 2] - prevPos[i * 3 + 2];

		tempXrX = tempXrX - (maxX * (round(tempXrX / maxX)));
		tempYrY = tempYrY - (maxY * (round(tempYrY / maxY)));
		tempZrZ = tempZrZ - (maxZ * (round(tempZrZ / maxZ)));

		xrX = (2 * particles->pos[i * 3 + 0]) - (particles->
		pos[i * 3 + 0] - tempXrX) + (pow(deltaT, 2) *
		force[i * 3 + 0]);
		yrY = (2 * particles->pos[i * 3 + 1]) - (particles->
		pos[i * 3 + 1] - tempYrY) + (pow(deltaT, 2) *
		force[i * 3 + 1]);
		zrZ = (2 * particles->pos[i * 3 + 2]) - (particles->
		pos[i * 3 + 2] - tempZrZ) + (pow(deltaT, 2) *
		force[i * 3 + 2]);

		vel[i * 3 + 0] = (xrX - (particles->pos[i * 3 + 0] -
		tempXrX)) / (2 * deltaT);
		vel[i * 3 + 1] = (yrY - (particles->pos[i * 3 + 1] -
		tempYrY)) / (2 * deltaT);
		vel[i * 3 + 2] = (zrZ - (particles->pos[i * 3 + 2] -
		tempZrZ)) / (2 * deltaT);

		double tempXrX = xrX - prevPos[i * 3 + 0];
		double tempYrY = yrY - prevPos[i * 3 + 1];
		double tempZrZ = zrZ - prevPos[i * 3 + 2];

		vel[i * 3 + 0] = (tempXrX) / (2 * deltaT);
		vel[i * 3 + 1] = (tempYrY) / (2 * deltaT);
		vel[i * 3 + 2] = (tempZrZ) / (2 * deltaT);

		*/

		/*
		prevPos[i * 3 + 0] = particles->pos[i * 3 + 0];
		prevPos[i * 3 + 1] = particles->pos[i * 3 + 1];
		prevPos[i * 3 + 2] = particles->pos[i * 3 + 2];
		*/
		double boundaryWidth = 0;

		if (c == 0) {
			boundaryWidth = maxX;
		} else if (c == 1) {
			boundaryWidth = maxY;
		} else {
			boundaryWidth = maxZ;
			double velX = deviceVel[i - 2];
			double velY = deviceVel[i - 1];
			double velZ = deviceVel[i];
			kinEn += (velX * velX + velY * velY + velZ * velZ);
		}

		newPos = fmod(newPos, boundaryWidth);
		newPos = newPos < 0 ? boundaryWidth + newPos : newPos;

		devicePos[i] = newPos;
		++c;

		//comVel[c++] += vel[i];
	}
}

__global__ void DeviceUpdateForces(unsigned long numParticles, double* deviceForce,
	double* devicePos, double maxX, double maxY, double maxZ, IForceEvaluator* forceEvaluator)
{
	for (int i = 0; i < 3 * numParticles; i++) {
		deviceForce[i] = 0;
	}

	for (int i = 0; i < 3 * numParticles - 1; i += 3) {
		for (int j = i + 3; j < 3 * numParticles; j += 3) {
			double xDist = devicePos[j] - devicePos[i];
			double yDist = devicePos[j + 1] - devicePos[i + 1];
			double zDist = devicePos[j + 2] - devicePos[i + 2];

			xDist = xDist - (maxX * round(xDist / maxX));
			yDist = yDist - (maxY * round(yDist / maxY));
			zDist = zDist - (maxZ * round(zDist / maxZ));

			if (forceEvaluator->CheckCutoff(xDist, yDist, zDist)) {

				double scaledForce = forceEvaluator->
					EvaluateScaledForce();
				deviceForce[i] -= scaledForce * xDist;
				deviceForce[j] += scaledForce * xDist;
				deviceForce[i + 1] -= scaledForce * yDist;
				deviceForce[j + 1] += scaledForce * yDist;
				deviceForce[i + 2] -= scaledForce * zDist;
				deviceForce[j + 2] += scaledForce * zDist;

				forceEvaluator->EvaluateEnergy();
			}
		}
	}

	for (int i = 0; i < 3 * numParticles; i++) {
		deviceForce[i] = deviceForce[i] * 24;
	}
}

__global__ void DeviceUpdateVelocitiesT(unsigned long numParticles, double* deviceForce,
	double deltaT, double* deviceVel)
{
	for (int i = 0; i < 3 * numParticles; ++i) {
		double forceT = deviceForce[i] * deltaT * 0.5;
		deviceVel[i] += forceT;
	}
}