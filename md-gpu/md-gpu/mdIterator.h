#ifndef MD_ITERATOR_H
#define MD_ITERATOR_H

#define cudaCheckErrors(msg) \
    do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
                msg, cudaGetErrorString(__err), \
                __FILE__, __LINE__); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            exit(1); \
	        } \
        } while (0)


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
#include "device_launch_parameters.h"

class MolDynIterator : public ISimulationIterator
{
public:
        MolDynIterator();
        ~MolDynIterator();
	void Print(ParticleSystem* particles);
	void Iterate(ParticleSystem* particles);
	void Iterate(ParticleSystem* particles, int numBlocks = 0, int numThreadsPerBlock = 0);
	void Initialise(unsigned long numberIterations)
	{};
	void Initialise(ParticleSystem* particles, unsigned long numberIterations
		, double temperature, double deltaT, double cutoff, double maxX,
		double maxY, double maxZ);

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
	ParticleSystem* deviceParticles;
	MolDynIterator* deviceIterator;
	double* devKinEn;
	double* devTotEn;
	float flKinEn;
	float flTotEn;

	double maxX;
	double maxY;
	double maxZ;
	double cutoff;
	double xDist;
	double yDist;
	double zDist;
	double thickness;

//protected:
	void UpdatePositions(ParticleSystem* particles);
	void UpdateForces(ParticleSystem* particles);
	void UpdateVelocitiesT(unsigned long numParticles);
	/*
	__device__ void DeviceUpdateVelocitiesT(unsigned long numParticles, double* deviceForce,
		double deltaT, double* deviceVel, int pid);
	__device__ void DeviceUpdateForces(unsigned long numParticles, double* deviceForce,
		double* devicePos, double maxX, double maxY, double maxZ, double* devTotEn, double cutoffSquared, int pid);
	__device__ void DeviceUpdatePositions(unsigned long numParticles,
		double* devicePos, double* deviceVel, double* deviceForce, double deltaT, double maxX,
		double maxY, double maxZ, double* devKinEn, int pid);
		*/

private:
	/*
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
	*/
};

#if 0
__global__ inline void DeviceIterate(MolDynIterator* iterator, unsigned long numParticles, double* devicePos, double* deviceVel, double* deviceForce, double deltaT, double maxX,
	double maxY, double maxZ, double* devKinEn, double cutoff, double* devTotEn)
{
	//double r[3] = {0, 0, 0};
	//Move this out!

	//glutPostRedisplay();
	//glutMainLoopEvent();

	int pid = (blockIdx.x * blockDim.x) + threadIdx.x;
	devTotEn[0] = 0;
	devKinEn[0] = 0;

	iterator->DeviceUpdatePositions(numParticles,
		devicePos, deviceVel, deviceForce, deltaT, maxX, maxY, maxZ, devKinEn, pid);

	iterator->DeviceUpdateForces(numParticles, deviceForce,
		devicePos, maxX, maxY, maxZ, devTotEn, cutoff, pid);

	iterator->DeviceUpdateVelocitiesT(numParticles, deviceForce,
		deltaT, deviceVel, pid);

	//double totEn = devTotEn[0];

	//Now integrate!
	//comVel[0] = comVel[1] = comVel[2] = 0;
	//

	if (pid == numParticles - 1) {
		double instantTemp = devKinEn[0] / (3 * numParticles);
		double energyPerParticle = (devTotEn[0] + (0.5 * devKinEn[0])) / numParticles;

		printf("Instant temp: ");
		printf("%f", instantTemp);
		printf("\n");
		printf("Energy per particle: ");
		printf("%f", energyPerParticle);
		printf("\n\n");
	}

	//std::ofstream output("partEn250005r4.txt", std::ios_base::app);
	//output << iterator->energyPerParticle << "\n";

	/*
	printf("Com Vel: ");
	printf("%f", comVel[0]);
	printf("\n");
	printf("%f", comVel[1]);
	printf("\n");
	printf("%f", comVel[2]);
	printf("\n\n");
	*/
};
#endif

__global__ inline void DeviceUpdatePositions(unsigned long numParticles,
	double* devicePos, double* deviceVel, double* deviceForce, double deltaT, double maxX,
	double maxY, double maxZ, double* devKinEn)
{

	/*

	double newPos = integrator->Evaluate(devicePos[i], deviceVel[i],
	deviceForce[i], deltaT);

	*/
	int pid = ((blockIdx.x * blockDim.x) + threadIdx.x) * 3;
	if (pid >= numParticles * 3) {
		return;
	}

	int c = 0;

	for (int i = 0; i < 3; i++) {

		double forceT = deviceForce[pid + i] * deltaT * 0.5;
		deviceVel[pid + i] += forceT;
		double newPos = devicePos[pid + i] + deviceVel[pid + i] * deltaT + deltaT * forceT;

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
			double velX = deviceVel[pid + i - 2];
			double velY = deviceVel[pid + i - 1];
			double velZ = deviceVel[pid + i];
			devKinEn[pid / 3] = velX * velX + velY * velY + velZ * velZ; // devKinEn[0] += (velX * velX + velY * velY + velZ * velZ);
			//devKinEn[pid / 3 + 1] = 0;
			//devKinEn[pid / 3 + 2] = 0;
		}

		++c;

		newPos = fmod(newPos, boundaryWidth);
		newPos = newPos < 0 ? boundaryWidth + newPos : newPos;

		devicePos[pid + i] = newPos;
		//comVel[c++] += vel[i];

	}
};

__global__ inline void DeviceUpdateForces(unsigned long numParticles, double* deviceForce,
	double* devicePos, double maxX, double maxY, double maxZ, double* devTotEn, double cutoff)
{
	int pid = ((blockIdx.x * blockDim.x) + threadIdx.x) * 3;
	if (pid >= numParticles * 3) {
		return;
	}

	for (int i = 0; i < 3; i++) {
		deviceForce[pid + i] = 0;
	}

	for (int j = 0; j < 3 * numParticles; j += 3) {
		if (pid != j) {
			double xDist = devicePos[j] - devicePos[pid];
			double yDist = devicePos[j + 1] - devicePos[pid + 1];
			double zDist = devicePos[j + 2] - devicePos[pid + 2];

			xDist = xDist - (maxX * round(xDist / maxX));
			yDist = yDist - (maxY * round(yDist / maxY));
			zDist = zDist - (maxZ * round(zDist / maxZ));
			double rSquared = xDist * xDist + yDist * yDist + zDist * zDist;
			double cutoffSquared = cutoff * cutoff;
			/*
			if (forceEvaluator->CheckCutoff(xDist, yDist, zDist)) {
			*/
			if (rSquared <= cutoffSquared) {

				/*
				double scaledForce = forceEvaluator->
				EvaluateScaledForce();
				*/

				double rSquaredInv = 1 / rSquared;
				double rSquaredInvCubed = rSquaredInv * rSquaredInv * rSquaredInv;
				double scaledForce = rSquaredInv * rSquaredInvCubed * (2 * rSquaredInvCubed - 1);

				deviceForce[pid] -= scaledForce * xDist;
				//deviceForce[j] += scaledForce * xDist;
				deviceForce[pid + 1] -= scaledForce * yDist;
				//deviceForce[j + 1] += scaledForce * yDist;
				deviceForce[pid + 2] -= scaledForce * zDist;
				//deviceForce[j + 2] += scaledForce * zDist;

				//forceEvaluator->EvaluateEnergy();
				double cutoffSquaredInv = 1 / cutoffSquared;
				double cutoffSquaredInvCubed = cutoffSquaredInv * cutoffSquaredInv *
					cutoffSquaredInv;
				double ecut = 4 * ((cutoffSquaredInvCubed * cutoffSquaredInvCubed) -
					cutoffSquaredInvCubed);
				devTotEn[pid/3] = (4 * rSquaredInvCubed *
					(rSquaredInvCubed - 1)) - ecut;
				//devTotEn[pid / 3 + 1] = 0;
				//devTotEn[pid / 3 + 2] = 0;
			}


		}
	}

	for (int i = 0; i < 3; i++) {
		deviceForce[pid + i] = deviceForce[pid + i] * 24;
	}

};

__global__ inline void DeviceUpdateVelocitiesT(unsigned long numParticles, double* deviceForce,
	double deltaT, double* deviceVel)
{
	int pid = ((blockIdx.x * blockDim.x) + threadIdx.x) * 3;
	if (pid >= numParticles * 3) {
		return;
	}

	for (int i = 0; i < 3; i++) {
		double forceT = deviceForce[pid + i] * deltaT * 0.5;
		deviceVel[pid + i] += forceT;
	}
};

__global__ inline void DeviceUpdateEnvironmentForces(unsigned long numParticles, double* deviceForce,
	double* devicePos, double maxX, double maxY, double maxZ, double thickness, double cutoff, double* deviceVel)
{
	int pid = ((blockIdx.x * blockDim.x) + threadIdx.x) * 3;
	if (pid >= numParticles * 3) {
		return;
	}

		if (devicePos[pid] <= 60 + thickness && devicePos[pid] >= 40 - thickness) {
			if (devicePos[pid + 1] <= 70 && devicePos[pid + 1] >= 20) {
				if (devicePos[pid + 2] <= 60 + thickness && devicePos[pid + 2] >= 40 - thickness) {
					//In box! center is 50,50

					//check if hitting cylinder
					double xPosRel = devicePos[pid] - 50;
					double zPosRel = devicePos[pid + 2] - 50;
					double rDispSquared = xPosRel * xPosRel + zPosRel * zPosRel;

					//radial disp on circle
					if (rDispSquared > (10 * 10) && rDispSquared < ( (10 + (thickness / double(2))) * (10 + (thickness / double(2))) ) )
					{
						//inner circle
						if (devicePos[pid + 1] < 70 - cutoff && deviceVel[pid + 1] < 0) {
							//apply up force
							/*
							double yDist = (70 - cutoff) - devicePos[pid + 1];
							double rSquared = yDist * yDist;
							double rSquaredInv = 1 / rSquared;
							double rSquaredInvCubed = rSquaredInv * rSquaredInv * rSquaredInv;
							double scaledForce = rSquaredInv * rSquaredInvCubed * (2 * rSquaredInvCubed - 1);
							deviceForce[pid + 1] -= scaledForce * yDist;
							*/
						} else {
							//work out particle it's hitting!
							double magPos = sqrt(rDispSquared);

							double xWallPos = (xPosRel / magPos) * (10 + cutoff);
							double zWallPos = (zPosRel / magPos) * (10 + cutoff);

							double xDist = xWallPos - xPosRel;
							double zDist = zWallPos - zPosRel;

							double rSquared = xDist * xDist + zDist * zDist;

							double rSquaredInv = 1 / rSquared;
							double rSquaredInvCubed = rSquaredInv * rSquaredInv * rSquaredInv;
							double scaledForce = rSquaredInv * rSquaredInvCubed * (2 * rSquaredInvCubed - 1);

							deviceForce[pid] -= scaledForce * xDist;
							deviceForce[pid + 2] -= scaledForce * zDist;
						}
					} else if (rDispSquared < ((10 + thickness) * (10 + thickness)) && rDispSquared > ( (10 + (thickness / double(2))) * (10 + (thickness / double(2))) ) ) {
						//outer circle
						if (devicePos[pid + 1] < 70 - cutoff && deviceVel[pid + 1] < 0) {
							//apply upward force!
							/*
							double yDist = (70 - cutoff) - devicePos[pid + 1];
							double rSquared = yDist * yDist;
							double rSquaredInv = 1 / rSquared;
							double rSquaredInvCubed = rSquaredInv * rSquaredInv * rSquaredInv;
							double scaledForce = rSquaredInv * rSquaredInvCubed * (2 * rSquaredInvCubed - 1);
							deviceForce[pid + 1] -= scaledForce * yDist;
							*/
						} else {
							//work out particle it's hitting!
							double magPos = sqrt(rDispSquared);

							double xWallPos = (xPosRel / magPos) * (10 + cutoff);
							double zWallPos = (zPosRel / magPos) * (10 + cutoff);

							double xDist = xWallPos - xPosRel;
							double zDist = zWallPos - zPosRel;

							double rSquared = xDist * xDist + zDist * zDist;

							double rSquaredInv = 1 / rSquared;
							double rSquaredInvCubed = rSquaredInv * rSquaredInv * rSquaredInv;
							double scaledForce = rSquaredInv * rSquaredInvCubed * (2 * rSquaredInvCubed - 1);

							deviceForce[pid] -= scaledForce * xDist;
							deviceForce[pid + 2] -= scaledForce * zDist;
						}
					}
					
					//check if hitting top of cylinder!
				}
			}
		}
	
}

#endif //MD_ITERATOR_H