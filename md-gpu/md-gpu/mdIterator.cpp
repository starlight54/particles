#include "mdIterator.h"

MolDynIterator::MolDynIterator()
{

}

MolDynIterator::~MolDynIterator()
{

}

/*
__global__
void MolDynIterator::UpdateForces(unsigned long numParticles)
{

}

__global__ 
void MolDynIterator::UpdatePositions(ParticleSystem* particles)
{

}

__global__ 
void MolDynIterator::UpdateVelocitiesT(unsigned long numParticles)
{

}
*/

void MolDynIterator::Print(ParticleSystem* particles)
{
	std::ofstream output("partVel.txt");
	for (int i = 0; i < particles->numParticles; i++) {
		output << sqrt(pow(vel[i * 3 + 0], 2) + pow(vel[i * 3 + 1], 2) +
			pow(vel[i * 3 + 2], 2)) << " \n"; // behaves like cout - cout is also a stream
	}
}

void MolDynIterator::Iterate(ParticleSystem* particles)
{
	//double r[3] = {0, 0, 0};
	//Move this out!

	glutPostRedisplay();
	glutMainLoopEvent();

	forceEvaluator->totEn = 0;
	kinEn = 0;

	UpdatePositions(particles);

	UpdateForces(particles->numParticles);

	UpdateVelocitiesT(particles->numParticles);
		
	double totEn = forceEvaluator->totEn;

	//Now integrate!
	//comVel[0] = comVel[1] = comVel[2] = 0;
	//

	instantTemp = kinEn / (3 * particles->numParticles);
	energyPerParticle = (totEn + (0.5 * kinEn)) / particles->numParticles;

	printf("Instant temp: ");
	printf("%f", instantTemp);
	printf("\n");
	printf("Energy per particle: ");
	printf("%f", energyPerParticle);
	printf("\n\n");

	/*
	printf("Com Vel: ");
	printf("%f", comVel[0]);
	printf("\n");
	printf("%f", comVel[1]);
	printf("\n");
	printf("%f", comVel[2]);
	printf("\n\n");
	*/

}

void MolDynIterator::Initialise(ParticleSystem* particles,
	unsigned long numberIterations, double temperature, double deltaT,
	double cutoff, double maxX, double maxY, double maxZ)
{
	forceEvaluator = ForceEvaluatorFactory::Get()->Create("LennardJones");
	forceEvaluator->Initialise(cutoff);
	integrator = IntegrationEvaluatorFactory::Get()->Create("VelocityVerlet");

	this->numberIterations = numberIterations;
	this->deltaT = deltaT;
	this->cutoff = cutoff;
	this->maxX = maxX;
	this->maxY = maxY;
	this->maxZ = maxZ;

	vel = (double*)malloc(sizeof(double) * particles->numParticles * 3);
	prevPos = (double*)malloc(sizeof(double) * particles->numParticles * 3);
	force = (double*)malloc(sizeof(double) * particles->numParticles * 3);

	Visualiser::Get()->SetData(particles, vel);

	comVel[0] = comVel[1] = comVel[2] = 0;
	for (int i = 0; i < 3 * particles->numParticles; i++) {
		force[i] = 0;
	}

	for (int i = 0; i < particles->numParticles; i++) {
		//Velocities between 0.5 and -0.5 in each dir
		double velX = ((double)rand() / RAND_MAX) - 0.5;
		double velY = ((double)rand() / RAND_MAX) - 0.5;
		double velZ = ((double)rand() / RAND_MAX) - 0.5;
		vel[i * 3 + 0] = velX;
		vel[i * 3 + 1] = velY;
		vel[i * 3 + 2] = velZ;

		//Better for accuracy and speed, although probably an even better method exists?
		comVel[0] += velX;
		comVel[1] += velY;
		comVel[2] += velZ;

		//Scalar so have to work it out now...
		//kinEn += pow(velX, 2) + pow(velY, 2) + pow(velZ, 2);
	}

	comVel[0] = comVel[0] / particles->numParticles;
	comVel[1] = comVel[1] / particles->numParticles;
	comVel[2] = comVel[2] / particles->numParticles;

	//kinEn = kinEn / particles->numberParticles;

	//comVel[0] = comVel[1] = comVel[2] = 0;
	double velScaleFactor = 20; // sqrt((3 * temperature) / kinEn);

	for (int i = 0; i < particles->numParticles; i++) {
		vel[i * 3 + 0] = (vel[i * 3 + 0] - comVel[0]) * velScaleFactor;
		vel[i * 3 + 1] = (vel[i * 3 + 1] - comVel[1]) * velScaleFactor;
		vel[i * 3 + 2] = (vel[i * 3 + 2] - comVel[2]) * velScaleFactor;

		/*
		double prevX = fmod((particles->pos[i * 3 + 0] + (vel[i * 3 + 0]
			* deltaT)), maxX);
		double prevY = fmod((particles->pos[i * 3 + 1] + (vel[i * 3 + 1]
			* deltaT)), maxY);
		double prevZ = fmod((particles->pos[i * 3 + 2] + (vel[i * 3 + 2]
			* deltaT)), maxZ);

		particles->pos[i * 3 + 0] = prevX < 0 ? maxX + prevX : prevX;
		particles->pos[i * 3 + 1] = prevY < 0 ? maxY + prevY : prevY;
		particles->pos[i * 3 + 2] = prevZ < 0 ? maxZ + prevZ : prevZ;
		*/
	}

	for (int i = 0; i < 3 * particles->numParticles - 1; i += 3) {
		for (int j = i + 3; j < 3 * particles->numParticles; j += 3) {
			xDist = particles->pos[j] - particles->pos[i];
			yDist = particles->pos[j + 1] - particles->pos[i + 1];
			zDist = particles->pos[j+ 2] - particles->pos[i + 2];

			xDist = xDist - (maxX * round(xDist / maxX));
			yDist = yDist - (maxY * round(yDist / maxY));
			zDist = zDist - (maxZ * round(zDist / maxZ));

			if (forceEvaluator->CheckCutoff(xDist, yDist, zDist)) {

				double scaledForce = forceEvaluator->
					EvaluateScaledForce();
				force[i] -= scaledForce * xDist;
				force[j] += scaledForce * xDist;
				force[i + 1] -= scaledForce * yDist;
				force[j + 1] += scaledForce * yDist;
				force[i + 2] -= scaledForce * zDist;
				force[j + 2] += scaledForce * zDist;

				forceEvaluator->EvaluateEnergy();
			}
		}
	}

	for (int i = 0; i < 3 * particles->numParticles; i++) {
		force[i] = force[i] * 24;
	}

	double totEn = forceEvaluator->totEn;

	//Copy to GPU
	cudaMalloc((void**)&devicePos, sizeof(double) * particles->numParticles * 3);
	cudaMalloc((void**)&deviceVel, sizeof(double) * particles->numParticles * 3);
	cudaMalloc((void**)&deviceForce, sizeof(double) * particles->numParticles * 3);

	cudaMemcpy(devicePos, particles->pos,
		sizeof(double) * particles->numParticles * 3, cudaMemcpyHostToDevice);
	cudaMemcpy(deviceVel, vel, sizeof(double) * particles->numParticles * 3, 
		cudaMemcpyHostToDevice);
	cudaMemcpy(deviceForce, force,
		sizeof(double) * particles->numParticles * 3, cudaMemcpyHostToDevice);
}

void MolDynIterator::CreateVelocities()
{

}

void MolDynIterator::CreateInitPrevPos()
{

}

void MolDynIterator::UpdatePositions(ParticleSystem* particles)
{
	this->DeviceUpdatePositions << < particles.numBlocks,
		particles->numThreadsPerBlock >> > (particles, integrator,
		devicePos, deviceVel, deviceForce, deltaT, maxX, maxY, maxZ, kinEn);
}

void MolDynIterator::UpdateForces(unsigned long numParticles)
{
	this->DeviceUpdateForces << < particles.numBlocks,
		particles->numThreadsPerBlock >> > (particles->numParticles, deviceForce,
		devicePos, maxX, maxY, maxZ, forceEvaluator);
}

void MolDynIterator::UpdateVelocitiesT(unsigned long numParticles)
{
	this->DeviceUpdateVelocitiesT << < particles.numBlocks,
		particles->numThreadsPerBlock >> > (particles->numParticles, deviceForce,
		deltaT, deviceVel);
}

__device__ void MolDynIterator::DeviceUpdatePositions(ParticleSystem* particles, IIntegrationEvaluator* integrator,
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

__device__ void MolDynIterator::DeviceUpdateForces(unsigned long numParticles, double* deviceForce,
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

__device__ void MolDynIterator::DeviceUpdateVelocitiesT(unsigned long numParticles, double* deviceForce,
	double deltaT, double* deviceVel)
{
	for (int i = 0; i < 3 * numParticles; ++i) {
		double forceT = deviceForce[i] * deltaT * 0.5;
		deviceVel[i] += forceT;
	}
}