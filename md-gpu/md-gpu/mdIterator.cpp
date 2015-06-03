#include "mdIterator.h"

MolDynIterator::MolDynIterator()
{

}

MolDynIterator::~MolDynIterator()
{

}

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

}

void MolDynIterator::Iterate(ParticleSystem* particles, int numBlocks, int numThreadsPerBlock)
{
	//system("pause");
	//printf("%i", numThreadsPerBlock);
	//printf("%i", numBlocks);
	//DeviceIterate << <numBlocks, numThreadsPerBlock>> >(particles, this);

	//DeviceIterate << <numBlocks, numThreadsPerBlock >> >(deviceIterator, particles->numParticles, devicePos, deviceVel, deviceForce, deltaT, maxX,
	//	maxY, maxZ, devKinEn, cutoff, devTotEn);
	glutPostRedisplay();
	glutMainLoopEvent();

	DeviceUpdatePositions << <numBlocks, numThreadsPerBlock >> >(particles->numParticles,
		devicePos, deviceVel, deviceForce, deltaT, maxX, maxY, maxZ, devKinEn);

	cudaDeviceSynchronize();

	DeviceUpdateForces << <numBlocks, numThreadsPerBlock >> >(particles->numParticles, deviceForce,
		devicePos, maxX, maxY, maxZ, devTotEn, cutoff);

	cudaDeviceSynchronize();

	thickness = 2 * cutoff;

	//DeviceUpdateEnvironmentForces << <numBlocks, numThreadsPerBlock >> >(particles->numParticles, deviceForce, devicePos, maxX, maxY, maxZ, thickness, cutoff, deviceVel);

	//cudaDeviceSynchronize();

	DeviceUpdateVelocitiesT << <numBlocks, numThreadsPerBlock >> >(particles->numParticles, deviceForce,
		deltaT, deviceVel);

	cudaDeviceSynchronize();



	/*REMOVE
	//double r[3] = {0, 0, 0};
	//Move this out!

	glutPostRedisplay();
	glutMainLoopEvent();

	forceEvaluator->totEn = 0;
	kinEn = 0;

	DeviceUpdatePositions(particles, integrator,
		particles->pos, vel, force, deltaT, maxX, maxY, maxZ, kinEn);

	DeviceUpdateForces(particles->numParticles, force,
		particles->pos, maxX, maxY, maxZ, forceEvaluator);

	DeviceUpdateVelocitiesT(particles->numParticles, force,
		deltaT, vel);
		
	double totEn = forceEvaluator->totEn;

	//Now integrate!
	//comVel[0] = comVel[1] = comVel[2] = 0;
	//

	/*
	printf("Com Vel: ");
	printf("%f", comVel[0]);
	printf("\n");
	printf("%f", comVel[1]);
	printf("\n");
	printf("%f", comVel[2]);
	printf("\n\n");
	*/

	cudaMemcpy(kinEn, devKinEn, sizeof(double) * particles->numParticles, cudaMemcpyDeviceToHost);
	cudaMemcpy(totEn, devTotEn, sizeof(double) * particles->numParticles, cudaMemcpyDeviceToHost);
	//cudaMemcpy(particles->pos, devicePos, sizeof(double) * particles->numParticles * 3, cudaMemcpyDeviceToHost);
	//cudaMemcpy(vel, deviceVel, sizeof(double) * particles->numParticles * 3, cudaMemcpyDeviceToHost);

	double sumKin = 0;
	double sumTot = 0;

	for (int i = 0; i < particles->numParticles; i++) {
		sumKin += kinEn[i];
		sumTot += totEn[i];
	}

	instantTemp = sumKin / (3 * particles->numParticles);
	energyPerParticle = (sumTot + (0.5 * sumKin)) / particles->numParticles;

	printf("Instant temp: ");
	printf("%f", instantTemp);
	printf("\n");
	printf("Energy per particle: ");
	printf("%f", energyPerParticle);
	printf("\n\n");
	
	std::ofstream output("partEn250005r4.txt", std::ios_base::app);
	output << energyPerParticle << "\n";
}

void MolDynIterator::Initialise(ParticleSystem* particles,
	unsigned long numberIterations, double temperature, double deltaT,
	double cutoff, double maxX, double maxY, double maxZ)
{
	//system("pause");
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
	kinEn = (double*)malloc(sizeof(double) * particles->numParticles);
	totEn = (double*)malloc(sizeof(double) * particles->numParticles);

	Visualiser::Get()->SetData(particles, vel);

	comVel[0] = comVel[1] = comVel[2] = 0;
	for (int i = 0; i < 3 * particles->numParticles; i++) {
		force[i] = 0;
	}

	for (int i = 0; i < particles->numParticles; i++) {
		//Velocities between 0.5 and -0.5 in each dir
		kinEn[i] = 0;
		totEn[i] = 0;
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
	double velScaleFactor = 10; // sqrt((3 * temperature) / kinEn);

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

	cudaMalloc((void**)&devicePos, sizeof(double) * particles->numParticles * 3);
	cudaMalloc((void**)&deviceVel, sizeof(double) * particles->numParticles * 3);
	cudaMalloc((void**)&deviceForce, sizeof(double) * particles->numParticles * 3);
	cudaMalloc((void**)&devKinEn, sizeof(double) * particles->numParticles);
	cudaMalloc((void**)&devTotEn, sizeof(double) * particles->numParticles);

	cudaMemcpy(devicePos, particles->pos, sizeof(double) * particles->numParticles * 3, cudaMemcpyHostToDevice);
	cudaMemcpy(deviceVel, vel, sizeof(double) * particles->numParticles * 3, cudaMemcpyHostToDevice);
	cudaMemcpy(deviceForce, force, sizeof(double) * particles->numParticles * 3, cudaMemcpyHostToDevice);
	cudaMemcpy(devKinEn, kinEn, sizeof(double) * particles->numParticles, cudaMemcpyHostToDevice);
	cudaMemcpy(devTotEn, totEn, sizeof(double) * particles->numParticles, cudaMemcpyHostToDevice);

	DeviceUpdateForces << <particles->numBlocks, particles->numThreadsPerBlock >> >(particles->numParticles, deviceForce,
		devicePos, maxX, maxY, maxZ, devTotEn, cutoff);

	cudaDeviceSynchronize();

	/*
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

	*/

	//double totEn = forceEvaluator->totEn;

	//Copy to GPU
	/*
	cudaMalloc((void**)&devicePos, sizeof(double) * particles->numParticles * 3);
	cudaMalloc((void**)&deviceVel, sizeof(double) * particles->numParticles * 3);
	cudaMalloc((void**)&deviceForce, sizeof(double) * particles->numParticles * 3);
	cudaMalloc((void**)&devKinEn, sizeof(double) * particles->numParticles);
	cudaMalloc((void**)&devTotEn, sizeof(double) * particles->numParticles);

	cudaMemcpy(devicePos, particles->pos, sizeof(double) * particles->numParticles * 3, cudaMemcpyHostToDevice);
	cudaMemcpy(deviceVel, vel, sizeof(double) * particles->numParticles * 3, cudaMemcpyHostToDevice);
	cudaMemcpy(deviceForce, force, sizeof(double) * particles->numParticles * 3, cudaMemcpyHostToDevice);
	cudaMemcpy(devKinEn, kinEn, sizeof(double) * particles->numParticles, cudaMemcpyHostToDevice);
	cudaMemcpy(devTotEn, totEn, sizeof(double) * particles->numParticles, cudaMemcpyHostToDevice);
	*/

	/*
	cudaMalloc((void**)&deviceParticles, sizeof(ParticleSystem));
	cudaCheckErrors("allocating first part pointer");
	cudaMalloc((void**)&deviceIterator, sizeof(MolDynIterator));
	cudaCheckErrors("allocating first iterator pointer");

	cudaMemcpy(deviceParticles, particles, sizeof(ParticleSystem), cudaMemcpyHostToDevice);
	cudaCheckErrors("copying particles fail");
	cudaMemcpy(deviceIterator, this, sizeof(MolDynIterator), cudaMemcpyHostToDevice);
	cudaCheckErrors("copying moldyn fail");
	system("pause");
	cudaMemcpy(deviceParticles->pos, particles->pos,
		sizeof(double) * particles->numParticles * 3, cudaMemcpyHostToDevice);
	cudaCheckErrors("device pos, part pos transfer fail");
	system("pause");
	cudaMemcpy(deviceIterator->vel, vel, sizeof(double) * particles->numParticles * 3,
		cudaMemcpyHostToDevice);
	system("pause");
	cudaMemcpy(deviceIterator->force, force,
		sizeof(double) * particles->numParticles * 3, cudaMemcpyHostToDevice);
	system("pause");
	*/

	//cudaMemcpy(deviceParticles->pos, devicePos, sizeof(devicePos), cudaMemcpyDeviceToDevice);
	//cudaMemcpy(deviceIterator->vel, deviceVel, sizeof(deviceVel), cudaMemcpyDeviceToDevice);
	//cudaMemcpy(deviceIterator->force, deviceForce, sizeof(deviceForce), cudaMemcpyDeviceToDevice);
	//system("pause");

}

void MolDynIterator::CreateVelocities()
{

}

void MolDynIterator::CreateInitPrevPos()
{

}

void MolDynIterator::UpdatePositions(ParticleSystem* particles)
{

}

void MolDynIterator::UpdateForces(ParticleSystem* particles)
{

}

void MolDynIterator::UpdateVelocitiesT(unsigned long numParticles)
{

}

#if 0

__device__ void MolDynIterator::DeviceUpdatePositions(unsigned long numParticles,
	double* devicePos, double* deviceVel, double* deviceForce, double deltaT, double maxX,
	double maxY, double maxZ, double* devKinEn, int pid)
{

		/*

		double newPos = integrator->Evaluate(devicePos[i], deviceVel[i],
			deviceForce[i], deltaT);

			*/
	pid = pid * 3;
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
				devKinEn[0] += (velX * velX + velY * velY + velZ * velZ);
			}

			++c;

			newPos = fmod(newPos, boundaryWidth);
			newPos = newPos < 0 ? boundaryWidth + newPos : newPos;

			devicePos[pid + i] = newPos;
			//comVel[c++] += vel[i];
		
	}
}

__device__ void MolDynIterator::DeviceUpdateForces(unsigned long numParticles, double* deviceForce,
	double* devicePos, double maxX, double maxY, double maxZ, double* devTotEn, double cutoff, int pid)
{
	pid = pid * 3;
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
					devTotEn[0] += 4 * rSquaredInvCubed *
						(rSquaredInvCubed - 1) - ecut;
				}


			}
		}

		for (int i = 0; i < 3; i++) {
			deviceForce[pid + i] = deviceForce[pid + i] * 24;
		}
	
}

__device__ void MolDynIterator::DeviceUpdateVelocitiesT(unsigned long numParticles, double* deviceForce,
	double deltaT, double* deviceVel, int pid)
{
	pid = pid * 3;
	if (pid >= numParticles * 3) {
		return;
	}

	for (int i = 0; i < 3; i++) {
		double forceT = deviceForce[pid + i] * deltaT * 0.5;
		deviceVel[pid + i] += forceT;
	}
}

*/

#endif