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
	

	//For text output
	//std::ofstream output("partEn.txt", std::ios_base::app);
	//output << energyPerParticle << "\n";
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
	double velScaleFactor = 1;  // sqrt((3 * temperature) / kinEn);

	for (int i = 0; i < particles->numParticles; i++) {
		vel[i * 3 + 0] = (vel[i * 3 + 0] - comVel[0]) * velScaleFactor;
		vel[i * 3 + 1] = (vel[i * 3 + 1] - comVel[1]) * velScaleFactor;
		vel[i * 3 + 2] = (vel[i * 3 + 2] - comVel[2]) * velScaleFactor;
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

}
