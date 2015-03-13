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
