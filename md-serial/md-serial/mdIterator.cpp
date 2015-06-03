#include "mdIterator.h"

MolDynIterator::MolDynIterator()
{

}

MolDynIterator::~MolDynIterator()
{

}

void EvaluateParticlesAxisDist()
{

}

void MolDynIterator::Print(ParticleSystem* particles)
{
	std::ofstream output("partVel.txt");
	for (int i = 0; i < particles->numberParticles; i++) {
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
	int c = 0;

	for (int i = 0; i < 3 * particles->numberParticles; i++) {

		if (c == 3) {
			c = 0;
		}

		double newPos = integrator->Evaluate(particles->pos[i], vel[i],
			force[i], deltaT);

		double boundaryWidth = 0;

		if (c == 0) {
			boundaryWidth = maxX;
		} else if (c == 1) {
			boundaryWidth = maxY;
		} else {
			boundaryWidth = maxZ;
			double velX = vel[i - 2];
			double velY = vel[i - 1];
			double velZ = vel[i];
			kinEn += (velX * velX + velY * velY + velZ * velZ);
		}

		newPos = fmod(newPos, boundaryWidth);
		newPos = newPos < 0 ? boundaryWidth + newPos : newPos;

		particles->pos[i] = newPos;
		++c;

		//comVel[c++] += vel[i];
	}

	for (int i = 0; i < 3 * particles->numberParticles; i++) {
		force[i] = 0;
	}

	for (int i = 0; i < 3 * particles->numberParticles - 1; i += 3) {
		for (int j = i + 3; j < 3 * particles->numberParticles; j += 3) {
			xDist = particles->pos[j] - particles->pos[i];
			yDist = particles->pos[j + 1] - particles->pos[i + 1];
			zDist = particles->pos[j + 2] - particles->pos[i + 2];

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

	for (int i = 0; i < 3 * particles->numberParticles; i++) {
		force[i] = force[i] * 24;
	}

	for (int i = 0; i < 3 * particles->numberParticles; ++i) {
		double forceT = force[i] * deltaT * 0.5;
		vel[i] += forceT;
	}

	double totEn = forceEvaluator->totEn;

	//Now integrate!
	//comVel[0] = comVel[1] = comVel[2] = 0;

	instantTemp = kinEn / (3 * particles->numberParticles);
	energyPerParticle = (totEn + (0.5 * kinEn)) / particles->numberParticles;

	printf("Instant temp: ");
	printf("%f", instantTemp);
	printf("\n");
	printf("Energy per particle: ");
	printf("%f", energyPerParticle);
	printf("\n\n");

	//std::ofstream output("partEn250005r4.txt", std::ios_base::app);
	//output << energyPerParticle << "\n";

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

	vel = (double*)malloc(sizeof(double) * particles->numberParticles * 3);
	prevPos = (double*)malloc(sizeof(double) * particles->numberParticles * 3);
	force = (double*)malloc(sizeof(double) * particles->numberParticles * 3);

	Visualiser::Get()->SetData(particles, vel);

	comVel[0] = comVel[1] = comVel[2] = 0;
	for (int i = 0; i < 3 * particles->numberParticles; i++) {
		force[i] = 0;
	}

	for (int i = 0; i < particles->numberParticles; i++) {
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

	comVel[0] = comVel[0] / particles->numberParticles;
	comVel[1] = comVel[1] / particles->numberParticles;
	comVel[2] = comVel[2] / particles->numberParticles;

	//kinEn = kinEn / particles->numberParticles;

	//comVel[0] = comVel[1] = comVel[2] = 0;
	double velScaleFactor = 10; // sqrt((3 * temperature) / kinEn);

	for (int i = 0; i < particles->numberParticles; i++) {
		vel[i * 3 + 0] = (vel[i * 3 + 0] - comVel[0]) * velScaleFactor;
		vel[i * 3 + 1] = (vel[i * 3 + 1] - comVel[1]) * velScaleFactor;
		vel[i * 3 + 2] = (vel[i * 3 + 2] - comVel[2]) * velScaleFactor;
	}

	for (int i = 0; i < 3 * particles->numberParticles - 1; i += 3) {
		for (int j = i + 3; j < 3 * particles->numberParticles; j += 3) {
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

	for (int i = 0; i < 3 * particles->numberParticles; i++) {
		force[i] = force[i] * 24;
	}

	double totEn = forceEvaluator->totEn;
}

void MolDynIterator::CreateVelocities()
{

}

void MolDynIterator::CreateInitPrevPos()
{

}
