#include "mdIterator.h"

MolDynIterator::MolDynIterator()
{

}

MolDynIterator::~MolDynIterator()
{

}

void MolDynIterator::Iterate(ParticleSystem* particles)
{
	//double r[3] = {0, 0, 0};
	double xrX;
	double yrY;
	double zrZ;
	double rSquared;
	//Move this out!
	double cutoffSquaredInvCubed = pow((1 / cutoffSquared), 3);
	double ecut = 4 * (pow(cutoffSquaredInvCubed, 2) - cutoffSquaredInvCubed);

	for (int i = 0; i < numberIterations; i++) {

		glutPostRedisplay();
		glutMainLoopEvent();

		scaledTotEn = 0;
		totEn = 0;

		for (int i = 0; i < particles->numberParticles; i++) {
			force[i * 3 + 0] = 0;
			force[i * 3 + 1] = 0;
			force[i * 3 + 2] = 0;
		}

		for (int i = 0; i < particles->numberParticles - 1; i++) {
			for (int j = i + 1; j < particles->numberParticles; j++) {
				xrX = (particles->pos[i * 3 + 0] -
					particles->pos[j * 3 + 0]);
				yrY = (particles->pos[i * 3 + 1] -
					particles->pos[j * 3 + 1]);
				zrZ = (particles->pos[i * 3 + 2] -
					particles->pos[j * 3 + 2]);
				xrX = xrX - (maxX * (round(xrX / maxX)));
				yrY = yrY - (maxY * (round(yrY / maxY)));
				zrZ = zrZ - (maxZ * (round(zrZ / maxZ)));

				//Too much nesting! Needs busting up....
				rSquared = pow(xrX, 2) + pow(yrY, 2) + pow(zrZ, 2);
				if (rSquared < cutoffSquared) {
					double rSquaredInv = 1 / rSquared;
					double rSquaredInvCubed = pow(rSquaredInv, 3);
					double scaledScalarForce = rSquaredInv *
						rSquaredInvCubed * (rSquaredInvCubed - 0.5);
					force[i * 3 + 0] = force[i * 3 + 0] +
						(scaledScalarForce * xrX);
					force[j * 3 + 0] = force[j * 3 + 0] -
						(scaledScalarForce * xrX);
					force[i * 3 + 1] = force[i * 3 + 1] +
						(scaledScalarForce * yrY);
					force[j * 3 + 1] = force[j * 3 + 1] -
						(scaledScalarForce * yrY);
					force[i * 3 + 2] = force[i * 3 + 2] +
						(scaledScalarForce * zrZ);
					force[j * 3 + 2] = force[j * 3 + 2] -
						(scaledScalarForce * zrZ);
					scaledTotEn += rSquaredInvCubed *
						(rSquaredInvCubed - 1) - ecut;
				}
			}
		}

		for (int i = 0; i < particles->numberParticles; i++) {
			force[i * 3 + 0] = force[i * 3 + 0] * 48;
			force[i * 3 + 1] = force[i * 3 + 1] * 48;
			force[i * 3 + 2] = force[i * 3 + 2] * 48;
		}

		totEn = scaledTotEn * 4;

		//Now integrate!
		comVel[0] = comVel[1] = comVel[2] = 0;
		kinEn = 0;

		for (int i = 0; i < particles->numberParticles; i++) {
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

				vel[i * 3 + 0] += force[i * 3 + 0] * deltaT;
				vel[i * 3 + 1] += force[i * 3 + 1] * deltaT;
				vel[i * 3 + 2] += force[i * 3 + 2] * deltaT;
			*/

			xrX = particles->pos[i * 3 + 0] + (vel[i * 3 + 0] *
				deltaT) + ((force[i * 3 + 0] / 2) * 
				pow(deltaT, 2));

			yrY = particles->pos[i * 3 + 1] + (vel[i * 3 + 1] *
				deltaT) + ((force[i * 3 + 1] / 2) *
				pow(deltaT, 2));

			zrZ = particles->pos[i * 3 + 2] + (vel[i * 3 + 2] *
				deltaT) + ((force[i * 3 + 2] / 2) *
				pow(deltaT, 2));

			double tempXrX = xrX - prevPos[i * 3 + 0];
			double tempYrY = yrY - prevPos[i * 3 + 1];
			double tempZrZ = zrZ - prevPos[i * 3 + 2];

			tempXrX = tempXrX - (maxX * (round(tempXrX / maxX)));
			tempYrY = tempYrY - (maxY * (round(tempYrY / maxY)));
			tempZrZ = tempZrZ - (maxZ * (round(tempZrZ / maxZ)));

			vel[i * 3 + 0] = (tempXrX) / (2 * deltaT);
			vel[i * 3 + 1] = (tempYrY) / (2 * deltaT);
			vel[i * 3 + 2] = (tempZrZ) / (2 * deltaT);

			comVel[0] += vel[i * 3 + 0];
			comVel[1] += vel[i * 3 + 1];
			comVel[2] += vel[i * 3 + 2];

			kinEn += pow(vel[i * 3 + 0], 2) + pow(vel[i * 3 + 1], 2) +
				pow(vel[i * 3 + 2], 2);

			prevPos[i * 3 + 0] = particles->pos[i * 3 + 0];
			prevPos[i * 3 + 1] = particles->pos[i * 3 + 1];
			prevPos[i * 3 + 2] = particles->pos[i * 3 + 2];

			
			xrX = fmod(xrX, maxX);
			yrY = fmod(yrY, maxY);
			zrZ = fmod(zrZ, maxZ);

			xrX = xrX < 0 ? maxX + xrX : xrX;
			yrY = yrY < 0 ? maxY + yrY : yrY;
			zrZ = zrZ < 0 ? maxZ + zrZ : zrZ;
			
			particles->pos[i * 3 + 0] = xrX;
			particles->pos[i * 3 + 1] = yrY;
			particles->pos[i * 3 + 2] = zrZ;
		}

		instantTemp = kinEn / (3 * particles->numberParticles);
		energyPerParticle = (totEn + (0.5 * kinEn)) / particles->numberParticles;

		printf("Instant temp: ");
		printf("%f", instantTemp);
		printf("\n");
		printf("Energy per particle: ");
		printf("%f", energyPerParticle);
		printf("\n\n");
	}

	std::ofstream output("partVel.txt");
	for (int i = 0; i < particles->numberParticles; i++) {
		output << sqrt(pow(vel[i * 3 + 0], 2) + pow(vel[i * 3 + 1], 2) + 
			pow(vel[i * 3 + 2], 2)) << " \n"; // behaves like cout - cout is also a stream
	}
}

void MolDynIterator::Initialise(ParticleSystem* particles, 
	unsigned long numberIterations, double temperature, double deltaT, 
	double cutoff, double maxX, double maxY, double maxZ)
{
	this->numberIterations = numberIterations;
	this->deltaT = deltaT;
	vel = (double*)malloc(sizeof(double) * particles->numberParticles * 3);
	prevPos = (double*)malloc(sizeof(double) * particles->numberParticles * 3);
	force = (double*)malloc(sizeof(double) * particles->numberParticles * 3);
	Visualiser::Get()->SetData(particles, vel);
	cutoffSquared = pow(cutoff, 2);
	this->maxX = maxX;
	this->maxY = maxY;
	this->maxZ = maxZ;
	comVel[0] = comVel[1] = comVel[2] = 0;
	for (int i = 0; i < particles->numberParticles; i++) {
		//Velocities between 0.5 and -0.5 in each dir
		double velX = ((double)rand() / RAND_MAX) -0.5;
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
		kinEn += pow(velX, 2) + pow(velY, 2) + pow(velZ, 2);
	}

	comVel[0] = comVel[0] / particles->numberParticles;
	comVel[1] = comVel[1] / particles->numberParticles;
	comVel[2] = comVel[2] / particles->numberParticles;

	kinEn = kinEn / particles->numberParticles;

	double velScaleFactor = sqrt((3 * temperature) / kinEn);

	for (int i = 0; i < particles->numberParticles; i++) {
		vel[i * 3 + 0] = (vel[i * 3 + 0] - comVel[0]) * velScaleFactor;
		vel[i * 3 + 1] = (vel[i * 3 + 1] - comVel[1]) * velScaleFactor;
		vel[i * 3 + 2] = (vel[i * 3 + 2] - comVel[2]) * velScaleFactor;

		double prevX = fmod((particles->pos[i * 3 + 0] - (vel[i * 3 + 0]
			* deltaT)), maxX);
		double prevY = fmod((particles->pos[i * 3 + 1] - (vel[i * 3 + 1]
			* deltaT)), maxY);
		double prevZ = fmod((particles->pos[i * 3 + 2] - (vel[i * 3 + 2]
			* deltaT)), maxZ);

		prevPos[i * 3 + 0] = prevX < 0 ? maxX + prevX : prevX;
		prevPos[i * 3 + 1] = prevY < 0 ? maxY + prevY : prevY;
		prevPos[i * 3 + 2] = prevZ < 0 ? maxZ + prevZ : prevZ;
	}

}

void MolDynIterator::CreateVelocities()
{

}

void MolDynIterator::CreateInitPrevPos()
{

}
