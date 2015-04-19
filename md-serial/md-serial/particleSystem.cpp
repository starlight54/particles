#include "particleSystem.h"
#include "simulationSpace.h"

ParticleSystem::ParticleSystem()
{

}

ParticleSystem::~ParticleSystem()
{

}

void ParticleSystem::AutoInit(double maxX, double maxY, double maxZ, unsigned long numberParticles)
{
	srand(315496841515);
	this->numberParticles = numberParticles;
        pos = (double*)malloc(sizeof(double) * numberParticles * 3);

        for (int i = 0; i < numberParticles; ++i) {
                //No overlap protection, should use lattice in final approach? And a proper random number generator
		
		
		double posX = ((double)rand() / RAND_MAX) * maxX;
		double posY = ((double)rand() / RAND_MAX) * maxY;
		double posZ = ((double)rand() / RAND_MAX) * maxZ;
		
                pos[i * 3 + 0] = posX;
                pos[i * 3 + 1] = posY;
                pos[i * 3 + 2] = posZ;
		
		
		/*
		pos[i * 3 + 0] = (i + 1) * 1;
		pos[i * 3 + 1] = 50;
		pos[i * 3 + 2] = 50;
		*/
        }
}