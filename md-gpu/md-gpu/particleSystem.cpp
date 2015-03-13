#include "particleSystem.h"

ParticleSystem::ParticleSystem()
{

}

ParticleSystem::~ParticleSystem()
{

}

void ParticleSystem::AutoInit(double maxX, double maxY, double maxZ, unsigned long numParticles)
{
	srand(315496841515);
	this->numParticles = numParticles;
        pos = (double*)malloc(sizeof(double) * numParticles * 3);

	//GPU Init
	numThreadsPerBlock = ConfigManager::Get()->device.maxThreadsPerBlock;
	numBlocks = ceil(numParticles / numThreadsPerBlock);

        for (int i = 0; i < numParticles; ++i) {
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