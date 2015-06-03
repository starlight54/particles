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
	numThreadsPerBlock = 256; // ConfigManager::Get()->device.maxThreadsPerBlock;
	numBlocks = ceil(double(numParticles) / double(numThreadsPerBlock));

	
	double axisSpace = std::cbrt(maxX * maxY * maxZ / double(numParticles));

	int c = 0;
	for (double i = axisSpace / 2; i < maxX; i += axisSpace) {
		for (double j = axisSpace / 2; j < maxY; j += axisSpace) {
			for (double k = axisSpace / 2; k < maxZ; k += axisSpace) {
				if (c < numParticles * 3 - 1) {
					pos[c] = i;
					pos[c + 1] = j;
					pos[c + 2] = k;
					c += 3;
				}
			}
		}
	}
	
	/*
        for (int i = 0; i < numParticles; ++i) {
                //No overlap protection, should use lattice in final approach? And a proper random number generator
		
		
		double posX = ((double)rand() / RAND_MAX) * maxX;
		double posY = ((double)rand() / RAND_MAX) * maxY;
		double posZ = ((double)rand() / RAND_MAX) * maxZ;

		//double posX = ((double)rand() / RAND_MAX) * 20 + 40;
		//double posY = ((double)rand() / RAND_MAX) * maxY;
		//double posZ = ((double)rand() / RAND_MAX) * 20 + 40;
		
                pos[i * 3 + 0] = posX;
                pos[i * 3 + 1] = posY;
                pos[i * 3 + 2] = posZ;
		
		
		//pos[i * 3 + 0] = (i + 1) * 1;
		//pos[i * 3 + 1] = 50;
		//pos[i * 3 + 2] = 50;
		
        }
	*/
	
}