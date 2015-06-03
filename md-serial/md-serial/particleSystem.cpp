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
	dsfmt_init_gen_rand(&dsfmtRand, 427285630591);
	this->numberParticles = numberParticles;
        pos = (double*)malloc(sizeof(double) * numberParticles * 3);
	double rand;

	double axisSpace = std::cbrt(maxX * maxY * maxZ / double(numberParticles));
	//double axisSapce = maxX / density;

	int c = 0;
	for (double i = 0; i < maxX; i += axisSpace) {
		for (double j = 0; j < maxY; j += axisSpace) {
			for (double k = 0; k < maxZ; k += axisSpace) {
				if (c < numberParticles * 3 - 1) {
					pos[c] = i;
					pos[c + 1] = j;
					pos[c + 2] = k;
					c += 3;
				}
			}
		}
	}

	/*
        for (int i = 0; i < numberParticles; ++i) {
                //No overlap protection, should use lattice in final approach? And a proper random number generator

		rand = dsfmt_genrand_close_open(&dsfmtRand);
		double posX = rand * maxX;
		rand = dsfmt_genrand_close_open(&dsfmtRand);
		double posY = rand * maxY;
		rand = dsfmt_genrand_close_open(&dsfmtRand);
		double posZ = rand * maxZ;
		
                pos[i * 3 + 0] = posX;
                pos[i * 3 + 1] = posY;
                pos[i * 3 + 2] = posZ;
		
		
		//pos[i * 3 + 0] = (i + 1) * 1;
		//pos[i * 3 + 1] = 50;
		//pos[i * 3 + 2] = 50;
		
        }
	*/
}