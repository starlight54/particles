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
        pos = (double*)malloc(sizeof(double) * numberParticles * 3);

        for (int i = 0; i < numberParticles; i++) { 
                //No overlap protection, should use lattice in final approach? And a proper random number generator
                double posX = (double)rand() * maxX;
                double posY = (double)rand() * maxY;
                double posZ = (double)rand() * maxZ;
                pos[i * 3 + 0] = posX;
                pos[i * 3 + 1] = posY;
                pos[i * 3 + 2] = posZ;
        }
}