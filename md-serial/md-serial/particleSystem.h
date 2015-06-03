#ifndef PARTICLE_SYSTEM_H
#define PARTICLE_SYSTEM_H

#include <stdlib.h>
#include "dSFMT.h"

class ParticleSystem
{
public:
        ParticleSystem();
        ~ParticleSystem();
        //Need to set up some sort of proper inputs for this, for second section after MD on GPU
        void AutoInit(double maxX, double maxY, double maxZ, unsigned long numberParticles);
        void CreateParticleData(double maxX, double maxY, double mazZ);
        void CreatePositions();
        unsigned long numberParticles;
        double* pos;
	float* mass;
	dsfmt_t dsfmtRand;
private:
};

#endif //PARTICLE_SYSTEM_H