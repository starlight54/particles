#ifndef FORCE_EVALUATOR_H
#define FORCE_EVALUATOR_H

#include <host_defines.h>

class IForceEvaluator
{
public:
        IForceEvaluator();
        virtual ~IForceEvaluator();
	virtual void Initialise(double cutoff) = 0;
        virtual void Evaluate() = 0;
	__host__ __device__ virtual void EvaluateEnergy() = 0;
	__host__ __device__ virtual double EvaluateScaledForce() = 0;
	__host__ __device__ virtual bool CheckCutoff(double xDist, double yDist, double zDist) = 0;
	double totEn;
};

#endif //FORCE_EVALUTAOR_H