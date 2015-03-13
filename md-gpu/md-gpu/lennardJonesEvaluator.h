#ifndef LJ_EVALUATOR_H
#define LJ_EVALUATOR_H

#include "iForceEvaluator.h"

class LennardJonesEvaluator : public IForceEvaluator
{
public:
        LennardJonesEvaluator();
	~LennardJonesEvaluator();
	void Initialise(double cutoff);
	__device__ bool CheckCutoff(double xDist, double yDist, double zDist);
	__device__ void EvaluateEnergy();
	void Evaluate()
	{};
	void EvaluateParticlePair(unsigned long i, unsigned long j);
	__device__ double EvaluateScaledForce();
	void ShouldEvaluate();
private:
	double* pos;
	double* force;
	double ecut;
	double cutoffSquared;
	double rSquared;
	double rSquaredInvCubed;
};

#endif //LJ_EVALUATOR_H