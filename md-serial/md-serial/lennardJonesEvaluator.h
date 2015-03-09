#ifndef LJ_EVALUATOR_H
#define LJ_EVALUATOR_H

#include "iForceEvaluator.h"

class LennardJonesEvaluator : public IForceEvaluator
{
public:
        LennardJonesEvaluator();
	~LennardJonesEvaluator();
	void Initialise(double cutoff);
	bool CheckCutoff(double xDist, double yDist, double zDist);
	void EvaluateEnergy();
	void Evaluate()
	{};
	void EvaluateParticlePair(unsigned long i, unsigned long j);
	double EvaluateScaledForce();
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