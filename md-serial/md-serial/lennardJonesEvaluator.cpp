#include "lennardJonesEvaluator.h"

LennardJonesEvaluator::LennardJonesEvaluator()
{

}

LennardJonesEvaluator::~LennardJonesEvaluator()
{

}

void LennardJonesEvaluator::Initialise(double cutoff)
{
	totEn = 0;
	cutoffSquared = cutoff * cutoff;
	double cutoffSquaredInv = 1 / cutoffSquared;
	double cutoffSquaredInvCubed = cutoffSquaredInv * cutoffSquaredInv *
		cutoffSquaredInv;
	ecut = 4 * ((cutoffSquaredInvCubed * cutoffSquaredInvCubed) -
		cutoffSquaredInvCubed);
}

bool LennardJonesEvaluator::CheckCutoff(double xDist, double yDist, double zDist)
{
	rSquared = xDist * xDist + yDist * yDist + zDist * zDist;
	if (rSquared <= cutoffSquared) {
		return true;
	}

	return false;
}

void LennardJonesEvaluator::EvaluateEnergy()
{
	totEn += 4 * rSquaredInvCubed *
		(rSquaredInvCubed - 1) - ecut;
}

void LennardJonesEvaluator::EvaluateParticlePair(unsigned long i, unsigned long j)
{

}

double LennardJonesEvaluator::EvaluateScaledForce()
{
	double rSquaredInv = 1 / rSquared;
	rSquaredInvCubed = rSquaredInv * rSquaredInv * rSquaredInv;
	return rSquaredInv * rSquaredInvCubed * (2 * rSquaredInvCubed - 1);
}
