#ifndef INTEGRATION_EVALUTAOR_H
#define INTEGRATION_EVALUTAOR_H

#include <cuda_runtime_api.h> //Particularly nasty, solely for IDE intellisense!

class IIntegrationEvaluator
{
public:
        IIntegrationEvaluator();
        virtual ~IIntegrationEvaluator();
        virtual void Evaluate() = 0;
	__device__ virtual double Evaluate(double pos, double &vel, double force,
		double deltaT) = 0;
};

#endif //INTEGRATION_EVALUTAOR_H