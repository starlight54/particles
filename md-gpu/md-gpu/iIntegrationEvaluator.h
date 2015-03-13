#ifndef INTEGRATION_EVALUTAOR_H
#define INTEGRATION_EVALUTAOR_H

class IIntegrationEvaluator
{
public:
        IIntegrationEvaluator();
        virtual ~IIntegrationEvaluator();
        virtual void Evaluate() = 0;
	virtual double Evaluate(double pos, double &vel, double force, 
		double deltaT) = 0;
};

#endif //INTEGRATION_EVALUTAOR_H