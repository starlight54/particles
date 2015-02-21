#ifndef INTEGRATION_EVALUTAOR_H
#define INTEGRATION_EVALUTAOR_H

class IntegrationEvaluator
{
public:
        IntegrationEvaluator();
        ~IntegrationEvaluator();
        virtual void Evaluate() = 0;
};

#endif //INTEGRATION_EVALUTAOR_H