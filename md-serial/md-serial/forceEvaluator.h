#ifndef FORCE_EVALUATOR_H
#define FORCE_EVALUATOR_H

class ForceEvaluator
{
public:
        ForceEvaluator();
        virtual ~ForceEvaluator();
        virtual void Evaluate() = 0;
};

#endif //FORCE_EVALUTAOR_H