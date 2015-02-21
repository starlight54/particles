#ifndef LJ_EVALUATOR_H
#define LJ_EVALUATOR_H

#include "forceEvaluator.h"

class LennardJonesEvaluator : public ForceEvaluator
{
public:
        LennardJonesEvaluator();
        ~LennardJonesEvaluator();
        void Evaluate();
private:

};

#endif //LJ_EVALUATOR_H