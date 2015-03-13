#ifndef FACTORY_REGISTER_H
#define FACTORY_REGISTER_H

#include "simulationIteratorFactory.h"
#include "forceEvaluatorFactory.h"
#include "integrationEvaluatorFactory.h"

class FactoryRegister
{
public:
	static void GlobalFactoryRegister();
};

#endif //FACTORY_REGISTER_H