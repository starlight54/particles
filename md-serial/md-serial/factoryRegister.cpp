#include "factoryRegister.h"

//Forward type declarations
// 
//Iterators
//class MolDynIterator;

//Force Evaluators

//Integral Evaluators

void FactoryRegister::GlobalFactoryRegister()
{
	SimulationIteratorFactory::Get()->
		Register<MolDynIterator>("MolDynIterator");
}