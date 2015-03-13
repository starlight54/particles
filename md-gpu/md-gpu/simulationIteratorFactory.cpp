/* */

#include "simulationIteratorFactory.h"

SimulationIteratorFactory* SimulationIteratorFactory::Get()
{
	static SimulationIteratorFactory iteratorFactory;
	return &iteratorFactory;
}

SimulationIteratorFactory::~SimulationIteratorFactory()
{
	typeMap.clear();
}

ISimulationIterator* SimulationIteratorFactory::Create(std::string iteratorTypeName)
{
	TypeMap::iterator kvp = typeMap.find(iteratorTypeName);

	if (kvp != typeMap.end()) {
		return kvp->second();
	}

	return NULL;
}