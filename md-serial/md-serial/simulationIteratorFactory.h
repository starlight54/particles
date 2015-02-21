#ifndef SIMULATION_ITERATOR_FACTORY_H
#define SIMULATION_ITERATOR_FACTORY_H

#include "ISimulationIterator.h"
#include <string>
#include <map>
//Need to figure out how to remove the iterator headers from here, nasty!
#include "mdIterator.h"

class SimulationIteratorFactory
{
public:
	static SimulationIteratorFactory* Get();
        ~SimulationIteratorFactory();
	template <typename T>
        void Register(const char* iteratorTypeName);
	ISimulationIterator* Create(std::string iteratorTypeName);
private:
        //SimulationIteratorFactory();
	template <typename T>
        static ISimulationIterator*  CreateType();
        typedef ISimulationIterator* (*MapCreateType)();
	typedef std::map<std::string, MapCreateType> TypeMap;
        TypeMap typeMap;
};

template <typename T>
void SimulationIteratorFactory::Register(const char* typeName)
{
	typeMap[typeName] = &CreateType < T > ;
}

template <typename T>
ISimulationIterator* SimulationIteratorFactory::CreateType()
{
	return new T;
}

#endif //SIMULATION_ITERATOR_FACTORY_H