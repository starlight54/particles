#ifndef INTEGRATION_EVALUATOR_FACTORY_H
#define INTEGRATION_EVALUATOR_FACTORY_H

#include "iIntegrationEvaluator.h"
#include <string>
#include <map>

//Nasty
#include "lennardJonesEvaluator.h"

class IntegrationEvaluatorFactory
{
public:
	static IntegrationEvaluatorFactory* Get();
	~IntegrationEvaluatorFactory();
	IIntegrationEvaluator* Create(std::string evaluatorTypeName);
	template <typename T>
	void Register(const char* iteratorTypeName);
private:
	template <typename T>
	static IIntegrationEvaluator*  CreateType();
	typedef IIntegrationEvaluator* (*MapCreateType)();
	typedef std::map<std::string, MapCreateType> TypeMap;
	TypeMap typeMap;
};

template <typename T>
void IntegrationEvaluatorFactory::Register(const char* typeName)
{
	typeMap[typeName] = &CreateType < T > ;
}

template <typename T>
IIntegrationEvaluator* IntegrationEvaluatorFactory::CreateType()
{
	return new T;
}

#endif