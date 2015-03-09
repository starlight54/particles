#ifndef FORCE_EVALUATOR_FACTORY_H
#define FORCE_EVALUATOR_FACTORY_H

#include "iForceEvaluator.h"
#include <string>
#include <map>

//Nasty
#include "velocityVerletIntegrator.h"

class ForceEvaluatorFactory
{
public:
	static ForceEvaluatorFactory* Get();
	~ForceEvaluatorFactory();
	IForceEvaluator* Create(std::string evaluatorTypeName);
	template <typename T>
	void Register(const char* iteratorTypeName);
private:
	template <typename T>
	static IForceEvaluator*  CreateType();
	typedef IForceEvaluator* (*MapCreateType)();
	typedef std::map<std::string, MapCreateType> TypeMap;
	TypeMap typeMap;
};

template <typename T>
void ForceEvaluatorFactory::Register(const char* typeName)
{
	typeMap[typeName] = &CreateType < T > ;
}

template <typename T>
IForceEvaluator* ForceEvaluatorFactory::CreateType()
{
	return new T;
}

#endif