#include "forceEvaluatorFactory.h"

ForceEvaluatorFactory* ForceEvaluatorFactory::Get()
{
	static ForceEvaluatorFactory forceFactory;
	return &forceFactory;
}

ForceEvaluatorFactory::~ForceEvaluatorFactory()
{
	typeMap.clear();
}

IForceEvaluator* ForceEvaluatorFactory::Create(std::string evaluatorTypeName)
{
	TypeMap::iterator kvp = typeMap.find(evaluatorTypeName);

	if (kvp != typeMap.end()) {
		return kvp->second();
	}

	return NULL;
}