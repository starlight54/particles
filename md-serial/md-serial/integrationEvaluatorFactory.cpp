#include "integrationEvaluatorFactory.h"

IntegrationEvaluatorFactory* IntegrationEvaluatorFactory::Get()
{
	static IntegrationEvaluatorFactory integratorFactory;
	return &integratorFactory;
}

IntegrationEvaluatorFactory::~IntegrationEvaluatorFactory()
{
	typeMap.clear();
}

IIntegrationEvaluator* IntegrationEvaluatorFactory::Create(std::string evaluatorTypeName)
{
	TypeMap::iterator kvp = typeMap.find(evaluatorTypeName);

	if (kvp != typeMap.end()) {
		return kvp->second();
	}

	return NULL;
}
