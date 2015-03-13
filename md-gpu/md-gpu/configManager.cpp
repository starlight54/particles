#include "configManager.h"

ConfigManager* ConfigManager::Get()
{
	static ConfigManager configManager;
	return &configManager;
}

ConfigManager::~ConfigManager()
{

}

void ConfigManager::Initialise()
{
	//GPU Init
	cudaGetDeviceProperties(&device, 0);
}