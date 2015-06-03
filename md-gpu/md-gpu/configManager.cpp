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
	//config.LoadFile("config.xml");
	//root = config.RootElement();
}

string ConfigManager::GetValue(string value)
{
	return "";
}