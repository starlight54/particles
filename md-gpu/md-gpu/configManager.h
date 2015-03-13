#ifndef CONFIG_MANAGER_H
#define CONFIG_MANAGER_H

//#include <tinyxml2.h>
#include <cuda_runtime_api.h>
//using namespace tinyxml2;

class ConfigManager
{
public:
	static ConfigManager* Get();
	~ConfigManager();
	void Initialise();
	cudaDeviceProp device;
private:
	//XMLDocument config;
};

#endif