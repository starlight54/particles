#ifndef CONFIG_MANAGER_H
#define CONFIG_MANAGER_H

//#include <tinyxml2.h>
#include <string>
#include <cuda_runtime_api.h>
//using namespace tinyxml2;
using namespace std;

class ConfigManager
{
public:
	static ConfigManager* Get();
	~ConfigManager();
	void Initialise();
	cudaDeviceProp device;
	string GetValue(string value);
private:
	//XMLDocument config;
	//XMLElement* root;
};

#endif