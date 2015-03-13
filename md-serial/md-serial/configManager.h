#ifndef CONFIG_MANAGER_H
#define CONFIG_MANAGER_H

#include <tinyxml2.h>
using namespace tinyxml2;

class ConfigManager
{
public:
	static ConfigManager* Get();
	~ConfigManager();
	void Initialise();
private:
	XMLDocument config;
};

#endif