
#pragma once

#include <vector>
#include <string>

#include "core/lights/Light.h"

using namespace std;

#define MAX_LIGHTS_PER_SCENE 32

class LightManager
{
public:
	vector<Light> lights;

public:
	static LightManager *GetInstance()
	{
		static LightManager instance;
		return &instance;
	};

	long Add(Light pLight);
};
