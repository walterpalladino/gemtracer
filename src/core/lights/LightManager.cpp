

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "core/lights/LightManager.h"

long LightManager::Add(Light light)
{
	lights.push_back(light);
	return (0);
}
