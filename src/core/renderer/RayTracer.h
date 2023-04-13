
#pragma once

#include "core/materials/Material.h"
#include "math/Vector3d.h"

class RayTracer
{
public:
	Vector3d ambientColor;
	Vector3d backgroundColor;
	check_type hit;

public:
	RayTracer(void);

	static RayTracer *GetInstance()
	{
		static RayTracer instance;
		return &instance;
	};

	void Shade(Vector3d loc,
			   Vector3d rayvec,
			   short depth,
			   short inside_object,
			   check_type *newcolor);

	void TraceRay(Vector3d rayvec,
				  Vector3d raystart,
				  short depth,
				  short inside_object,
				  check_type *Check);

	Vector3d DoDither(Vector3d color, double percentage);
};
