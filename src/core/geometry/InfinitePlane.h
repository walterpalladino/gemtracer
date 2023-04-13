
#pragma once

#include "core/geometry/Geometry.h"
#include "math/Vector3d.h"

class InfinitePlane : public Geometry
{
public:
	Vector3d surface_normal;
	double distance;

public:
	InfinitePlane(void);
	InfinitePlane(Vector3d Normal, double Distance);
	void Init(Vector3d Normal, double Distance);

	short WhoIAm(void) { return (INFINITE_PLANE); };
	void GetNormal(Vector3d Point, Vector3d *normal);
	double Intersect(Vector3d rayvec, Vector3d raystart);
};
