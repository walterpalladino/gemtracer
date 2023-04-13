
#pragma once

#include "core/geometry/Geometry.h"
#include "math/Vector3d.h"

class Sphere : public Geometry
{
public:
	Vector3d center; // Center of sphere
	double radius;	 // Radius of sphere

public:
	Sphere(void);
	Sphere(Vector3d pCenter, double pradius);
	void Init(Vector3d pCenter, double pRadius);

	short WhoIAm(void) { return (SPHERE); };
	void GetNormal(Vector3d Point, Vector3d *normal);
	double Intersect(Vector3d rayvec, Vector3d raystart);
};
