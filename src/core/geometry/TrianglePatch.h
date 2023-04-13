

#pragma once

#include <string>
#include <vector>

using namespace std;

#include "core/geometry/Geometry.h"
#include "core/geometry/Triangle.h"

class Patch : public Geometry
{
public:
	vector<Triangle> triangles;
	Vector3d center;
	Triangle actual_hit;

public:
	long AddTriangle(Vector3d Vertex1,
					 Vector3d Vertex2,
					 Vector3d Vertex3);
	short WhoIAm(void) { return TRIANGLE_PATCH; };
	double Intersect(Vector3d rayvec, Vector3d raystart);
	void GetNormal(Vector3d Point, Vector3d *normal);
	void SetCenter(Vector3d CENTER);
	void Init(void);
};
