

#pragma once

#include "math/Vector3d.h"

#include <string>
#include <vector>

using namespace std;

class Triangle
{
public:
	Vector3d vertex[3];
	Vector3d normal[3];
	Vector3d transf[3];
	double u_hit,
		v_hit;

public:
	Triangle(void){};
	Triangle(Vector3d Vertex1, Vector3d Vertex2, Vector3d Vertex3);
	void Init(Vector3d Vertex1, Vector3d Vertex2, Vector3d Vertex3);
};
