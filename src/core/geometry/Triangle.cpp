//
//
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "core/geometry/Triangle.h"

#include "math/Vector3d.h"
#include "core/renderer/Ray.h"

// Constructor
void Triangle::Init(Vector3d Vertex1,
					Vector3d Vertex2,
					Vector3d Vertex3)
{
	double k;
	Vector3d temp[3];
	Vector3d temp1,
		temp2;

	vertex[0] = Vertex1;
	vertex[1] = Vertex2;
	vertex[2] = Vertex3;

	/*
		normal[0]	= crossV( vertex[2] - vertex[0], vertex[1] - vertex[0] );
		normal[1]	= crossV( vertex[0] - vertex[1], vertex[2] - vertex[1] );
		normal[2]	= crossV( vertex[1] - vertex[2], vertex[0] - vertex[2] );

		normal[0] = normalize( normal[0] );
		normal[1] = normalize( normal[1] );
		normal[2] = normalize( normal[2] );
	*/

	temp1 = vertex[2] - vertex[0];
	temp2 = vertex[1] - vertex[0];
	normal[0] = temp1.Cross(temp2);

	temp1 = vertex[0] - vertex[1];
	temp2 = vertex[2] - vertex[1];
	normal[1] = temp1.Cross(temp2);

	temp1 = vertex[1] - vertex[2];
	temp2 = vertex[0] - vertex[2];
	normal[2] = temp1.Cross(temp2);

	normal[0].Normalize();
	normal[1].Normalize();
	normal[2].Normalize();

	temp[0] = vertex[1] - vertex[0];
	temp[1] = vertex[2] - vertex[0];

	temp[2] = temp[0].Cross(temp[1]);
	temp[2].Normalize();

	k = 1.0 /
		(temp[0].x * temp[1].y * temp[2].z + temp[0].y * temp[1].z * temp[2].x +
		 temp[0].z * temp[1].x * temp[2].y - temp[0].z * temp[1].y * temp[2].x -
		 temp[0].x * temp[1].z * temp[2].y - temp[0].y * temp[1].x * temp[2].z);

	transf[0].x = (temp[1].y * temp[2].z - temp[1].z * temp[2].y);
	transf[1].x = -(temp[0].y * temp[2].z - temp[0].z * temp[2].y);
	transf[2].x = (temp[0].y * temp[1].z - temp[0].z * temp[1].y);
	transf[0].y = -(temp[1].x * temp[2].z - temp[1].z * temp[2].x);
	transf[1].y = (temp[0].x * temp[2].z - temp[0].z * temp[2].x);
	transf[2].y = -(temp[0].x * temp[1].z - temp[0].z * temp[1].x);
	transf[0].z = (temp[1].x * temp[2].y - temp[1].y * temp[2].x);
	transf[1].z = -(temp[0].x * temp[2].y - temp[0].y * temp[2].x);
	transf[2].z = (temp[0].x * temp[1].y - temp[0].y * temp[1].x);

	transf[0] = transf[0] * k;
	transf[1] = transf[1] * k;
	transf[2] = transf[2] * k;
}
