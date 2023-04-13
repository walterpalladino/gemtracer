//
//
//

#include <math.h>
#include <float.h>

/*
#include    "math.hpp"
#include    "vec_math.hpp"

#include    "types.hpp"
#include    "ray.hpp"

#include    "raytracer.hpp"
#include    "sphere.hpp"
*/

#include "core/geometry/Sphere.h"
// #include	"math.h"
#include "core/renderer/Ray.h"

//  Constructor
Sphere::Sphere(void)
{
	//	setV ( &center, 0.0, 0.0, 0.0 );
	center.Set(0.0, 0.0, 0.0);
	radius = 100.0;
}

Sphere::Sphere(Vector3d pCenter, double pRadius)
{
	Init(pCenter, pRadius);
}

//  Initialize Data
void Sphere::Init(Vector3d pCenter, double pRadius)
{
	center = pCenter;
	radius = pRadius;
}

//  Identify object type
// short sphere_who_i_am ( void )
//{
//	return ( SPHERE );
//}

//  Calc Intersection Point
double Sphere::Intersect(Vector3d rayvec,
						 Vector3d raystart)
{
	double Xc, Yc, Zc, Sr,
		X0, Y0, Z0,
		Xd, Yd, Zd,
		b, c,
		discriminant,
		sd,
		t, t0, t1;

	t = RAY_INFINITE;

	Xc = center.x;
	Yc = center.y;
	Zc = center.z;
	Sr = radius;

	X0 = raystart.x;
	Y0 = raystart.y;
	Z0 = raystart.z;

	Xd = rayvec.x;
	Yd = rayvec.y;
	Zd = rayvec.z;

	b = 2.0 * (Xd * (X0 - Xc) + Yd * (Y0 - Yc) + Zd * (Z0 - Zc));
	c = (X0 - Xc) * (X0 - Xc) + (Y0 - Yc) * (Y0 - Yc) + (Z0 - Zc) * (Z0 - Zc) - Sr * Sr;
	discriminant = b * b - 4 * c;

	if (discriminant >= 0)
	{
		sd = sqrt(discriminant);
		t0 = (-b - sd) / 2;
		t1 = (-b + sd) / 2;
		if (((t0 > 0.0) || (t1 > 0.0)) && (t0 != t1))
		{
			if (t0 > 0.0)
				t = t0;
			if ((t1 < t0) && (t1 >= 0.0))
				t = t1;
		}
	}

	return (t);
}

void Sphere::GetNormal(Vector3d Point, Vector3d *normal)
{

	normal->x = (Point.x - center.x) / radius;
	normal->y = (Point.y - center.y) / radius;
	normal->z = (Point.z - center.z) / radius;

	//	normalize ( *normal, normal ) ;
	(*normal).Normalize();
}
