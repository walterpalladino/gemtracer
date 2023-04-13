//
//
//

#include <math.h>
#include <float.h>

#include "core/geometry/InfinitePlane.h"
#include "core/renderer/Ray.h"

/*
//  Initialize Data
void plane_init (	TPLANE		*plane,
					Vector3d		Center,
					double		Radius )
{
	sphere->center	= Center;
	sphere->radius	= Radius;
}

//  Identify object type
short sphere_who_i_am ( void )
{
	return ( SPHERE );
}
*/
/*
//  Calc Intersection Point
double  sphere_intersect (	TSPHERE		*sphere,
							Vector3d		rayvec,
							Vector3d		raystart )
{
double	Xc, Yc, Zc, Sr,
		X0, Y0, Z0,
		Xd, Yd, Zd,
		b, c,
		discriminant,
		sd,
		t, t0, t1;


	t = RAY_INFINITE;

	Xc = sphere->center.x;
	Yc = sphere->center.y;
	Zc = sphere->center.z;
	Sr = sphere->radius;

	X0 = raystart.x;
	Y0 = raystart.y;
	Z0 = raystart.z;

	Xd = rayvec.x;
	Yd = rayvec.y;
	Zd = rayvec.z;

	b = 2.0*(Xd*(X0-Xc)+Yd*(Y0-Yc)+Zd*(Z0-Zc));
	c = (X0-Xc)*(X0-Xc)+(Y0-Yc)*(Y0-Yc)+(Z0-Zc)*(Z0-Zc)-Sr*Sr;
	discriminant = b*b-4*c;

	if (discriminant >= 0)
	{
		sd = sqrt(discriminant);
		t0 = (-b-sd)/2;
		t1 = (-b+sd)/2;
		if (((t0 > 0.0) || (t1 > 0.0)) && (t0 != t1 ))
		{
			if (t0 > 0.0)
				t = t0;
			if ((t1 < t0) && (t1 >= 0.0))
				t = t1;
		}
	}

	return ( t );
}

void sphere_get_Normal (	TSPHERE		*sphere,
							Vector3d		Point,
							Vector3d		*normal )
{

	normal->x = (Point.x - sphere->center.x) / sphere->radius;
	normal->y = (Point.y - sphere->center.y) / sphere->radius;
	normal->z = (Point.z - sphere->center.z) / sphere->radius;

	normalize ( *normal, normal ) ;

}
*/

//  Constructor
InfinitePlane::InfinitePlane(void)
{
	surface_normal.Set(1.0, 1.0, 1.0);
	distance = 100.0;
}

InfinitePlane::InfinitePlane(Vector3d pNormal,
							 double pDistance)
{
	Init(pNormal, pDistance);
}

//  Initialize Data
void InfinitePlane::Init(Vector3d pNormal,
						 double pDistance)
{
	surface_normal = pNormal;
	distance = pDistance;
}

void InfinitePlane::GetNormal(Vector3d Point, Vector3d *normal)
{

	*normal = surface_normal;
	//	normalize ( *normal, normal ) ;
	(*normal).Normalize();
}

//  Calc Intersection Point
double InfinitePlane::Intersect(Vector3d rayvec, Vector3d raystart)
{
	double v0, vd, t;

	t = RAY_INFINITE;

	//	vd = dotV ( surface_normal, rayvec ) ;
	vd = surface_normal.Dot(rayvec);
	if (vd != 0)
	{
		//		v0 = - ( dotV( surface_normal, raystart ) + distance ) ;
		v0 = -(surface_normal.Dot(raystart) + distance);
		t = (v0 / vd);
	}

	return (t);
}
