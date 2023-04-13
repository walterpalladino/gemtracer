//
//
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "core/geometry/TrianglePatch.h"

#include "math/Vector3d.h"
#include "core/renderer/Ray.h"

//
long Patch::AddTriangle(Vector3d Vertex1,
						Vector3d Vertex2,
						Vector3d Vertex3)
{
	Triangle triangle;

	//	triangle		= (TTRIANGLE *)malloc ( sizeof ( TTRIANGLE ) ) ;
	//	triangle		= new ( TTRIANGLE ) ;

	triangle.Init(Vertex1, Vertex2, Vertex3);

	triangles.push_back(triangle);

	return 0;
}

//  Calc Intersection Point
double Patch::Intersect(Vector3d rayvec, Vector3d raystart)
{
	double t = RAY_INFINITE;
	double distance, k, k0;
	Vector3d p;
	Triangle triangle;
	// printf("Entre a Intersect\n");

	for (size_t tp = 0; tp < triangles.size(); tp++)
	{
		triangle = triangles[tp];

		//		k = dotV( rayvec, triangle.transf[2]);
		k = rayvec.Dot(triangle.transf[2]);
		// printf("k %f ", k );
		if (fabs(k) >= RAY_EPSILON)
		{
			//			p = vertex[0] - raystart;
			//			subV ( triangle.vertex[0], raystart, &p ) ;
			p = triangle.vertex[0] - raystart;

			//			distance = fabs(dotV( p, triangle.transf[2]) / k);
			distance = fabs(p.Dot(triangle.transf[2]) / k);
			// printf("distance %f ", distance );
			if (distance > RAY_EPSILON)
			{
				//				p = (rayvec * distance) + raystart - vertex[0];
				// mulS( rayvec, distance, &p ) ;
				p = rayvec * distance;
				// addV( p, raystart, &p ) ;
				p = p + raystart;
				// subV( p, triangle.vertex[0], &p ) ;
				p = p - triangle.vertex[0];

				// k  = dotV( p, triangle.transf[0] );
				k = p.Dot(triangle.transf[0]);
				// k0 = dotV( p, triangle.transf[1] );
				k0 = p.Dot(triangle.transf[1]);
				// printf("k %f k0 %f k+k0 %f\n", k, k0, k+k0 );
				if ((k >= 0.0) && (k0 >= 0.0) && ((k + k0) <= 1.0))
				{
					triangle.u_hit = k;
					triangle.v_hit = k0;
					if (distance < t)
					{
						t = distance;
						actual_hit = triangle;
					}
				}
			}
		}
	}
	// printf("t : %e\n",t);
	return (t);
}

void Patch::GetNormal(Vector3d Point, Vector3d *normal)
{
	double t,
		u,
		v;

	u = actual_hit.u_hit;
	v = actual_hit.v_hit;
	t = 1.0 - u - v;

	normal->x = t * actual_hit.normal[0].x + u * actual_hit.normal[1].x +
				v * actual_hit.normal[2].x;
	normal->y = t * actual_hit.normal[0].y + u * actual_hit.normal[1].y +
				v * actual_hit.normal[2].y;
	normal->z = t * actual_hit.normal[0].z + u * actual_hit.normal[1].z +
				v * actual_hit.normal[2].z;

	//	return ( normalize( normal ) );
	//    normal = normal / lenV( normal, ORIGIN );

	//    return ( normal );

	//	normalize( *normal, normal ) ;
	(*normal).Normalize();
}

void Patch::SetCenter(Vector3d CENTER)
{
	center.x = CENTER.x;
	center.y = CENTER.y;
	center.z = CENTER.z;
}

void Patch::Init(void)
{

	SetCenter(Vector3d(0.0, 0.0, 0.0));
	//	actual_hit	= NULL ;
}
