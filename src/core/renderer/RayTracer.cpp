
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "core/geometry/Geometry.h"
#include "core/geometry/Sphere.h"
#include "core/geometry/InfinitePlane.h"

#include "core/lights/LightManager.h"

#include "core/renderer/Ray.h"
#include "core/renderer/RayTracer.h"

// #include	"math.h"
#include "core/materials/Material.h"
#include "core/materials/Bump.h"

// Vector3d	BACKGROUND  ( 0.5f, 0.5f, 1.0f );
Vector3d NOTHING(-1.0f, -1.0f, -1.0f);
// Vector3d	AMBIENT     ( 0.6f, 0.6f, 0.6f );
Vector3d ZERO(0.0f, 0.0f, 0.0f);
Vector3d ORIGIN(0.0f, 100.0f, -800.0f);
Vector3d BLACK(0.0f, 0.0f, 0.0f);
Vector3d WHITE(1.0f, 1.0f, 1.0f);
Vector3d RED(1.0f, 0.0f, 0.0f);
Vector3d GREEN(0.0f, 1.0f, 0.0f);
Vector3d BLUE(0.0f, 0.0f, 1.0f);
Vector3d YELLOW(1.0f, 1.0f, 0.0f);
Vector3d MAGENTA(1.0f, 0.0f, 1.0f);
Vector3d CYAN(0.0f, 1.0f, 1.0f);

// Frame frame;

RayTracer::RayTracer(void)
{
	rand();

	ambientColor.Set(0.6f, 0.6f, 0.6f);
	backgroundColor.Set(0.5f, 0.5f, 1.0f);
}

void RayTracer::TraceRay(Vector3d rayvec,
						 Vector3d raystart,
						 short depth,
						 short inside_object,
						 check_type *Check)
{
	double t;
	double prev_t;
	Vector3d norm, loc; // plane_loc,pixel,diff,temp,plane;

	Vector3d bnormal;

	Check->color = backgroundColor;
	Check->object_hit = NO_HIT;
	prev_t = RAY_INFINITE;

	for (std::vector<Geometry *>::const_iterator it = GeometryManager::GetInstance()->objects.begin(); it != GeometryManager::GetInstance()->objects.end(); ++it)
	{
		t = (*it)->Intersect(rayvec, raystart);

		if ((t < prev_t) && (t > RAY_EPSILON)) //  && (t<RAY_INFINITE))
		{

			prev_t = t;
			//			addmulV ( raystart, rayvec, t, &loc ) ;
			loc = raystart + (rayvec * t);

			(*it)->GetNormal(loc, &norm);

			// Check if surface is affected by any type of BUMP
			if ((*it)->material->bump_type != NO_BUMP)
			{
				//				normalize ( norm, &norm ) ;
				norm.Normalize();
				get_Bump_Normal((TWAVE *)(*it)->material->bump_data, loc, &bnormal);
				//				addV ( norm, bnormal, &norm ) ;
				norm = norm + bnormal;
			}

			//			normalize ( norm, &norm ) ;
			norm.Normalize();
			Check->normal = norm;
			Check->material = (*it)->material;

			Check->object_hit = (*it)->WhoIAm();
		}

		if (prev_t < RAY_INFINITE)
		// if (( prev_t < RAY_INFINITE) &&
		//	(depth > 0))
		{
			Shade(loc,
				  rayvec,
				  depth,
				  inside_object,
				  Check);
			Check->distance = prev_t;
		}
	}

	hit = *Check;
}

Vector3d RayTracer::DoDither(Vector3d color, double percentage)
{

	double randval, integer;
	Vector3d newcolor;

	if (percentage == 0.0)
	{
		newcolor.x = color.x;
		newcolor.y = color.y;
		newcolor.z = color.z;
	}
	else
	{
		randval = (10000 * rand()) / RAND_MAX;

		if (modf((255.0f / percentage) * color.x, &integer) * 10000 > randval)
			newcolor.x = ceil((255.0f / percentage) * color.x) / (255.0f / percentage);
		else
			newcolor.x = floor((255.0f / percentage) * color.x) / (255.0f / percentage);
		if (modf((255.0 / percentage) * color.y, &integer) * 10000 > randval)
			newcolor.y = ceil((255.0f / percentage) * color.y) / (255.0f / percentage);
		else
			newcolor.y = floor((255.0f / percentage) * color.y) / (255.0f / percentage);
		if (modf((255.0 / percentage) * color.z, &integer) * 10000 > randval)
			newcolor.z = ceil((255.0f / percentage) * color.z) / (255.0f / percentage);
		else
			newcolor.z = floor((255.0f / percentage) * color.z) / (255.0f / percentage);
	}

	return newcolor;
}

void RayTracer::Shade(Vector3d loc,
					  Vector3d rayvec,
					  short depth,
					  short inside_object,
					  check_type *newcolor)
{

	check_type reflected,
		refracted,
		shadow;

	Vector3d color;
	double dot;
	double dot_reflected;
	Vector3d reflected_dir;
	Vector3d refracted_dir;
	Vector3d lightvec_dir;
	Vector3d light_color;
	Vector3d vtemp;

	double intensity;
	Vector3d kh;

	Material *MATERIAL;
	Vector3d normal;

	Light light;

	normal = newcolor->normal;
	MATERIAL = newcolor->material;

	newcolor->color = BLACK;
	reflected.color = BLACK;
	refracted.color = BLACK;

	MATERIAL->GetColor(loc, &color);

	dot = rayvec.Dot(normal);

	if (dot > 0.0)
	{
		dot = -dot;
		normal = -normal;
	}

	// Compute reflected direction
	vtemp = normal * (-2.0f * dot);
	reflected_dir = rayvec + vtemp;
	reflected_dir.Normalize();

	//	Light Sources -->
	Vector3d difusse_impact;
	Vector3d specular_impact;
	Vector3d shadow_impact;

	for (size_t l = 0; l < LightManager::GetInstance()->lights.size(); l++)
	{
		difusse_impact.Set(0.0f, 0.0f, 0.0f);
		specular_impact.Set(0.0f, 0.0f, 0.0f);
		shadow_impact.Set(0.0f, 0.0f, 0.0f);

		light = LightManager::GetInstance()->lights[l];

		shadow.color = WHITE;

		lightvec_dir = light.position - loc;
		lightvec_dir.Normalize();

		intensity = normal.Dot(lightvec_dir);

		if (intensity < 0.0)
			intensity = 0;
		if (intensity > 1.0)
			intensity = 1.0;

		dot_reflected = -reflected_dir.Dot(rayvec);

		// if ((dot_reflected > 0.0 ) && ( testNullV ( MATERIAL->specular) == 0))
		if ((dot_reflected > 0.0) &&
			((MATERIAL->specular.x != 0.0) ||
			 (MATERIAL->specular.y != 0.0) ||
			 (MATERIAL->specular.z != 0.0)))
			kh = MATERIAL->specular * pow(dot_reflected, MATERIAL->ES);
		else
			kh.Set(0.0, 0.0, 0.0);

		// Calc effect of SHADOW
		if (depth)
		{

			TraceRay(lightvec_dir,
					 loc,
					 0,
					 inside_object,
					 &shadow);

			if (shadow.object_hit == NO_HIT)
				shadow.color = WHITE;
			//				shadow.color = light.color ;
			else if (shadow.material->ir == 0.0)
				shadow.color = BLACK;
			// else
			// if (testNullV(shadow.material->reflected) == 0.0)
			//	shadow.color = BLACK ;
		}

		//	Shadow
		shadow_impact = light.color * shadow.color;
		//	Difusse
		difusse_impact = shadow_impact * MATERIAL->difusse;
		difusse_impact = difusse_impact * intensity;
		//	Specular
		//		if ((dot_reflected > 0.0 ) && ( testNullV ( MATERIAL->specular) == 0))
		if ((dot_reflected > 0.0) &&
			((MATERIAL->specular.x != 0.0) ||
			 (MATERIAL->specular.y != 0.0) ||
			 (MATERIAL->specular.z != 0.0)))
			specular_impact = MATERIAL->specular * pow(dot_reflected, MATERIAL->ES);
		else
			specular_impact.Set(0.0, 0.0, 0.0);
		//		mulV ( light.color, MATERIAL->specular, &specular_impact ) ;
		//		mulS ( specular_impact, dot_reflected, &specular_impact ) ;

		// Effect of SHADOW
		light_color = light.color * shadow.color;

		// Effect of LIGHT SOURCE
		light_color = light_color * intensity;

		// Affect diffuse factor
		light_color = light_color * MATERIAL->difusse;

		// Affect specular factor
		if (MATERIAL->METAL)
		{
			vtemp = color * kh;
			light_color = light_color + vtemp;

			specular_impact = color * specular_impact;
		}
		else
		{
			vtemp = light.color * kh;
			light_color = light_color + vtemp;

			vtemp = light.color * specular_impact;
			specular_impact = specular_impact + vtemp;
		}

		newcolor->color = newcolor->color + light_color;

		// newcolor->color	= newcolor->color + difusse_impact ;
		// newcolor->color	= newcolor->color + specular_impact ;
	}
	//	<-- Light Sources

	newcolor->color = newcolor->color * color;

	// Compute refracted light
	if (MATERIAL->ir != 0)
	{
		double eta,
			c1, cs2;

		if (inside_object == OUT_OBJECT)
		{
			eta = AIR_IR / MATERIAL->ir;
		}
		else
		{
			eta = MATERIAL->ir / AIR_IR;
		}
		if (eta == 1.0)
		{
			if (inside_object == OUT_OBJECT)
				TraceRay(rayvec, loc, depth, IN_OBJECT, &refracted);
			else
				TraceRay(rayvec, loc, depth, OUT_OBJECT, &refracted);
		}
		else
		{
			if (dot < 0.0)
			{
				dot = -dot;
				normal = -normal;
			}

			c1 = dot;
			cs2 = 1.0 - eta * eta * (1.0 - c1 * c1);

			if (cs2 >= 0.0)
			{
				cs2 = eta * c1 - sqrt(cs2);
				vtemp = normal * cs2;
				refracted_dir = rayvec * eta;
				refracted_dir = vtemp + refracted_dir;

				refracted_dir.Normalize();

				if (inside_object == OUT_OBJECT)
				{
					TraceRay(vtemp, loc, depth, IN_OBJECT, &refracted);
				}
				else
				{
					TraceRay(vtemp, loc, depth, OUT_OBJECT, &refracted);
				}
			}
			else if (depth)
			{
				// TOTAL INTERNAL REFRACTION
				if (dot > 0.0)
				{
					dot = -dot;
					// negV ( normal, &normal ) ;
					normal = -normal;
				}

				vtemp = normal * (-2 * dot);
				refracted_dir = rayvec + vtemp;
				vtemp.Normalize();
				TraceRay(refracted_dir, loc, --depth, inside_object, &refracted);
			}
		}
	}

	// Compute reflected light
	//	if (!testNullV(MATERIAL->reflected) && (depth != 0))// && (inside_object == OUT_OBJECT))
	if ((depth != 0) &&
		((MATERIAL->reflected.x != 0.0) ||
		 (MATERIAL->reflected.y != 0.0) ||
		 (MATERIAL->reflected.z != 0.0)))
		TraceRay(reflected_dir, loc, --depth, inside_object, &reflected);

	// Add effect of Reflected light
	vtemp = reflected.color * MATERIAL->reflected;
	newcolor->color = newcolor->color + vtemp;

	// Add effect of Refracted light
	vtemp = refracted.color * MATERIAL->refracted;
	newcolor->color = newcolor->color + vtemp;

	// Add effect of LUMINOSITY
	vtemp = color * MATERIAL->kl;
	newcolor->color = newcolor->color + vtemp;

	// Add AMBIENT color
	//	addV ( newcolor->color, AMBIENT, &newcolor->color ) ;
	vtemp = color * ambientColor;
	newcolor->color = newcolor->color + vtemp;

	//  DITHERING
	//	newcolor->color	= DoDither ( newcolor->color, 0.0 ) ;

	// Clip values of COLOR
	newcolor->color.Clip();

	return;
}
