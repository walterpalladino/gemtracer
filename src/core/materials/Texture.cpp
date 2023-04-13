

#include <stdio.h>
#include <math.h>
#include <float.h>

#include "core/materials/Material.h"
#include "core/materials/Perlin.h"

// extern Vector3d ZERO;

double Chaos(Vector3d V, int octaves)
{
	double s, t;
	Vector3d tp;

	s = 1.0;
	t = 0.0;
	tp = V;

	while (octaves--)
	{
		t += fabsf(float(noise(tp.x, tp.y, tp.z))) * s;
		s *= 0.5;
		tp.x *= 2.;
		tp.y *= 2.;
		tp.z *= 2.;
	}

	return t;
}

/*
double Chaos ( Vector3d V, int octaves )
{
double		s, t;
float		tp [3] ;

	s = 1.0;
	t = 0.0;
	tp[0] = V.x ;
	tp[1] = V.y ;
	tp[2] = V.z ;

	while (octaves--)
	{
		t += Noise3( tp ) * s;
		s *= 0.5;
		tp[0] *= 2.;
		tp[1] *= 2.;
		tp[2] *= 2.;
	}

	return t;
}*/

void Flat::GetColor(Vector3d POINT, Vector3d *pColor)
{
	pColor->x = color.x;
	pColor->y = color.y;
	pColor->z = color.z;
}

//
//  MOSAIC
//
void Mosaic::Init(void)
{
	center.x = 10.0;
	center.y = 10.0;
	center.z = 10.0;

	color1.x = 1.0;
	color1.y = 0.0;
	color1.z = 1.0;

	color2.x = 0.0;
	color2.y = 0.0;
	color2.z = 1.0;

	scale.x = 50.0;
	scale.y = 50.0;
	scale.z = 50.0;

	mortar = 0.15;
}

void Mosaic::SetCenter(Vector3d CENTER)
{
	center = CENTER;
}

void Mosaic::SetColor(Vector3d COLOR1, Vector3d COLOR2)
{
	color1 = COLOR1;
	color2 = COLOR2;
}

void Mosaic::SetScale(Vector3d SCALE)
{
	scale = SCALE;
}

void Mosaic::SetMortar(double MORTAR)
{
	mortar = MORTAR;
}

void Mosaic::GetColor(Vector3d POINT, Vector3d *color)
{

	POINT.x = POINT.x - center.x;
	POINT.y = POINT.y - center.y;
	POINT.z = POINT.z - center.z;

	POINT.x = fabs(fmod(POINT.x + .001f, scale.x));
	POINT.y = fabs(fmod(POINT.y + .001f, scale.y));
	POINT.z = fabs(fmod(POINT.z + .001f, scale.z));

	if (((POINT.x > 0.0) && (POINT.x < mortar * scale.x)) ||
		((POINT.y > 0.0) && (POINT.y < mortar * scale.y)) ||
		((POINT.z > 0.0) && (POINT.z < mortar * scale.z)))
	{
		color->x = color1.x;
		color->y = color1.y;
		color->z = color1.z;
	}
	else
	{
		color->x = color2.x;
		color->y = color2.y;
		color->z = color2.z;
	}

	return;
}

//
//  CHECKER
//
void Checker::Init(void)
{

	center.x = 0.0;
	center.y = 0.0;
	center.z = 0.0;

	color1.x = 0.0;
	color1.y = 0.0;
	color1.z = 1.0;

	color2.x = 1.0;
	color2.y = 1.0;
	color2.z = 0.0;

	scale.x = 200.00;
	scale.y = 200.00;
	scale.z = 200.00;
}

void Checker::SetScale(Vector3d pScale)
{

	scale.x = pScale.x;
	scale.y = pScale.y;
	scale.z = pScale.z;
}

void Checker::SetColor(Vector3d pColor1, Vector3d pColor2)
{

	color1.x = pColor1.x;
	color1.y = pColor1.y;
	color1.z = pColor1.z;

	color2.x = pColor2.x;
	color2.y = pColor2.y;
	color2.z = pColor2.z;
}

void Checker::GetColor(Vector3d POINT, Vector3d *color)
{
	int checkindex;

	POINT.x = POINT.x - center.x;
	POINT.y = POINT.y - center.y;
	POINT.z = POINT.z - center.z;

	if (scale.x != 0)
		POINT.x = (POINT.x + .001f) / scale.x;
	else
		POINT.x = 0;
	if (scale.y != 0)
		POINT.y = (POINT.y + .001f) / scale.y;
	else
		POINT.y = 0;
	if (scale.z != 0)
		POINT.z = (POINT.z + .001f) / scale.z;
	else
		POINT.z = 0;

	checkindex = (int)(floor(POINT.x) + floor(POINT.y) + floor(POINT.z));

	if (checkindex & 1)
	{
		color->x = color1.x;
		color->y = color1.y;
		color->z = color1.z;
	}
	else
	{
		color->x = color2.x;
		color->y = color2.y;
		color->z = color2.z;
	}
}

void Brick::Init(void)
{

	center.x = 0.0f;
	center.y = 0.0f;
	center.z = 0.0f;

	color1.x = 0.6f;
	color1.y = 0.31f;
	color1.z = 0.11f;

	color2.x = 0.0f;
	color2.y = 0.0f;
	color2.z = 0.0f;

	scale.x = 50.0f;
	scale.y = 50.0f;
	scale.z = 250.0f;

	mortar = 0.10f;
}

void Brick::SetCenter(Vector3d CENTER)
{
	center = CENTER;
}

void Brick::SetColor(Vector3d COLOR1, Vector3d COLOR2)
{
	color1 = COLOR1;
	color2 = COLOR2;
}

void Brick::SetScale(Vector3d SCALE)
{
	scale = SCALE;
}

void Brick::SetMortar(double MORTAR)
{
	mortar = MORTAR;
}

void Brick::GetColor(Vector3d POINT, Vector3d *color)
{
	// double xr, yr, zr;

	color->x = color2.x;
	color->y = color2.y;
	color->z = color2.z;

	POINT.x = POINT.x - center.x;
	POINT.y = POINT.y - center.y;
	POINT.z = POINT.z - center.z;

	if (scale.x != 0)
		POINT.x = (POINT.x + .001) / scale.x;
	else
		POINT.x = 0;
	if (scale.y != 0)
		POINT.y = (POINT.y + .001) / scale.y;
	else
		POINT.y = 0;
	if (scale.z != 0)
		POINT.z = (POINT.z + .001) / scale.z;
	else
		POINT.z = 0;

	POINT.x = fabs(fmod(POINT.x, 1.0f));
	POINT.y = fabs(fmod(POINT.y, 1.0f));
	POINT.z = fabs(fmod(POINT.z, 1.0f));

	if (POINT.y < (0.5f - mortar))
	{
		if ((POINT.x > (mortar / 2.0f)) && (POINT.x < (1.0f - mortar / 2.0f)))
		{
			if (((POINT.z > mortar) && (POINT.z < 0.5f)) ||
				((POINT.z > (0.5f + mortar)) && (POINT.z < 1.0f)))
			{
				color->x = color1.x;
				color->y = color1.y;
				color->z = color1.z;
			}
		}
	}
	else if ((POINT.y > 0.5) && (POINT.y < (1.0 - mortar)))
	{
		if ((POINT.z > (mortar / 2.0)) && (POINT.z < (1.0 - mortar / 2.0)))
		{
			if (((POINT.x > 0.0) && (POINT.x < (0.5 - mortar / 2.0))) ||
				((POINT.x > (0.5 + mortar / 2.0)) && (POINT.x < 1.0)))
			{
				color->x = color1.x;
				color->y = color1.y;
				color->z = color1.z;
			}
		}
	}
}

//
//  SPHERIC
//
void Spheric::Init(void)
{
	//	setV( &center, 0.0, 0.0, 0.0 );
	center.Set(0.0, 0.0, 0.0);
	radius1 = 10;
	radius2 = 20;
	//	setV( &color1, 1.0, 1.0, 1.0 );
	color1.Set(1.0, 1.0, 1.0);
	//	setV( &color2, 0.0, 0.0, 0.0 );
	color2.Set(0.0, 0.0, 0.0);
	//	setV( &scale, 1.0, 1.0, 1.0 );
	scale.Set(1.0, 1.0, 1.0);
	blend = 0.1;
}

void Spheric::SetCenter(Vector3d CENTER)
{
	center = CENTER;
}

void Spheric::SetColor(Vector3d COLOR1, Vector3d COLOR2)
{
	color1 = COLOR1;
	color2 = COLOR2;
}

void Spheric::SetRadius(double RADIUS1, double RADIUS2)
{
	radius1 = RADIUS1;
	radius2 = RADIUS2;
}

void Spheric::SetScale(Vector3d SCALE)
{
	scale = SCALE;
}

void Spheric::SetBlend(double BLEND)
{
	blend = BLEND;
}

void Spheric::GetColor(Vector3d POINT, Vector3d *COLOR)
{

	double r, dist, value;
	Vector3d point;

	//    subV ( POINT, center, &point ) ;
	point = POINT - center;

	if (scale.x != 0)
		point.x = (point.x + .001f) / scale.x;
	else
		point.x = 0;
	if (scale.y != 0)
		point.y = (point.y + .001f) / scale.y;
	else
		point.y = 0;
	if (scale.z != 0)
		point.z = (point.z + .001f) / scale.z;
	else
		point.z = 0;

	//    dist = lenV( point, ZERO );
	dist = (point - Vector3d::Zero()).Magnitude();

	r = radius1 + radius2;

	dist += (radius1 / 2.0f);
	dist = fmod(dist, r);

	if (dist < radius1)
	{
		dist = dist / radius1;
		if (dist > 0.5)
		{
			dist = 1.0 - dist;
		}
		if (dist >= blend / 2.0)
		{
			value = 1.0;
		}
		else
		{
			value = 0.5 + dist / blend;
		}
	}
	else
	{
		dist = (dist - radius1) / radius2;
		if (dist > 0.5)
		{
			dist = 1.0 - dist;
		}
		if (dist >= blend / 2.0)
		{
			value = 0.0;
		}
		else
		{
			value = 0.5 - dist / blend;
		}
	}

	if (value > 1.0)
		value = 1.0;
	if (value < 0.0)
		value = 0.0;

	//	subV ( color2, color1, COLOR ) ;
	*COLOR = color2 - color1;
	//	addmulV ( color1, *COLOR, value, COLOR ) ;
	*COLOR = color1 + (*COLOR * value);

	// clipZEROV( *COLOR, COLOR );
	// clipONEV( *COLOR, COLOR );
	(*COLOR).Clip();
}

//
//  LEOPARD
//
void Leopard::Init(void)
{
	// setV( &center, 0.0, 0.0, 0.0 );
	center.Set(0.0, 0.0, 0.0);
	// setV( &color1, 1.0, 1.0, 1.0 );
	color1.Set(1.0, 1.0, 1.0);
	// setV( &color2, 0.0, 0.0, 0.0 );
	color2.Set(0.0, 0.0, 0.0);
	// setV( &scale, 1.0, 1.0, 1.0 );
	scale.Set(1.0, 1.0, 1.0);
	blend = 0.1;
}

void Leopard::SetCenter(Vector3d CENTER)
{
	center = CENTER;
}

void Leopard::SetColor(Vector3d COLOR1, Vector3d COLOR2)
{
	color1 = COLOR1;
	color2 = COLOR2;
}

void Leopard::SetScale(Vector3d SCALE)
{
	scale = SCALE;
}

void Leopard::SetBlend(double BLEND)
{
	blend = BLEND;
}

void Leopard::GetColor(Vector3d point, Vector3d *COLOR)
{

	double r, dist, value;

	Vector3d POINT;

	//    subV ( point, center, &POINT ) ;
	POINT = point - center;

	if (scale.x != 0)
		POINT.x = (POINT.x + .001f) / scale.x;
	else
		POINT.x = 0;
	if (scale.y != 0)
		POINT.y = (POINT.y + .001f) / scale.y;
	else
		POINT.y = 0;
	if (scale.z != 0)
		POINT.z = (POINT.z + .001f) / scale.z;
	else
		POINT.z = 0;

	value = (sin(POINT.x) + sin(POINT.y) + sin(POINT.z)) / 3.0;
	value = value * value;

	if (blend == 0.0)
	{
		if (value > 0.25)
			value = 1.0;
		else
			value = 0.0;
	}
	else
	{
		if (value > 0.75)
			value = 1.0;
		else if (value < 0.25)
			value = 0.0;
		else
			//            value = 0.25 + (value - 0.25) * blend / 0.5;
			value = 0.25 + (value - 0.25) * blend;
	}

	//	subV ( color2, color1, COLOR ) ;
	*COLOR = color2 - color1;
	//	addmulV ( color1, *COLOR, value, COLOR ) ;
	*COLOR = color1 + (*COLOR * value);

	// clipZEROV( *COLOR, COLOR );
	// clipONEV( *COLOR, COLOR );
	(*COLOR).Clip();
}

//
//	GRADIENT
//

void Gradient::Init(void)
{
	// setV( &center, 0.0, 0.0, 0.0 );
	center.Set(0.0, 0.0, 0.0);
	// setV( &scale, 1.0, 1.0, 1.0 );
	scale.Set(1.0, 1.0, 1.0);
	// setV( &gradient, 0.1, 0.0, 0.0 );
	gradient.Set(0.1, 0.0, 0.0);
	// setV( &color1, 0.0, 1.0, 0.0 );
	color1.Set(0.0, 1.0, 0.0);
	// setV( &color2, 0.0, 0.0, 1.0 );
	color2.Set(0.0, 0.0, 1.0);
	blend = 0.5;
}

void Gradient::SetGradient(Vector3d GRADIENT)
{
	gradient = GRADIENT;
}

void Gradient::SetCenter(Vector3d CENTER)
{
	center = CENTER;
}

void Gradient::SetColor(Vector3d COLOR1, Vector3d COLOR2)
{
	color1 = COLOR1;
	color2 = COLOR2;
}

void Gradient::SetScale(Vector3d SCALE)
{
	scale = SCALE;
}

void Gradient::SetBlend(double BLEND)
{
	blend = BLEND;
}

void Gradient::GetColor(Vector3d point, Vector3d *COLOR)
{

	double value;
	Vector3d POINT;

	value = 0.0;
	//	setZ ( COLOR ) ;
	(*COLOR).Set(0.0, 0.0, 0.0);
	POINT = point;

	//    subV ( point, center, &POINT ) ;
	POINT = point - center;

	if (scale.x != 0)
		POINT.x = (POINT.x + .001) / scale.x;
	else
		POINT.x = 0;
	if (scale.y != 0)
		POINT.y = (POINT.y + .001) / scale.y;
	else
		POINT.y = 0;
	if (scale.z != 0)
		POINT.z = (POINT.z + .001) / scale.z;
	else
		POINT.z = 0;

	if (gradient.x != 0.0)
	{
		POINT.x = fabs(POINT.x);
		// obtain fractional X component
		value += POINT.x - floor(POINT.x);
	}
	if (gradient.y != 0.0)
	{
		POINT.y = fabs(POINT.y);
		// obtain fractional Y component
		value += POINT.y - floor(POINT.y);
	}
	if (gradient.z != 0.0)
	{
		POINT.z = fabs(POINT.z);
		// obtain fractional Z component
		value += POINT.z - floor(POINT.z);
	}

	if (blend == 0.0)
	{
		if (value > 0.25)
			value = 1.0;
		else
			value = 0.0;
	}
	else
	{
		if (value > 0.75)
			value = 1.0;
		else if (value < 0.25)
			value = 0.0;
		else
			//            value = 0.25 + (value - 0.25) * blend / 0.5;
			value = 0.25 + (value - 0.25) * blend;
	}

	// subV ( color2, color1, COLOR ) ;
	// addmulV ( color1, *COLOR, value, COLOR ) ;
	*COLOR = color2 - color1;
	*COLOR = color1 + (*COLOR * value);
	// clipZEROV( *COLOR, COLOR );
	// clipONEV( *COLOR, COLOR );
	(*COLOR).Clip();
}

/*
static DBL onion (VECTOR EPoint)
{
	// The variable noise is not used as noise in this function

  register DBL noise;


	// This ramp goes 0-1, 0-1, 0-1, 0-1 ...

  noise = (fmod(sqrt(Sqr(EPoint[X])+Sqr(EPoint[Y])+Sqr(EPoint[Z])), 1.0));

  return(noise);
}
*/
/*
static DBL radial (VECTOR EPoint)
{
  register DBL value;

  if ((fabs(EPoint[X])<0.001) && (fabs(EPoint[Z])<0.001))
  {
	value = 0.25;
  }
  else
  {
	value = 0.25 + (atan2(EPoint[X],EPoint[Z]) + M_PI) / TWO_M_PI;
  }

  return(value);
}

*/

void Wood::SetColor(Vector3d COLOR1, Vector3d COLOR2)
{
	color1 = COLOR1;
	color2 = COLOR2;
}

void Wood::SetScale(Vector3d SCALE)
{
	scale = SCALE;
}

//
void Wood::GetColor(Vector3d point, Vector3d *COLOR)
{
	float perturb;
	float noiseCoef;

	perturb = fabsf(float(noise(point.x * scale.x,
								point.y * scale.y,
								point.z * scale.z)));
	noiseCoef = perturb * grain_scale;
	noiseCoef = noiseCoef - int(noiseCoef);

	*COLOR = (color1 * noiseCoef) + (color2 * (1.0f - noiseCoef));
	/*
		{
	double		s, t;
	Vector3d	tp;
	int			octaves	= 7 ;
	float		red, blu, grn ;
	float		brownPerturb,
				greenPerturb,
				grnPerturb ;
	float		greenLayer,
				brownLayer ;

		float	chaos ;

		tp		= point * 0.05f ;

		chaos	= Chaos ( tp, 7 ) ;

		t = sin(sin(8.*chaos + 7*point.x +3.*point.y));

		greenLayer = brownLayer = fabs(t);

		perturb = sin(40.*chaos + 50.*point.z);
		perturb = fabs(perturb);
	//noiseCoef	= perturb ;
		brownPerturb = .6*perturb + 0.3;
		greenPerturb = .2*perturb + 0.8;
		grnPerturb = .15*perturb + 0.85;
		grn = 0.5 * pow(fabs(brownLayer), 0.3f);
		brownLayer = pow(0.5 * (brownLayer+1.0), 0.6) * brownPerturb;
		greenLayer = pow(0.5 * (greenLayer+1.0), 0.6) * greenPerturb;

		red = (0.5*brownLayer + 0.35*greenLayer)*2.*grn;
		blu = (0.25*brownLayer + 0.35*greenLayer)*2.0*grn;
	//	grn *= max(brownLayer, greenLayer) * grnPerturb;
		if (brownLayer > greenLayer)
			grn *= brownLayer * grnPerturb;
		else
			grn *= greenLayer
			* grnPerturb;

		//surf->diff.r *= red;
		//surf->diff.g *= grn;
		//surf->diff.b *= blu;
		//surf->amb.r *= red;
		//surf->amb.g *= grn;
		//surf->amb.b *= blu;

		COLOR->x	= red ;
		COLOR->y	= grn ;
		COLOR->z	= blu ;
		}
	*/
	*COLOR = (color1 * noiseCoef) + (color2 * (1.0f - noiseCoef));
}

void Wood::Init(void)
{

	// center.x	= 0.0 ;
	// center.y	= 0.0 ;
	// center.z	= 0.0 ;

	color1.Set(0.73, 0.46, 0.25);
	color2.Set(0.29, 0.19, 0.10);

	scale.Set(0.01f, 0.01f, 0.01f);

	grain_scale = 5.0f;
}

void Marble::Init(void)
{

	// center.x	= 0.0 ;
	// center.y	= 0.0 ;
	// center.z	= 0.0 ;

	color1.x = 0.0;
	color1.y = 0.0;
	color1.z = 1.0;

	color2.x = 1.0;
	color2.y = 1.0;
	color2.z = 0.0;

	scale.x = 0.01f;
	scale.y = 0.01f;
	scale.z = 0.01f;
}

void Marble::SetScale(Vector3d SCALE)
{
	scale = SCALE;
}

void Marble::SetColor(Vector3d COLOR1, Vector3d COLOR2)
{
	color1 = COLOR1;
	color2 = COLOR2;
}

void Marble::GetColor(Vector3d point, Vector3d *COLOR)
{
	float noiseCoef;
	Vector3d temp;

	noiseCoef = 0.0f;
	for (int level = 1; level < 10; level++)
	{
		noiseCoef += (1.0f / level) * fabsf(float(noise(level * 0.05 * point.x,
														level * 0.05 * point.y,
														level * 0.05 * point.z)));
	};

	noiseCoef = 0.5f * sinf((point.x + point.y) * 0.05f + noiseCoef) + 0.5f;
	/*
		noiseCoef	= 0.0f ;
		for (int level = 1; level < 10; level ++)
		{
			noiseCoef +=  (1.0f / level)
					* fabsf ( float( noise ( level * scale.x * point.x,
											 level * scale.y * point.y,
											 level * scale.z * point.z
										   )
								  )
							) ;
		};

		noiseCoef = 0.5f * sinf( (point.x + point.y) * scale.x + noiseCoef) + 0.5f ;

		noiseCoef	= 0.0f ;
		for (int level = 1; level < 4; level ++)
		{
			noiseCoef +=  (1.0f / level)
					* fabsf ( float( noise ( level * scale.x * point.x,
											 level * scale.y * point.y,
											 level * scale.z * point.z
										   )
								  )
							) ;
		} ;

		noiseCoef = 0.5f * sinf( 0.01f * noiseCoef + point.z ) + 0.5f ;
	*/
	float perturb;

	noiseCoef = 0.0f;
	perturb = 0.0f;
	for (int level = 1; level < 3; level++)
	{
		perturb += (1.0f / level) * fabsf(float(noise(level * 0.01 * point.x,
													  level * 0.01 * point.y,
													  level * 0.01 * point.z)));
	};

	int octaves = 7;
	float s = 1.0f;
	Vector3d tp;

	perturb = 0.0f;
	s = 1.0;
	tp = point * scale;

	while (octaves--)
	{
		perturb += fabsf(float(noise(tp.x, tp.y, tp.z))) * s;
		s *= 0.5;
		tp.x *= 2.;
		tp.y *= 2.;
		tp.z *= 2.;
	}
	//	perturb	= Chaos ( point, 7 ) ;
	noiseCoef = 0.5f * sinf(point.x * 0.08f + perturb * 5.0f) + 0.5f;

	//	noiseCoef = 0.5f * sinf( (point.x + point.y) * 0.05f + perturb) + 0.5f ;
	//	Marmol
	//	noiseCoef = 0.5f * sinf( perturb * 2.0f ) + 0.5f ;

	// perturb	= fabsf( float( noise(	0.105f * point.x,
	//								0.105f * point.y,
	//								0.105f * point.z)));
	//	noiseCoef = 0.5f * cos ( point.x * 0.02f + point.y * 0.05f + perturb * 7.0f ) + 0.5f ;

	/*
	{

		Vector3d	v ;
		float		parameter[12] ;
			v.Set(0,1,0);

			// normalized vector in the direction of the cylinder texture
			parameter[0] = (float)v.x;
			parameter[1] = (float)v.y;
			parameter[2] = (float)v.z;

			// point on the axis of cylinder
			parameter[3] = 20;
			parameter[4] = 0;
			parameter[5] = 5;

			// frequency and noise multipliers
			parameter[6] = 4;//3;
			parameter[7] = 1;//3;

			// exponent
			parameter[8] = 8;//3;

			// color multipliers
			parameter[9]  = 1.0; //1;
			parameter[10] = 0.7; //1;
			parameter[11] = 0.2; //1;

		Vector3d	color_result;
		Vector3d	vec ;
		float amp, dist ;
		Vector3d	p, n, o ;

		vec		= point ;

		// get vector in the direction of the trunk
		n.Set(parameter[0], parameter[1], parameter[2]);

		// get the origin of the trunk
		o.Set(parameter[3], parameter[4], parameter[5]);

		// calculate vector from center of tree to p
		p = p-o;

		// get shortest vector to center of trunk
		p = p - n*(p*n);

		dist = (float)p.Magnitude();
	float scaler = 0.02f ;
		amp = (1 + cosf ( parameter[6]*dist + parameter[7]* noise(vec.x*scaler, vec.y*scaler, vec.z*scaler) ) ) / 2 ;

		amp = (float)pow(amp, (int)parameter[8]);
		amp = 1-amp;

		color_result.x = amp*parameter[9];
		color_result.y = amp*parameter[10];
		color_result.z = amp*parameter[11];

			*COLOR	= color_result ;
	}
	*/

	*COLOR = (color1 * noiseCoef) + (color2 * (1.0f - noiseCoef));
}

void Turbulence::Init(void)
{

	// center.x	= 0.0 ;
	// center.y	= 0.0 ;
	// center.z	= 0.0 ;

	color1.x = 0.0;
	color1.y = 0.0;
	color1.z = 1.0;

	color2.x = 1.0;
	color2.y = 1.0;
	color2.z = 0.0;

	scale.x = 200.00;
	scale.y = 200.00;
	scale.z = 200.00;
}

void Turbulence::SetScale(Vector3d SCALE)
{
	scale = SCALE;
}

void Turbulence::SetColor(Vector3d COLOR1, Vector3d COLOR2)
{
	color1 = COLOR1;
	color2 = COLOR2;
}

void Turbulence::GetColor(Vector3d point, Vector3d *COLOR)
{
	float noiseCoef = 0.0f;
	Vector3d temp;

	for (int level = 1; level < 10; level++)
	{
		noiseCoef += (1.0f / level) * fabsf(float(noise(level * 0.05 * point.x,
														level * 0.05 * point.y,
														level * 0.05 * point.z)));
	};

	*COLOR = (color1 * noiseCoef) + (color2 * (1.0f - noiseCoef));
}

/*
double		s, t;
Vector3d	tp;
int			octaves	= 10 ;

	s = 1.0;
	t = 0.0;
	tp = point ;

	while (octaves--)
	{
		t += fabsf(float(noise ( tp.x, tp.y, tp.z ))) * s ;
		s *= 0.5;
		tp.x *= 2.;
		tp.y *= 2.;
		tp.z *= 2.;
	}

	noiseCoef	= sin (1.0f * t + 1.0f * point.z ) ;
	noiseCoef	= fabs(noiseCoef);
*/
/*
	{
	float red, grn, blu;
	float chaos, brownLayer, greenLayer;
	float perturb, brownPerturb, greenPerturb, grnPerturb;
	float t;

//	chaos = Chaos(pos, 7);
	{
double		s, t;
Vector3d	tp;
int			octaves	= 7 ;

	s = 1.0;
	t = 0.0;
	tp = point * scale ;

	while (octaves--)
	{
		t += fabsf(float(noise ( tp.x, tp.y, tp.z ))) * s ;
		s *= 0.5;
		tp.x *= 2.;
		tp.y *= 2.;
		tp.z *= 2.;
	}
	}

	t = sin(sin(8.*chaos + 7*point.x +3.*point.y));

	greenLayer = brownLayer = fabs(t);

	perturb = sin(40.*chaos + 50.*point.z);
	perturb = fabs(perturb);

	brownPerturb = .6*perturb + 0.3;
	greenPerturb = .2*perturb + 0.8;
	grnPerturb = .15*perturb + 0.85;
	grn = 0.5 * pow(fabs(brownLayer), 0.3f);
	brownLayer = pow(0.5 * (brownLayer+1.0), 0.6) * brownPerturb;
	greenLayer = pow(0.5 * (greenLayer+1.0), 0.6) * greenPerturb;

	red = (0.5*brownLayer + 0.35*greenLayer)*2.*grn;
	blu = (0.25*brownLayer + 0.35*greenLayer)*2.0*grn;
//	grn *= max(brownLayer, greenLayer) * grnPerturb;
	if (brownLayer > greenLayer)
		grn *= brownLayer * grnPerturb;
	else
		grn *= greenLayer
		* grnPerturb;

	//surf->diff.r *= red;
	//surf->diff.g *= grn;
	//surf->diff.b *= blu;
	//surf->amb.r *= red;
	//surf->amb.g *= grn;
	//surf->amb.b *= blu;

	COLOR->x	= red ;
	COLOR->y	= grn ;
	COLOR->z	= blu ;
	}
*/
