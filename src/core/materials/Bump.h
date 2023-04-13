

#pragma once

#include "math/Vector3d.h"

typedef enum
{
	NO_BUMP = 0,
	WAVE
} BUMPS;

typedef struct wave_type
{
	Vector3d center;	//  Center of the Wave
	double wavelength;	//  Wavelenght
	double phase;		//  Startpoint
	double attenuation; //  Attenuation
	double amplitude;	//  Amplitude

} TWAVE;

void get_Bump_Normal(TWAVE *bump, Vector3d POINT, Vector3d *R);
