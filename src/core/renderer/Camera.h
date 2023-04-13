
#pragma once

#include    <stdio.h>
#include    <stdlib.h>
#include    <math.h>

#include	"math/Vector3d.h"


class Camera
{
public:
	long			id ;
	long			status ;
	double			angle ;

	Vector3d			location ;
	Vector3d			direction ;
	Vector3d			right ;
	Vector3d			up ;

} ;

