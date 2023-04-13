
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "core/geometry/Geometry.h"

long GeometryManager::Add(Geometry *pObject)
{
	/*
	object_type	*tobject ;



	//	tobject						= (object_type *)malloc ( sizeof ( object_type ) ) ;
		tobject						= new ( object_type ) ;

		if (tobject==NULL)
			return ( -1 ) ;

		tobject->type_of_object		= object->type_of_object ;
		tobject->obj.sphere			= object->obj.sphere ;
		tobject->obj.plane			= object->obj.plane ;
		tobject->material			= object->material ;
		tobject->next				= NULL ;

		if (Lista->first == NULL)
		{
			Lista->first			= tobject ;
		}
		else
		{
			Lista->last->next	= tobject ;
		}

		Lista->last				= tobject ;
		Lista->actual			= tobject ;
	*/
	objects.push_back(pObject);

	return (0);
}
