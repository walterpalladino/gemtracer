
#pragma once

#include <string>
#include <vector>

using namespace std;

#include "math/Math.h"
#include "core/materials/Material.h"

#define OBJECT_GENERIC 0
#define RECTANGLE 1
#define INFINITE_PLANE 2
#define SPHERE 3
#define TRIANGLE_PATCH 4

class Geometry
{
public:
	int type_of_object;
	char name[64];
	//	Vector3d		color;
	Material *material;

public:
	//	Object ( void ) {} ;
	virtual ~Geometry(void){};
	void SetMaterial(Material *pMaterial) { material = pMaterial; };
	Material *GetMaterial(void) { return material; };

	virtual short WhoIAm(void) { return (OBJECT_GENERIC); };
	virtual void GetNormal(Vector3d Point, Vector3d *normal){};
	virtual double Intersect(Vector3d rayvec, Vector3d raystart) { return 0.0; };
};

/*
typedef struct sphere_type {
  double l,m,n; // Center of sphere
  double r;     // Radius of sphere
} TSPHERE ;
*/

class GeometryManager
{
public:
	vector<Geometry *> objects;

public:
	static GeometryManager *GetInstance()
	{
		static GeometryManager instance;
		return &instance;
	};

	long Add(Geometry *pObject);
};

// #define	NUMBER_OF_OBJECTS	5
