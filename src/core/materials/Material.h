
#pragma once

#include "core/materials/Bump.h"
#include "math/Vector3d.h"

typedef enum
{
	FLAT = 0,
	MOSAIC,
	CHECKER,
	BRICK,
	SPHERIC,
	LEOPARD,
	GRADIENT,
	WOOD
} TEXTURES;

class Material
{
public:
	double				 // kd,     //  Coef. Difusse
						 //					ks,     //  Coef. Specular
		ES,				 //  Specular Strenght
						 //					kt,     //  Coef. Transmitted
		kl,				 //  Coef. Luminosity
		ir;				 //  Index of Refraction
	unsigned char METAL; //  PHONG / METAL Flag
	Vector3d color;		 // Ambient Color
	Vector3d difusse;	 //	Difusse Color
	Vector3d specular;	 //	Specular Color
	Vector3d reflected;
	Vector3d refracted;
	//	TEXTURES		text_type;
	//	void			*text_data;
	BUMPS bump_type;
	void *bump_data;

public:
	virtual ~Material(void){};
	virtual void Init(void){};
	virtual void GetColor(Vector3d POINT, Vector3d *pColor){};
};

class Flat : public Material
{
public:
	void Init(void){};
	void GetColor(Vector3d POINT, Vector3d *pColor);
};

class Mosaic : public Material
{
public:
	Vector3d center;
	Vector3d color1;
	Vector3d color2;
	Vector3d scale;
	double mortar;

public:
	void Init(void);

	void SetColor(Vector3d COLOR1, Vector3d COLOR2);
	void GetColor(Vector3d POINT, Vector3d *color);
	void SetCenter(Vector3d CENTER);
	void SetScale(Vector3d SCALE);
	void SetMortar(double MORTAR);
};

class Checker : public Material
{
public:
	Vector3d center;
	Vector3d color1;
	Vector3d color2;
	Vector3d scale;

public:
	void Init(void);

	void GetColor(Vector3d POINT, Vector3d *color);
	void SetScale(Vector3d scale);
	void SetColor(Vector3d color1, Vector3d color2);
};

class Brick : public Material
{
public:
	Vector3d center;
	Vector3d color1;
	Vector3d color2;
	Vector3d scale;
	double mortar;

public:
	void Init(void);

	void GetColor(Vector3d POINT, Vector3d *color);
	void SetCenter(Vector3d CENTER);
	void SetColor(Vector3d COLOR1, Vector3d COLOR2);
	void SetScale(Vector3d SCALE);
	void SetMortar(double MORTAR);
};

class Spheric : public Material
{
public:
	Vector3d center;
	double radius1;
	double radius2;
	Vector3d color1;
	Vector3d color2;
	Vector3d scale;
	double blend;

public:
	void Init(void);

	void SetCenter(Vector3d CENTER);
	void SetColor(Vector3d COLOR1, Vector3d COLOR2);
	void SetRadius(double RADIUS1, double RADIUS2);
	void SetScale(Vector3d SCALE);
	void GetColor(Vector3d POINT, Vector3d *COLOR);
	void SetBlend(double BLEND);
};

class Leopard : public Material
{
public:
	Vector3d center;
	Vector3d color1;
	Vector3d color2;
	Vector3d scale;
	double blend;

public:
	void Init(void);

	void SetCenter(Vector3d CENTER);
	void SetColor(Vector3d COLOR1, Vector3d COLOR2);
	void SetScale(Vector3d SCALE);
	void SetBlend(double BLEND);
	void GetColor(Vector3d point, Vector3d *COLOR);
};

class Gradient : public Material
{
public:
	Vector3d center;
	Vector3d scale;
	Vector3d gradient;
	Vector3d color1;
	Vector3d color2;
	double blend;

public:
	void Init(void);

	void SetGradient(Vector3d GRADIENT);
	void SetCenter(Vector3d CENTER);
	void SetColor(Vector3d COLOR1, Vector3d COLOR2);
	void SetScale(Vector3d SCALE);
	void SetBlend(double BLEND);
	void GetColor(Vector3d point, Vector3d *COLOR);
};

class Marble : public Material
{
public:
	Vector3d scale;
	Vector3d color1;
	Vector3d color2;

public:
	void Init(void);

	void SetColor(Vector3d COLOR1, Vector3d COLOR2);
	void SetScale(Vector3d SCALE);
	void GetColor(Vector3d point, Vector3d *COLOR);
};

class Turbulence : public Material
{
public:
	Vector3d scale;
	Vector3d color1;
	Vector3d color2;

public:
	void Init(void);

	void SetColor(Vector3d COLOR1, Vector3d COLOR2);
	void SetScale(Vector3d SCALE);
	void GetColor(Vector3d point, Vector3d *COLOR);
};

class Wood : public Material
{
public:
	Vector3d scale;
	Vector3d color1;
	Vector3d color2;
	float grain_scale;

public:
	void Init(void);

	void SetColor(Vector3d COLOR1, Vector3d COLOR2);
	void SetScale(Vector3d SCALE);
	void GetColor(Vector3d point, Vector3d *COLOR);

	void SetGrainScale(float pGrainScale) { grain_scale = pGrainScale; };
	float GetGrainScale(void) { return grain_scale; };
};

/*
//	Atributos del material
Ke			: emission ;
Ka			: ambient ;
Kd			: diffuse ;
Ks			: specular ;
shininess	: shininess ;

// Compute the emissive term
emissive = Ke;

// Compute the ambient term
ambient = Ka * globalAmbient;

// Compute the diffuse term
LightDirection	= normalize ( IN.lightDirection ) ;
//LightDirection	= normalize ( lightPosition - eyePosition ) ;
//LightDirection	= normalize ( lightPosition - P ) ;
diffuseLight	= max ( dot ( N, LightDirection ), 0.0 ) ;
diffuse		= Kd * lightColor * diffuseLight ;


// Compute the specular term
ViewerDirection	= normalize( eyePosition.xyz - P );
ViewerDirection	= normalize( - P );
//ViewerDirection	= normalize( float3( 0.0, 0.0, 0.0 ) - P );
H		= normalize( LightDirection + ViewerDirection );
specularLight	= pow(max(dot(N, H), 0), shininess);

if (diffuseLight <= 0) specularLight = 0;
specular	= Ks * lightColor * specularLight;

diffuse		= Kd * lightColor * diffuseLight ;

color.xyz	= emissive + ambient + diffuse + specular ;
color.w		= 1 ;

*/

struct check_type
{
	Vector3d color;
	short object_hit;
	double distance;
	Material *material;
	Vector3d normal;
};
