
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "core/renderer/View.h"
#include "core/renderer/Ray.h"
#include "core/renderer/RayTracer.h"

#include "core/geometry/Geometry.h"
#include "core/geometry/InfinitePlane.h"
#include "core/geometry/Sphere.h"
#include "core/geometry/TrianglePatch.h"

#include "core/lights/LightManager.h"

#include "core/materials/Material.h"
#include "core/materials/bump.h"

#include "math/Vector3d.h"

long NUMBER_OF_LIGHTS;

Vector3d TColor;

Vector3d AMBIENT_INTENSITY(0.3, 0.3, 0.3);
Vector3d BACKGROUND_COLOR(0.65, 0.65, 0.65);

#ifndef M_PI_360
#define M_PI_360 0.00872664625997164788
#endif
Vector3d VZERO(0.0, 0.0, 0.0);

View::View(void)
{
	width = 1920;
	height = 1080;
}

long View::Init(void)
{

	//	Objects
	Vector3d scale;
	Vector3d color1,
		color2;
	Vector3d center;

	Patch *new_patch;
	Vector3d vertex[3];

	Light light;

	InfinitePlane *new_plane;
	Sphere *new_sphere;

	/*
	Sphere	sphere ( Vector3d ( 500.0, 300.0, 1000.0 ),   // Coordinates of center
						200.0
					) ;       // Radius
	Sphere	sphere2 ( Vector3d ( -400.0, 500.0, 800.0 ),   // Coordinates of center
							170.0 ) ;       // Radius
	Sphere	sphere3 ( Vector3d ( 0.0, 400.0 ,700.0 ),   // Coordinates of center
							100.0 ) ;       // Radius
	*/
	/*
	TSPHERE sphere4={-50,100,800,   // Coordinates of center
							150};       // Radius
	*/
	Sphere sphere4(Vector3d(-50.0, -25.0, 1500.0), // Coordinates of center
				   150.0);						   // Radius
												   /*
												   Plane	plane( Vector3d (0.0, 1.0, 0.0), 	 // Surface normal
																	   300			   // Distance from origin
																	   ) ;
												   */

	//	Video Buffer
	imageBuffer.Create(width, height);

	///////////////////////////////////////////
	//	20091216

	TWAVE *onda01;

	onda01 = new (TWAVE);
	onda01->center = Vector3d(500.0, 0.0, 2000.0); //  Center of the Wave
	onda01->wavelength = 500;					   //  Wavelenght
	onda01->phase = 0;							   // Startpoint
	onda01->attenuation = 0.10;					   //  Attenuation
	onda01->amplitude = 10.0;					   //  Amplitude

	Leopard *material_pintas;
	Brick *material_wall;
	Spheric *material_stripes;
	Checker *material_floor;
	Flat *material_flat_mirror;
	Flat *material06;
	Marble *material_marble;
	Wood *material_wood;

	material_pintas = new (Leopard);

	material_pintas->Init();

	material_pintas->difusse = Vector3d(0.5, 0.5, 0.5);
	material_pintas->specular = Vector3d(0.5, 0.5, 0.5);
	material_pintas->ES = 100.0;							 //  Specular Strenght
	material_pintas->kl = 0.0;								 //  Coef. Luminosity
	material_pintas->ir = 0.10;								 //  Index of Refraction
	material_pintas->METAL = 0;								 //  PHONG / METAL Flag
	material_pintas->color = Vector3d(0.80, 0.80, 0.80);	 // color
	material_pintas->reflected = Vector3d(0.1, 0.1, 0.1);	 // reflected
	material_pintas->refracted = Vector3d(0.40, 0.40, 0.50); // refracted
	material_pintas->bump_type = NO_BUMP;
	material_pintas->bump_data = NULL;

	material_pintas->SetBlend(1.0);
	center.Set(100.0, -150.0, 600.0);
	material_pintas->SetCenter(center);
	color1.Set(0.5, 0.5, 1.0);
	color2.Set(1.0, 1.0, 1.0);
	material_pintas->SetColor(color1, color2);
	scale.Set(10.0, 10.0, 10.0);
	material_pintas->SetScale(scale);

	material_wall = new (Brick);

	material_wall->Init();

	material_wall->difusse = Vector3d(1.0, 1.0, 1.0);
	material_wall->specular = Vector3d(0.0, 0.0, 0.0);
	material_wall->ES = 0.0;							//  Specular Strenght
	material_wall->kl = 0.0;							//  Coef. Luminosity
	material_wall->ir = 0.0;							//  Index of Refraction
	material_wall->METAL = 0;							//  PHONG / METAL Flag
	material_wall->color = Vector3d(0.60, 0.80, 0.80);	// color
	material_wall->reflected = Vector3d(0.0, 0.0, 0.0); // reflected
	material_wall->refracted = Vector3d(0.0, 0.0, 0.0); // refracted
	material_wall->bump_type = NO_BUMP;
	material_wall->bump_data = NULL;

	material_stripes = new (Spheric);

	material_stripes->Init();

	material_stripes->difusse = Vector3d(0.8, 0.8, 0.8);
	material_stripes->specular = Vector3d(0.4, 0.4, 0.4);
	material_stripes->ES = 100.0;							  //  Specular Strenght
	material_stripes->kl = 0.0;								  //  Coef. Luminosity
	material_stripes->ir = 1.30;							  //  Index of Refraction
	material_stripes->METAL = 0;							  //  PHONG / METAL Flag
	material_stripes->color = Vector3d(0.60, 0.40, 0.40);	  // color
	material_stripes->reflected = Vector3d(0.0, 0.0, 0.0);	  // reflected
	material_stripes->refracted = Vector3d(0.20, 0.20, 0.20); // refracted
	material_stripes->bump_type = NO_BUMP;
	material_stripes->bump_data = NULL;

	color1.Set(1.0, 1.0, 1.0);
	color2.Set(1.0, 0.0, 0.0);
	material_stripes->SetColor(color1, color2);
	material_stripes->SetRadius(10.0, 40.0);
	material_stripes->SetBlend(0.3);
	center.Set(0, 0, 0);
	material_stripes->SetCenter(center);
	scale.Set(1.0, 0.0, 0.0);
	material_stripes->SetScale(scale);

	material_floor = new (Checker);

	material_floor->Init();

	material_floor->difusse = Vector3d(0.8, 0.8, 0.8);
	material_floor->specular = Vector3d(0.6, 0.6, 0.6);
	material_floor->ES = 50.0;							 //  Specular Strenght
	material_floor->kl = 0.0;							 //  Coef. Luminosity
	material_floor->ir = 0.0;							 //  Index of Refraction
	material_floor->METAL = 0;							 //  PHONG / METAL Flag
	material_floor->color = Vector3d(0.30, 0.30, 0.30);	 // color
	material_floor->reflected = Vector3d(0.0, 0.0, 0.0); // reflected
	material_floor->refracted = Vector3d(0.0, 0.0, 0.0); // refracted
	// material_floor->bump_type	= WAVE ;
	// material_floor->bump_data	= onda01 ;
	material_floor->bump_type = NO_BUMP;
	material_floor->bump_data = NULL;

	scale.Set(100.0, 100.0, 100.0);
	material_floor->SetScale(scale);
	color1.Set(1.0, 1.0, 1.0);
	color2.Set(0.2, 0.2, 0.6);
	material_floor->SetColor(color1, color2);

	material_flat_mirror = new (Flat);

	material_flat_mirror->Init();

	material_flat_mirror->difusse = Vector3d(0.60, 0.60, 0.60);
	material_flat_mirror->specular = Vector3d(1.0, 1.0, 1.0);
	material_flat_mirror->ES = 10.0;							  //  Specular Strenght
	material_flat_mirror->kl = 0.0;								  //  Coef. Luminosity
	material_flat_mirror->ir = 0.0;								  //  Index of Refraction
	material_flat_mirror->METAL = 1;							  //  PHONG / METAL Flag
	material_flat_mirror->color = Vector3d(0.40, 0.40, 0.40);	  // color
	material_flat_mirror->reflected = Vector3d(0.60, 0.60, 0.60); // reflected
	material_flat_mirror->refracted = Vector3d(0.0, 0.0, 0.0);	  // refracted
	// material_flat_mirror->text_type	= FLAT ;					//	text_type
	// material_flat_mirror->text_data	= NULL ;					//	text_data
	material_flat_mirror->bump_type = NO_BUMP;
	material_flat_mirror->bump_data = NULL;

	material06 = new (Flat);

	// material06->kd	= 0.50 ;			     	//  Coef. Difusse
	material06->difusse = Vector3d(0.5, 0.5, 0.5);
	// material06->ks	= 0.20 ;    				//  Coef. Specular
	material06->specular = Vector3d(0.2, 0.2, 0.2);
	material06->ES = 100.0; //  Specular Strenght
	// material06->kt	= 0.5 ; 				    //  Coef. Transmitted
	material06->kl = 0.0;								//  Coef. Luminosity
	material06->ir = 1.5;								//  Index of Refraction
	material06->METAL = 0;								//  PHONG / METAL Flag
	material06->color = Vector3d(0.70, 0.70, 0.70);		// color
	material06->reflected = Vector3d(0.10, 0.10, 0.10); // reflected
	material06->refracted = Vector3d(0.70, 0.70, 0.70); // refracted
	// material06->text_type	= FLAT ;					//	text_type
	// material06->text_data	= NULL ;					//	text_data
	material06->bump_type = NO_BUMP;
	material06->bump_data = NULL;

	material_marble = new (Marble);

	material_marble->Init();

	material_marble->difusse = Vector3d(0.5, 0.5, 0.5);
	material_marble->specular = Vector3d(0.5, 0.5, 0.5);
	material_marble->ES = 100.0;						  //  Specular Strenght
	material_marble->kl = 0.0;							  //  Coef. Luminosity
	material_marble->ir = 0.0;							  //  Index of Refraction
	material_marble->METAL = 0;							  //  PHONG / METAL Flag
	material_marble->color = Vector3d(0.80, 0.80, 0.80);  // color
	material_marble->reflected = Vector3d(0.0, 0.0, 0.0); // reflected
	material_marble->refracted = Vector3d(0.0, 0.0, 0.0); // refracted
	material_marble->bump_type = NO_BUMP;
	material_marble->bump_data = NULL;

	// material_marble->SetBlend ( 1.0 ) ;
	center.Set(100.0, -150.0, 600.0);
	// material_marble->SetCenter ( center ) ;
	color1.Set(0.9, 0.9, 0.9);
	// color2.Set ( 1.0, 1.0, 1.0 ) ;
	color2.Set(0.0, 0.0, 0.0);
	material_marble->SetColor(color1, color2);
	scale.Set(0.005, 0.02, 0.05);
	material_marble->SetScale(scale);

	material_wood = new (Wood);

	material_wood->Init();

	material_wood->difusse = Vector3d(0.5, 0.5, 0.5);
	material_wood->specular = Vector3d(0.5, 0.5, 0.5);
	material_wood->ES = 5.0;							//  Specular Strenght
	material_wood->kl = 0.0;							//  Coef. Luminosity
	material_wood->ir = 0.0;							//  Index of Refraction
	material_wood->METAL = 0;							//  PHONG / METAL Flag
	material_wood->color = Vector3d(0.80, 0.80, 0.80);	// color
	material_wood->reflected = Vector3d(0.0, 0.0, 0.0); // reflected
	material_wood->refracted = Vector3d(0.0, 0.0, 0.0); // refracted
	material_wood->bump_type = NO_BUMP;
	material_wood->bump_data = NULL;

	// material_marble->SetBlend ( 1.0 ) ;
	center.Set(100.0, -150.0, 600.0);
	// material_marble->SetCenter ( center ) ;
	color1.Set(0.9, 0.9, 0.9);
	// color2.Set ( 1.0, 1.0, 1.0 ) ;
	color2.Set(0.0, 0.0, 0.0);
	// material_wood->SetColor ( color1, color2 ) ;
	scale.Set(0.005, 0.02, 0.05);
	// material_wood->SetScale ( scale ) ;

	///////////////////////////////////////////

	new_sphere = new (Sphere);
	new_sphere->Init(Vector3d(200.0, 200.0, 500.0), // Coordinates of center
					 150.0);						// Radius
													//	new_sphere->SetMaterial ( material_pintas ) ;
	new_sphere->SetMaterial(material_wood);
	GeometryManager::GetInstance()->Add(new_sphere);

	new_sphere = new (Sphere);
	new_sphere->Init(Vector3d(-300.0, 150.0, 800.0), // Coordinates of center
					 170.0);						 // Radius
	new_sphere->SetMaterial(material_stripes);
	GeometryManager::GetInstance()->Add(new_sphere);

	new_sphere = new (Sphere);
	new_sphere->Init(Vector3d(0.0, 400.0, 900.0), // Coordinates of center
					 200.0);					  // Radius
	new_sphere->SetMaterial(material_flat_mirror);
	GeometryManager::GetInstance()->Add(new_sphere);

	new_plane = new (InfinitePlane);
	new_plane->Init(Vector3d(0.0, 1.0, 0.0), // Surface normal
					150						 // Distance from origin
	);
	new_plane->SetMaterial(material_floor);
	GeometryManager::GetInstance()->Add(new_plane);

	/*
	 */

	// new_patch	= (TPATCH*)malloc ( sizeof (TPATCH) ) ;
	new_patch = new (Patch);
	new_patch->Init();

	/*
	setV( &vertex[0],  -60.0,    0.0, 475.0 );
	setV( &vertex[1],   60.0,    0.0, 475.0 );
	setV( &vertex[2],  -60.0,  120.0, 475.0 );
	patch.AddTriangle ( vertex[0], vertex[1], vertex[2] );

	setV( &vertex[0],   60.0,    0.0, 475.0 );
	setV( &vertex[1],   60.0,  120.0, 475.0 );
	setV( &vertex[2],  -60.0,  120.0, 475.0 );
	patch.AddTriangle ( vertex[0], vertex[1], vertex[2] );

	setV( &vertex[0],   60.0,  120.0, 475.0 );
	setV( &vertex[1],    0.0,  180.0, 475.0 );
	setV( &vertex[2],  -60.0,  120.0, 475.0 );
	patch.AddTriangle ( vertex[0], vertex[1], vertex[2] );
	*/

	// setV( &vertex[0],  -60.0,    0.0, 800.0 );
	vertex[0].Set(-60.0, 0.0, 800.0);
	// setV( &vertex[1],   60.0,    0.0, 800.0 );
	vertex[1].Set(60.0, 0.0, 800.0);
	// setV( &vertex[2],  -60.0,  120.0, 800.0 );
	vertex[2].Set(-60.0, 120.0, 800.0);
	new_patch->AddTriangle(vertex[0], vertex[1], vertex[2]);

	// setV( &vertex[0],   60.0,    0.0, 800.0 );
	vertex[0].Set(60.0, 0.0, 800.0);
	// setV( &vertex[1],   60.0,  120.0, 800.0 );
	vertex[1].Set(60.0, 120.0, 800.0);
	// setV( &vertex[2],  -60.0,  120.0, 800.0 );
	vertex[2].Set(-60.0, 120.0, 800.0);
	new_patch->AddTriangle(vertex[0], vertex[1], vertex[2]);

	// setV( &vertex[0],   60.0,  120.0, 800.0 );
	vertex[0].Set(60.0, 120.0, 800.0);
	// setV( &vertex[1],    0.0,  180.0, 800.0 );
	vertex[1].Set(0.0, 180.0, 800.0);
	// setV( &vertex[2],  -60.0,  120.0, 800.0 );
	vertex[2].Set(-60.0, 120.0, 800.0);
	new_patch->AddTriangle(vertex[0], vertex[1], vertex[2]);

	new_patch->SetCenter(Vector3d(100.0, 0.0, 0.0));

	new_patch->material = material_wall;

	GeometryManager::GetInstance()->Add(new_patch);

	/*
	strcpy( new_patch.name, "PATCH_1");
	new_brick = new brick;
	new_brick->set_Color( RED, WHITE );
	new_brick->set_Mortar( .15 );
	new_brick->set_Scale( setV( 15.0, 15.0, 15.0 ) );
	new_triangle_patch->set_Texture_Data( new_brick );
	new_triangle_patch->set_Texture( BRICK );
	objects.add_object( new_triangle_patch );
	*/

	/*
	 */

	//	Definicion de Luces
	//	tlight						= ( TLIGHT *)malloc ( sizeof ( TLIGHT ) ) ;

	sprintf(light.name, "LIGHT01");
	//	setV ( &light.position, 500.0,1000.0,-500.0 ) ;
	light.position.Set(500.0, 1000.0, -500.0);
	//	setV ( &light.color, 0.40,0.40,0.40 ) ;
	light.color.Set(0.30, 0.30, 0.30);

	LightManager::GetInstance()->Add(light);

	sprintf(light.name, "LIGHT02");
	light.position.Set(-700.0, 2000.0, -500.0);
	light.color.Set(0.60, 0.60, 0.60);

	LightManager::GetInstance()->Add(light);

	sprintf(light.name, "LIGHT03");
	light.position.Set(700.0, 2000.0, -500.0);
	light.color.Set(0.30, 0.30, 0.30);

	camera.location.Set(0.0, 100.0, -200.0);
	camera.direction.Set(0.0, 0.0, 1.0);
	camera.up.Set(0.0, 1.0, 0.0);
	camera.right.Set(1.33, 0.0, 0.0);
	camera.angle = 70;

	{

		double Right_Length,
			Direction_Length;

		//	normalize( camera.direction, &camera.direction ) ;
		camera.direction.Normalize();
		//	Right_Length	= lenV ( camera.right, VZERO ) ;
		Right_Length = (camera.right - VZERO).Magnitude();
		Direction_Length = Right_Length / tan(camera.angle * M_PI_360) / 2.0;
		//	mulS ( camera.direction, Direction_Length, &camera.direction ) ;
		camera.direction = camera.direction * Direction_Length;
	}

	return (0);
}

void View::RenderScene(float fAngle,
					   float fZPos,
					   float fXPos,
					   float fYPos)
{

	Vector3d rayvec, raystart;

	long yscreen;
	long xscreen;

	check_type check;

	double x0, y0;

	Vector3d color;

	for (yscreen = -1; yscreen <= height; yscreen += 1)
	{
		for (xscreen = -1; xscreen <= width; xscreen += 1)
		{
			raystart.x = camera.location.x;
			raystart.y = camera.location.y;
			raystart.z = camera.location.z;

			color.Set(0.0, 0.0, 0.0);

			float fragmentx = xscreen;
			float fragmenty = yscreen;

			// for (int n = 0; n < 4; n++)

			// for (float fragmentx = xscreen; fragmentx < xscreen + 1.0f; fragmentx += 0.5f)
			//	for (float fragmenty = yscreen; fragmenty < yscreen + 1.0f; fragmenty += 0.5f)
			{

				// fragmentx = xscreen + ((rand() % 200) - 100) / 100000.0f;
				// fragmenty = yscreen + ((rand() % 200) - 100) / 100000.0f;

				//	Convert the x coordinate to be a DBL from -0.5 to 0.5.
				x0 = (double)fragmentx / (double)width - 0.5;
				//	Convert the y coordinate to be a DBL from -0.5 to 0.5.
				y0 = ((double)(height - 1) - (double)fragmenty) / (double)height - 0.5;

				rayvec.x = camera.direction.x +
						   x0 * camera.right.x +
						   y0 * camera.up.x;
				rayvec.y = camera.direction.y +
						   x0 * camera.right.y +
						   y0 * camera.up.y;
				rayvec.z = camera.direction.z +
						   x0 * camera.right.z +
						   y0 * camera.up.z;

				rayvec.Normalize();

				RayTracer::GetInstance()->TraceRay(rayvec, raystart, DEPTH, OUT_OBJECT, &check);

				if (check.object_hit == NO_HIT)
					TColor = RayTracer::GetInstance()->backgroundColor;
				else
					TColor = check.color;

				// color = color + (TColor * 0.25f);
				color = TColor;
			}

			//			if ((yscreen > 0 && yscreen < 601) &&
			//				(xscreen > 0 && xscreen < 801))
			if ((yscreen > 0 && yscreen < height + 1) &&
				(xscreen > 0 && xscreen < width + 1))
			{
				//			Vector3d		AAColor ;

				//				AAColor	= TColor ;
				//				addV ( AAColor, AABuffer [ xscreen ], &AAColor ) ;
				//				addV ( AAColor, AABuffer [ xscreen - 1 ], &AAColor ) ;
				//				mulS ( AAColor, 1.0/3.0, &AAColor ) ;

				imageBuffer.Plot(
					xscreen - 1,
					// height - (yscreen - 1) - 1,
					yscreen - 1,
					(long)(color.x * 255.0),
					(long)(color.y * 255.0),
					(long)(color.z * 255.0),
					255);
			}

			//			AABuffer [ xscreen ]	= TColor ;
		}
	}
}
