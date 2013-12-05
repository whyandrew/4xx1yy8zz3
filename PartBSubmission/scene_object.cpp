/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements scene_object.h

***********************************************************/

#include <cmath>
#include <iostream>
#include "scene_object.h"

bool SceneObject::solveT(double A, double B, double C, double* p_tValue)
{
	// Solve and store correct value of t in quadratic equation
	bool b_soln = false;
	double t_value = 0.0, t_value2;
	double D = (B * B) - (4 * A * C);

	if (D == 0.0 && A != 0.0)
	{
		// 1 intersect. 
		b_soln = true;
		t_value = - B / (2 * A);
	}
	else if (D > 0.0 && A != 0.0)
	{
		// 2 intersects, have to pick correct one, t_value > t_value2
		t_value = - (B / (2*A) ) + (sqrt(D) / (2*A));
		t_value2 = - (B / (2*A) ) - (sqrt(D) / (2*A));
		if (t_value > 0.0 && t_value2 < 0.0)
		{
			b_soln = true;
		}
		else if (t_value2 > 0.0)
		{
			t_value = t_value2;
			b_soln = true;
		}
	}
	else if (A == 0.0)
	{
		t_value = -C / B;
		b_soln = true;
	}

	*p_tValue = t_value;
	return b_soln;
}

bool CompoundObject::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld, bool b_shadowRay ) 
{
	// Loop through p_objList to find the closest intersection with ray

	bool b_isHit = false;
	SceneObject *p_obj;
	SceneObject **objPList = p_objPList;

	for (int i = 0; i < numObj; i++)
	{
		p_obj = *(objPList + i);
		b_isHit |= p_obj->intersect(ray, worldToModel, modelToWorld, false);

	}

	return b_isHit;
}

bool UnitSquare::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld, bool b_shadowRay ) {
	// TODO: implement intersection code for UnitSquare, which is
	// defined on the xy-plane, with vertices (0.5, 0.5, 0), 
	// (-0.5, 0.5, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
	// (0, 0, 1).
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.

	bool b_isHit = false; // if ray actually hit this first
	Vector3D ray_dir = worldToModel * ray.dir;
	Point3D ray_orig = worldToModel * ray.origin;
	double t_value;

	// Sqaure is on xy-plane, so z==0
	// If ray parellel to xy-plane, then either it misses the square
	// or it technically intersect the "thickness" of the 2D square, which will be ingored
	if (ray_dir[2] != 0)
	{
		// Find ray intersection point to xy-plane
		t_value = -ray_orig[2] / ray_dir[2];
		if (t_value > 0.0) 
		{
			Point3D intersection = ray_orig + t_value * ray_dir;
			if (abs(intersection[0]) <= 0.50 && abs(intersection[1]) <= 0.50)
			{
				// within +/- 0.5 is an intersection
				// Now compare to other intersections, if any
				if (ray.intersection.none || ray.intersection.t_value > t_value)
				{
					b_isHit = true;
					ray.intersection.none = false;
					ray.intersection.t_value = t_value;
					//ray.intersection.fp_textureMapping = &SceneObject::textureMapping;
					ray.intersection.p_sceneObj = this;
					ray.intersection.M_m2w = modelToWorld;
					ray.intersection.M_w2m = worldToModel;

					if (!b_shadowRay)
					{
						Vector3D normal = transNorm(worldToModel, Vector3D(0.0, 0.0, 1.0));
						normal.normalize();
						ray.intersection.normal = normal;
						ray.intersection.point = modelToWorld * intersection;
					}
				}
			}
		}
	}
	
	return b_isHit;
}

bool UnitSphere::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld, bool b_shadowRay ) {
	// TODO: implement intersection code for UnitSphere, which is centred 
	// on the origin.  
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.
	bool b_isHit = false;
	Vector3D ray_dir = worldToModel * ray.dir;
	Point3D ray_orig = worldToModel * ray.origin;
	double t_value, t_value2;

	// equation: point^2 = 1
	// Solve: A*t^2 + 2B*t + C = 0, where
	// A = d.d ; B = (ray_orig - center).d ; 
	//C = (ray_orig - center).(ray_orig - center) - 1
	double d_A = ray_dir.dot(ray_dir);
	Vector3D ray_orig_vect = ray_orig - Point3D(0.0, 0.0, 0.0);
	double d_B = ray_orig_vect.dot(ray_dir);
	double d_C = ray_orig_vect.dot(ray_orig_vect) - 1;
	// Calculate D = B^2 - AC; 
	// t = -B/A +/- sqrt(D)/A
	double d_D = d_B * d_B - d_A * d_C;

	// A should not == 0, unless ray_dir is 0, which should not happen
	if (d_D == 0.0 && d_A != 0.0)
	{
		// 1 intersect. 
		b_isHit = true;
		t_value = - d_B / d_A;
	}
	else if (d_D > 0.0 && d_A != 0.0)
	{
		// 2 intersects, have to pick correct one, t_value > t_value2
		t_value = - (d_B / d_A ) + (sqrt(d_D) / d_A);
		t_value2 = - (d_B / d_A ) - (sqrt(d_D) / d_A);
		if (t_value > 0.0 && t_value2 < 0.0)
		{
			b_isHit = true;
		}
		else if (t_value2 > 0.0)
		{
			t_value = t_value2;
			b_isHit = true;
		}
	}

	// Populate ray.intersection 
	if (b_isHit)
	{
		if (ray.intersection.none || ray.intersection.t_value > t_value)
		{
			Point3D hitPt = ray_orig + (t_value * ray_dir);
			ray.intersection.none = false;
			ray.intersection.t_value = t_value;
			//ray.intersection.fp_textureMapping = &SceneObject::textureMapping;
			ray.intersection.p_sceneObj = this;
			ray.intersection.M_m2w = modelToWorld;
			ray.intersection.M_w2m = worldToModel;
			if (!b_shadowRay)
			{
				ray.intersection.point = modelToWorld * hitPt;
				Vector3D normal = transNorm( worldToModel, (hitPt - Point3D(0.0, 0.0, 0.0)) );
				normal.normalize();
				ray.intersection.normal = normal;
			}
		}
		else
		{
			// Has intersection but not the first one ray hits
			b_isHit = false;
		}
	}

	return b_isHit;
}

bool _Hyperboloid::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld, bool b_shadowRay ) 
{
	// "Draw" a hyperboloid (1 surface) with equation x^2 + y^2 - z^2 = 1
	// Object centered at 0,0,0. , height of z

	bool b_isHit = false;
	Vector3D D = worldToModel * ray.dir;
	Point3D O = worldToModel * ray.origin;
	double t_value;
	double zRange = _zRange>0.0? _zRange: -(_zRange);
	Point3D hitPt;

	// A = D.z*D.z-D.x*D.x-D.y*D.y		O = origin, D = dir
	// B = 2.0*(O.z*D.z - O.x*D.x - O.y*D.y)
	// C = O.z*O.z + k - O.x*O.x - O.y*O.y
	double d_A = (D[0]*D[0]) + (D[1]*D[1]) - (D[2]*D[2]);
	double d_B = 2.0 * ( (O[0]*D[0]) + (O[1]*D[1]) - (O[2]*D[2]) );
	double d_C = (O[0]*O[0]) + (O[1]*O[1]) - (O[2]*O[2]) - 1;

	b_isHit = solveT(d_A, d_B, d_C, &t_value);
	
	// Populate ray.intersection 
	if (b_isHit)
	{
		if (ray.intersection.none || ray.intersection.t_value > t_value)
		{
			hitPt = O + (t_value * D);
		}
		else
		{ 
			b_isHit = false;
		}
	}
	
	if (b_isHit)
	{
		// Restrict hyperboloid to z = [+/- zRange]
		if (hitPt[2] >= -zRange && hitPt[2] <= zRange)
		{
			ray.intersection.none = false;
			ray.intersection.t_value = t_value;
			if (!b_shadowRay)
			{
				ray.intersection.point = modelToWorld * hitPt;
				// f:  P.x^2 + P.y^2 - P.z^2 -1 = 0
				// gradient(f) = ( 2P.x, 2P.y, -2Pz)
				Vector3D normal(2 * hitPt[0], 2 * hitPt[1], -2 * hitPt[2]); 
				normal.normalize();
				//ray.intersection.normal = transNorm( worldToModel, normal);
				normal = transNorm( worldToModel, normal);
				//normal.normalize();
				ray.intersection.normal = normal;
			}
		}
		else 
		{
			b_isHit = false;
		}
	}

	return b_isHit;
}

bool _Circle::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld, bool b_shadowRay ) 
{
	// Draw circular plane parellal to xy plane.

	bool b_isHit = false;
	Vector3D D = worldToModel * ray.dir;
	Point3D O = worldToModel * ray.origin;
	double t_value;
	double radius = _radius>0? _radius: -(_radius);
	Point3D hitPt;
	double zvalue = _zvalue;

	// Find (x,y) value when ray hits z=zvalue
	if (D[2] != 0.0)
	{
		t_value = (zvalue - O[2]) / D[2];
		if (t_value > 0.0)
		{
			b_isHit = true;
			hitPt = O + (t_value * D);
		}
	}

	// See if (x,y) is within circle of radius = _radius
	if (b_isHit)
	{
		if ( (hitPt[0]*hitPt[0] + hitPt[1]*hitPt[1]) > radius*radius)
		{
			b_isHit = false;
		}
	}

	// Populate ray.intersection 
	if (b_isHit)
	{
		if (ray.intersection.none || ray.intersection.t_value > t_value)
		{
			ray.intersection.none = false;
			ray.intersection.t_value = t_value;
			if (!b_shadowRay)
			{
				ray.intersection.point = modelToWorld * hitPt;
				// Normal is facing +z by default, flip if (_flipNormal)
				Vector3D normal = _flipNormal? Vector3D(0.0, 0.0, -1.0): Vector3D(0.0, 0.0, 1.0);
				ray.intersection.normal = transNorm( worldToModel, normal);
				ray.intersection.normal.normalize();
			}
		}
		else
		{
			// Has intersection but not the first one ray hits
			b_isHit = false;
		}
	}

	return b_isHit;
}

void Hyperboloid::construct()
{
	// Construct a close hyperboloid with 2 circular planes at the ends
	// Object centered at 0,0,0. , height of 2*zRange

	double radius;
	double zValue = _zRange>0? _zRange : -(_zRange);
	static SceneObject* p_List[3];
	
	p_List[0] = new _Hyperboloid(_zRange*2.0);
	// radius of circle = max value of x,y = sqrt(1 + z^2)
	radius = sqrt( 1 + (zValue * zValue) );
	p_List[1] = new _Circle(radius, zValue, false);
	p_List[2] = new _Circle(radius, -zValue, true);

	p_objPList = p_List;
}

Colour UnitSquare::textureMapping(const Ray3D& ray)
{
	// Intersection pt (x,y) should be between -0,5 to 0.5
	Point3D pt = ray.intersection.M_w2m * ray.intersection.point;
	unsigned long int picWidth = ray.intersection.mat->txt_width;
	unsigned long int picHeight = ray.intersection.mat->txt_height;
	// Range from 0 to 1.0
	double x = pt[0] + 0.5;
	double y = pt[1] + 0.5;

	// Range over size of texture e.g. 0 to 1024
	unsigned long int xCoord;
	unsigned long int yCoord;
	xCoord = (unsigned long int)(x * picWidth);
	yCoord = (unsigned long int)(y * picHeight);
	
	unsigned long int pixel = (xCoord + yCoord * picWidth);
	Material *mat = ray.intersection.mat;
	Colour pixelCol = Colour(mat->txt_rbuffer[pixel] / 255.0,
		mat->txt_gbuffer[pixel] / 255.0,
		mat->txt_bbuffer[pixel] / 255.0);

	//ray.col = ray.col + pixelCol;
	return pixelCol;
}

Colour UnitSphere::textureMapping(const Ray3D& ray)
{
	//printf("UnitShpere Texture\n");

	// Map texture to sphere as vertical and horizontal components
	// Ray must contain intersection.point 

	Vector3D vert(0.0, 1.0, 0.0);
	Vector3D horiz(1.0, 0.0, 0.0);
	Point3D pt = ray.intersection.M_w2m * ray.intersection.point;
	Vector3D ptVect(pt[0], pt[1], pt[2]);
	//ptVect.normalize();
	double x, y;

	//Angle between pt and vertical
	double phi = acos( -ptVect.dot(vert) );
	y = phi / M_PI;

	// angle to horizontal
	double theta = ( acos( ptVect.dot(horiz)/ sin(phi)) )/(2 * M_PI);
	if ( horiz.dot( vert.cross(horiz) ) > 0)
	{
		x = theta;
	}
	else
	{
		x = 1 - theta;
	}

	//printf("(%f, %f)  ", x, y);
	// Now get the corresponding pixel from texture
	Material *p_mat = ray.intersection.mat;
	unsigned long int xCoord;
	unsigned long int yCoord;
	xCoord = (unsigned long int)(x * p_mat->txt_width);
	yCoord = (unsigned long int)(y * p_mat->txt_height);

	unsigned long int pixel = xCoord + yCoord * p_mat->txt_width;

	Colour pixelCol = Colour((double)p_mat->txt_rbuffer[pixel] / 255.0 ,
		(double)p_mat->txt_gbuffer[pixel] / 255.0, 
		(double)p_mat->txt_bbuffer[pixel] / 255.0);
	
	//ray.col = ray.col + pixelCol;
	//ray.col.clamp();
	return pixelCol;
}

Colour _Hyperboloid::textureMapping(const Ray3D& ray)
{
	Colour pixelColor;
	return pixelColor;
}

Colour _Circle::textureMapping(const Ray3D& ray)
{
	Colour pixelColor;
	return pixelColor;
}

Colour Hyperboloid::textureMapping(const Ray3D& ray)
{
	Colour pixelColor;
	return pixelColor;
}