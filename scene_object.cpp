/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements scene_object.h

***********************************************************/

#include <cmath>
#include <iostream>
#include "scene_object.h"

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
			Point3D intersection(t_value * ray_dir[0] + ray_orig[0], 
				t_value * ray_dir[1] + ray_orig[1],
				0);
			if (abs(intersection[0]) <= 0.50 && abs(intersection[1]) <= 0.50)
			{
				// within +/- 0.5 is an intersection
				// Now compare to other intersections, if any
				if (ray.intersection.none || ray.intersection.t_value > t_value)
				{
					b_isHit = true;
					ray.intersection.none = false;
					ray.intersection.t_value = t_value;
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
			Point3D intersectPt = ray_orig + (t_value * ray_dir);
			ray.intersection.none = false;
			ray.intersection.t_value = t_value;
			if (!b_shadowRay)
			{
				ray.intersection.point = modelToWorld * intersectPt;
				Vector3D normal = transNorm( worldToModel, (intersectPt - Point3D(0.0, 0.0, 0.0)) );
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

bool hyperboloid::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld, bool b_shadowRay ) 
{
	// "Draw" a hyperboloid (1 surface) with equation x^2 + y^2 - z^2 = 1
	// Top and bottom are enclosed with circle of x^2 + y^2 = 2 @ z = +/- 1
	// Object centered at 0,0,0. , height of 2, and outermost diameter of 2

	bool b_isHit = false;
	Vector3D D = worldToModel * ray.dir;
	//D.normalize();
	Point3D O = worldToModel * ray.origin;
	double t_value, t_value2;

	// A = D.z*D.z-D.x*D.x-D.y*D.y		O = origin, D = dir
	// B = 2.0*(O.z*D.z - O.x*D.x - O.y*D.y)
	// C = O.z*O.z + k - O.x*O.x - O.y*O.y
	//double d_A = (D[2]*D[2]) - (D[0]*D[0]) - (D[1]*D[1]);
	double d_A = 2*D[0]*D[2] -(D[0]*D[0]) - 2*D[0]*D[1] ;
	//double d_B = 2.0 * ((O[2]*D[2]) - (O[0]*D[0]) - (O[1]*D[1]));
	double d_B = 2.0*(-(O[0]*D[0]) - (O[0]*D[1] + D[0]*O[1]) + (O[0]*D[2] + D[0]*O[2]));
	//double d_C = (O[2]*O[2]) + 1 - (O[0]*O[0]) - (O[1]*O[1]);
	double d_C = -(O[0]*O[0]) - 2.0*(O[0]*O[1]) + 2.0*O[0]*O[2] + 1;
	// Calculate d_D = B^2 - 4AC; 
	double d_D = (d_B * d_B) - (4 * d_A * d_C);

	if (d_D == 0.0 && d_A != 0.0)
	{
		// 1 intersect. 
		b_isHit = true;
		t_value = - d_B / (2 * d_A);
	}
	else if (d_D > 0.0 && d_A != 0.0)
	{
		// 2 intersects, have to pick correct one, t_value > t_value2
		t_value = - (d_B / (2*d_A) ) + (sqrt(d_D) / (2*d_A));
		t_value2 = - (d_B / (2*d_A) ) - (sqrt(d_D) / (2*d_A));
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
	else if (d_A == 0.0)
	{
		t_value = -d_C / d_B;
		b_isHit = true;
	}

	// Populate ray.intersection 
	if (b_isHit)
	{
		if (ray.intersection.none || ray.intersection.t_value > t_value)
		{
			Point3D intersectPt = O + (t_value * D);
			ray.intersection.none = false;
			ray.intersection.t_value = t_value;
			if (!b_shadowRay)
			{
				ray.intersection.point = modelToWorld * intersectPt;
				// f = P.z^2 + k - P.x^2 - P.y^2 = 0
				// gradient(f) = ( -2P.x, -2P.y, 1 + 2Pz)
				Point3D outPt = Point3D(-2 * intersectPt[0],
					-2 * intersectPt[1], 1 + 2 * intersectPt[2]);
				Vector3D normal = transNorm( worldToModel, (Point3D(0.0, 0.0, 0.0) - outPt) );
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