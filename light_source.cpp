/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

#include <cmath>
#include "light_source.h"

void PointLight::shade( Ray3D& ray ) {
	// TODO: implement this function to fill in values for ray.col 
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray 
	// is available.  So be sure that traverseScene() is called on the ray 
	// before this function.  

	// Scene Signature only
	if (_render_mode == SIGNATURE)
	{
		// Just use plain diffuse color, disregard lights
		ray.col = ray.intersection.mat->diffuse;
		//ray.col = ray.intersection.mat->specular;
		return;
	}

	Colour newColor = Colour(0.0, 0.0, 0.0);
	Material* objMat = ray.intersection.mat;
	Vector3D objNormal = ray.intersection.normal;
	Point3D objPoint = ray.intersection.point;

	// Ambient
	if (_render_mode & (AMBIENT_ONLY | FULL_PHONG | NO_SPECULAR))
	{
		// Add ambient, uniform non-directional
		newColor = newColor + (_col_ambient * objMat->ambient);
	}

	Vector3D vec_light = _pos - objPoint;
	vec_light.normalize();

	// Diffuse
	if (_render_mode & (NO_SPECULAR | FULL_PHONG | DIFFUSE_ONLY))
	{
		// Add diffuse, mat*light*max(0,factor), factor = normal(dot)light
		double factor = objNormal.dot( vec_light );
		factor = factor>0.0? factor: 0.0;
		newColor = newColor + ( factor * (_col_diffuse * objMat->diffuse));
	}

	// Specular
	if (_render_mode & (FULL_PHONG | SPECULAR_ONLY))
	{
		// Add specular, mat*light*max(0,factor), factor = -ray(dot)reflect
		// reflect = 2(light(dot)normal)*normal - light
		Vector3D vec_reflect = (2 * objNormal.dot( vec_light ) * objNormal) - vec_light;
		vec_reflect.normalize();
		Vector3D vect_backRay = -ray.dir;
		vect_backRay.normalize();
		double factor = vect_backRay.dot(vec_reflect);
		factor = factor>0.0? factor: 0.0;
		newColor = newColor + ( pow(factor, objMat->specular_exp) * _col_specular * objMat->specular);
		//newColor = newColor + ( factor * _col_specular * objMat->specular);
	}

	newColor.clamp();

	ray.col = newColor;

}

