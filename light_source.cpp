/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

#include <cmath>
#include "light_source.h"

bool checkShadow( Raytracer *obj, Ray3D& ray, Point3D lightPOS,
							 Matrix4x4* modelToWorld, Matrix4x4* worldToModel)
{
	// Vector from intersect pt to light source
	// DO NOT normalize this vect, so light pos is at ~ t=0.99
	Vector3D vec_toLight = (lightPOS - ray.intersection.point);
	// Advance the intersection point by a small delta
	Point3D intersectPlusPt = ray.intersection.point + (0.01 * vec_toLight);
	// Create a new ray to light source
	Ray3D shadowRay(intersectPlusPt, vec_toLight);
	// Check for any intersection in scene
	SceneDagNode *local_root = obj->getSceneRoot();
	obj->traverseScene(local_root, shadowRay, modelToWorld, worldToModel, true);

	return !shadowRay.intersection.none;
}

void PointLight::shade( Raytracer *obj, Ray3D& ray, Matrix4x4* modelToWorld, Matrix4x4* worldToModel ) {
	// TODO: implement this function to fill in values for ray.col 
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray 
	// is available.  So be sure that traverseScene() is called on the ray 
	// before this function.  

	// Scene Signature only
	if (_render_mode == MODE_SIGNATURE)
	{
		// Just use plain diffuse color, disregard lights
		ray.col = ray.intersection.mat->diffuse;
		//ray.col = ray.intersection.mat->specular;
		return;
	}

	// Check if intersect pt is under shadow
	bool b_inShadow = false;
	if (_render_mode & MODE_SHADOW)
	{
		b_inShadow = checkShadow( obj, ray, _pos, modelToWorld, worldToModel);
	}

	Colour newColor = Colour(0.0, 0.0, 0.0);
	Material* objMat = ray.intersection.mat;
	Vector3D objNormal = ray.intersection.normal;
	Point3D objPoint = ray.intersection.point;

	// Ambient
	if (_render_mode & MODE_AMBIENT)
	{
		// Add ambient, uniform non-directional
		newColor = newColor + (_col_ambient * objMat->ambient);
	}

	Vector3D vec_light = _pos - objPoint;
	vec_light.normalize();

	// Diffuse
	if (_render_mode & MODE_DIFFUSE && !b_inShadow)
	{
		// Add diffuse, mat*light*max(0,factor), factor = normal(dot)light
		double factor = objNormal.dot( vec_light );
		factor = factor>0.0? factor: 0.0;
		newColor = newColor + ( factor * (_col_diffuse * objMat->diffuse));
	}

	// Specular
	if (_render_mode & MODE_SPECULAR && !b_inShadow)
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
	ray.col = ray.col + newColor;
}

void SpotLight::setupLights()
{
	//Setup an array of pointlights to simulate spotlight
	_lights[0] = &PointLight(_centerpos, _col_ambient, _col_diffuse, _col_specular);
	int light_index = 1;
	float section = _size / (float)_numLight;
	float half = _size / 2.0;

	for (int i = 0; i < _numLight; i++)
	{
		for (int j = 0; j< _numLight; j++)
		{
			Point3D pos;
			if (_axis = 'x')
			{
				pos = Point3D(0, section * (i + 0.5) - half, section * (j + 0.5) - half);
			}
			else if (_axis = 'y')
			{
				pos = Point3D(section * (i + 0.5) - half, 0, section * (j + 0.5) - half);
			}
			else
			{
				pos = Point3D(section * (i + 0.5) - half, section * (j + 0.5) - half, 0);
			}

			_lights[light_index] = &PointLight(pos, _col_ambient, _col_diffuse, _col_specular); 
			light_index++;
		}
	}
}

void SpotLight::shade( Raytracer *obj, Ray3D& ray, Matrix4x4* modelToWorld, Matrix4x4* worldToModel)
{

}

Point3D SpotLight::get_position()
{
	return Point3D();
}
