/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		classes defining primitives in the scene

***********************************************************/

#ifndef _UTIL_
#include "util.h"
#endif

// All primitives should provide a intersection function.  
// To create more primitives, inherit from SceneObject.
// Namely, you can create, Sphere, Cylinder, etc... classes
// here.
class SceneObject {
public:
	// Returns true if an intersection occured, false otherwise.
	// b_shadowRay indicates if it's a shadow ray, save unneeded calculations
	virtual bool intersect( Ray3D&, const Matrix4x4&, 
		const Matrix4x4&, bool b_shadowRay) = 0;

	// Solve quadratic for t-value, given A,B,C coefficients
	// Return true if there is a solution
	// Store the smaller t solution to p_tValue
	static bool solveT(double A, double B, double C, double* p_tValue);

	// Return mapped textural pixel color.
	// Need to be called by each light shade function
	// Does not set ray.col
	virtual Colour textureMapping(const Ray3D& ray) = 0;
};

// Subclass for object composing of multiple objects
class CompoundObject : public SceneObject
{
public:
	// All compound object can use the same intersect function
	// Loop through intersect function of each object in construction
	bool intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld, bool b_shadowRay );
		//virtual void construct();
protected:
	CompoundObject(  ): numObj(0) {  };
	~CompoundObject() 
	{
		if (*p_objPList != NULL)
		{ 
			for (int i = 0; i < numObj; i++)
				delete p_objPList[i];
			delete [] *p_objPList;
		}
	}
	// Function to construct the compound object
	// Save object points to list
	//virtual void construct();
	// Pointer to a list of SceneObject
	SceneObject **p_objPList;
	// Total number of object used to construct compound object
	// Inheriting class must set this number
	int numObj;
};

// Example primitive you can create, this is a unit square on 
// the xy-plane.
class UnitSquare : public SceneObject {
public:
	bool intersect( Ray3D& ray, const Matrix4x4& worldToModel,
			const Matrix4x4& modelToWorld, bool b_shadowRay );
	Colour textureMapping(const Ray3D& ray);
};

class UnitSphere : public SceneObject {
public:
	bool intersect( Ray3D& ray, const Matrix4x4& worldToModel,
			const Matrix4x4& modelToWorld, bool b_shadowRay );
	Colour textureMapping(const Ray3D& ray);
};

// Hyperboloid (open-ended) lays flat along z-axis
class _Hyperboloid : public SceneObject {
public :
	_Hyperboloid(): _zRange(0.5) {}
	// zRange determines how long the hyperboloid is, default = 1
	_Hyperboloid(double height): _zRange(height/2.0) {}

	bool intersect( Ray3D& ray, const Matrix4x4& worldToModel,
			const Matrix4x4& modelToWorld, bool b_shadowRay );

	Colour textureMapping(const Ray3D& ray);
private :
	double _zRange;
};

// Circle is parellel to xy-plane, cut across z=0
// Normal is facing to +z by default.
class _Circle: public SceneObject {
public :
	_Circle(): _radius(1.0), _flipNormal(false), _zvalue(0.0) {}
	_Circle(double radius): _radius(radius), _zvalue(0.0), _flipNormal(false) {}
	_Circle(double radius, double zvalue, bool flipNormal): 
		_radius(radius), _zvalue(zvalue),_flipNormal(flipNormal) {}

	bool intersect( Ray3D& ray, const Matrix4x4& worldToModel,
			const Matrix4x4& modelToWorld, bool b_shadowRay );

	Colour textureMapping(const Ray3D& ray);
private :
	double _radius;
	double _zvalue;
	bool _flipNormal;
};

// Hyperboloid close-ended with 2 circlar planes
class Hyperboloid : public CompoundObject {
public :
	Hyperboloid(): _zRange(0.5) 
	{
		numObj = 3; construct(); 
	}
	// zRange determines how long the hyperboloid is, default = 1
	Hyperboloid(double height): _zRange(height/2.0)
	{ 
		numObj = 3;
		construct(); 
	}

	//bool intersect( Ray3D& ray, const Matrix4x4& worldToModel,
			//const Matrix4x4& modelToWorld, bool b_shadowRay );
private :
	double _zRange;
public:
	// Setup a open hyperboloid and 2 circle planes
	void construct();

	Colour textureMapping(const Ray3D& ray);
};