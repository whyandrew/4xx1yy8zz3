/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		light source classes

***********************************************************/

#ifndef _UTIL_
#include "util.h"
#endif

#include "scene_object.h"

#ifndef H_RAY
#include "raytracer.h"
#endif

#ifndef H_LIGHTSOURCE
#define H_LIGHTSOURCE

// Base class for a light source.  You could define different types
// of lights here, but point light is sufficient for most scenes you
// might want to render.  Different light sources shade the ray 
// differently.
class LightSource {
public:
	//bool checkShadow( Raytracer*, Ray3D&, Point3D,
	//	Matrix4x4*, Matrix4x4*) {return false;}
	virtual void shade( Raytracer *obj, Ray3D& ray,
		Matrix4x4* modelToWorld, Matrix4x4* worldToModel) = 0;
	virtual Point3D get_position() const = 0; 
	enum lightType {
		LIGHT_POINT,
		LIGHT_SPOT,

		LIGHT_INVALID
	};
	lightType type; 
};

// Linked list containing light sources in the scene.
struct LightListNode {
	LightListNode() : light(NULL), next(NULL) {}
	LightListNode( LightSource* light, LightListNode* next = NULL ) : 
		light(light), next(next) {}
	~LightListNode() { 
		if (!light) delete light; 
	}
	LightSource* light;
	LightListNode* next;
};

// A point light is defined by its position in world space and its
// colour.
class PointLight : public LightSource {
public:
	PointLight( Point3D pos, Colour col ) : _pos(pos), _col_ambient(col), 
		_col_diffuse(col), _col_specular(col) , type(LIGHT_POINT) {}

	PointLight( Point3D pos, Colour ambient, Colour diffuse, Colour specular ):
		_pos(pos), _col_ambient(ambient), _col_diffuse(diffuse), 
		_col_specular(specular) , type(LIGHT_POINT) {}

	void PointLight::shade( Raytracer *obj, Ray3D& ray, Matrix4x4* modelToWorld, Matrix4x4* worldToModel );
	Point3D get_position() const { return _pos; }
	lightType type;

private:
	Point3D _pos;
	Colour _col_ambient;
	Colour _col_diffuse; 
	Colour _col_specular; 
};

// A spot light is unit square, modeled by a number of separate point lights
// Its square is perpendicular to "axis", but otherwise non-directional
class SpotLight : public LightSource {
public:
	SpotLight( Point3D pos, char axis, float size, Colour col ) : 
		_centerpos(pos), _col_ambient(col), _size(size),
		_col_diffuse(col), _col_specular(col) , type(LIGHT_POINT), _axis(axis) 
	{
		setupLights();
	}

	SpotLight( Point3D pos, char axis, float size, Colour ambient, 
		Colour diffuse, Colour specular):
		_centerpos(pos), _col_ambient(ambient), _col_diffuse(diffuse), 
		_col_specular(specular) , type(LIGHT_SPOT), _axis(axis), _size(size) 
	{
		setupLights();
	}

	void setupLights();
	void shade( Raytracer *obj, Ray3D& ray, Matrix4x4* modelToWorld, Matrix4x4* worldToModel );
	Point3D get_position();
	lightType type;

private:
	// num of point lights per side of its square shape
	static const int _numLight = 4; 
	//length of each side of square
	float _size;
	Point3D _centerpos;
	PointLight *_lights[_numLight * _numLight + 1];
	char _axis;
	Colour _col_ambient;
	Colour _col_diffuse; 
	Colour _col_specular; 
};

// Rendering mode, one bit per mode
enum mode {
	MODE_SIGNATURE = 0,	// Scene signature
	MODE_SPECULAR = 1 << 2,
	MODE_AMBIENT= 1 << 3,
	MODE_DIFFUSE= 1 << 4,
	MODE_FULL_PHONG = MODE_SPECULAR | MODE_AMBIENT | MODE_DIFFUSE,
	MODE_MULTITHREAD = 1 << 5,
	MODE_SSAA4 = 1 << 6,
	MODE_SSAA16 = 1 << 7,
	MODE_SSAA36 = 1 << 8,
	MODE_SSAA64 = 1 << 9,
	MODE_REFLECT = 1 << 10,
	MODE_REFRACT = 1 << 11,
	MODE_SHADOW = 1 << 12,

};
	
// Render mode
extern mode _render_mode;
extern int _reflect_rays;
extern int _reflect_depth;
extern int _refract_rays;

#endif