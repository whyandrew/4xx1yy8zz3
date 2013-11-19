/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		light source classes

***********************************************************/

#include "util.h"

// Base class for a light source.  You could define different types
// of lights here, but point light is sufficient for most scenes you
// might want to render.  Different light sources shade the ray 
// differently.
class LightSource {
public:
	virtual void shade( Ray3D& ) = 0;
	virtual Point3D get_position() const = 0; 
};

// A point light is defined by its position in world space and its
// colour.
class PointLight : public LightSource {
public:
	PointLight( Point3D pos, Colour col ) : _pos(pos), _col_ambient(col), 
	_col_diffuse(col), _col_specular(col) {}
	PointLight( Point3D pos, Colour ambient, Colour diffuse, Colour specular ) 
	: _pos(pos), _col_ambient(ambient), _col_diffuse(diffuse), 
	_col_specular(specular) {}
	void shade( Ray3D& ray );
	Point3D get_position() const { return _pos; }
	
private:
	Point3D _pos;
	Colour _col_ambient;
	Colour _col_diffuse; 
	Colour _col_specular; 
};

// Rendering mode, one bit per mode
enum mode {
	SIGNATURE = 0,	// Scene signature
	NO_SPECULAR = 1,	// Ambient + diffuse lighting using Phong model
	FULL_PHONG = 1 << 1,			// Ambient + diffuse + specular using Phong Model
	SPECULAR_ONLY = 1 << 2,
	AMBIENT_ONLY = 1 << 3,
	DIFFUSE_ONLY = 1 << 4

};
	
// Render mode
extern mode _render_mode;