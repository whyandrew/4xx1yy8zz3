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
	virtual void shade( Ray3D&, bool b_inShadow) = 0;
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
	void shade( Ray3D& ray, bool b_inShadow );
	Point3D get_position() const { return _pos; }
	
private:
	Point3D _pos;
	Colour _col_ambient;
	Colour _col_diffuse; 
	Colour _col_specular; 
};

// Rendering mode, one bit per mode
enum mode {
	MODE_SIGNATURE = 1,	// Scene signature
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
