/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		Implementations of functions in raytracer.h, 
		and the main function which specifies the 
		scene to be rendered.	

***********************************************************/


#include "raytracer.h"
#include "bmp_io.h"
#include <cmath>
#include <iostream>
#include <thread>
#include <ctime>
#include "material.h"

mode _render_mode = MODE_SIGNATURE;
int _reflect_rays;
int _reflect_depth;
double _reflect_fudge_factor;
double _refract_global_factor;

Raytracer::Raytracer() : _lightSource(NULL) {
	_root = new SceneDagNode();
	if (_render_mode & (mode)(MODE_REFLECT | MODE_REFRACT))
	{
		_reflect_depth = 6;
		_reflect_fudge_factor = 4;
		_refract_global_factor = 2;
	}
}

Raytracer::~Raytracer() {
	delete _root;
}

SceneDagNode* Raytracer::addObject( SceneDagNode* parent, 
		SceneObject* obj, Material* mat ) {
	SceneDagNode* node = new SceneDagNode( obj, mat );
	node->parent = parent;
	node->next = NULL;
	node->child = NULL;
	
	// Add the object to the parent's child list, this means
	// whatever transformation applied to the parent will also
	// be applied to the child.
	if (parent->child == NULL) {
		parent->child = node;
	}
	else {
		parent = parent->child;
		while (parent->next != NULL) {
			parent = parent->next;
		}
		parent->next = node;
	}
	
	return node;;
}

LightListNode* Raytracer::addLightSource( LightSource* light ) {
	LightListNode* tmp = _lightSource;
	_lightSource = new LightListNode( light, tmp );
	return _lightSource;
}

void Raytracer::rotate( SceneDagNode* node, char axis, double angle ) {
	Matrix4x4 rotation;
	double toRadian = 2*M_PI/360.0;
	int i;
	
	for (i = 0; i < 2; i++) {
		switch(axis) {
			case 'x':
				rotation[0][0] = 1;
				rotation[1][1] = cos(angle*toRadian);
				rotation[1][2] = -sin(angle*toRadian);
				rotation[2][1] = sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'y':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][2] = sin(angle*toRadian);
				rotation[1][1] = 1;
				rotation[2][0] = -sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'z':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][1] = -sin(angle*toRadian);
				rotation[1][0] = sin(angle*toRadian);
				rotation[1][1] = cos(angle*toRadian);
				rotation[2][2] = 1;
				rotation[3][3] = 1;
			break;
		}
		if (i == 0) {
		    node->trans = node->trans*rotation; 	
			angle = -angle;
		} 
		else {
			node->invtrans = rotation*node->invtrans; 
		}	
	}
}

void Raytracer::translate( SceneDagNode* node, Vector3D trans ) {
	Matrix4x4 translation;
	
	translation[0][3] = trans[0];
	translation[1][3] = trans[1];
	translation[2][3] = trans[2];
	node->trans = node->trans*translation; 	
	translation[0][3] = -trans[0];
	translation[1][3] = -trans[1];
	translation[2][3] = -trans[2];
	node->invtrans = translation*node->invtrans; 
}

void Raytracer::scale( SceneDagNode* node, Point3D origin, double factor[3] ) {
	Matrix4x4 scale;
	
	scale[0][0] = factor[0];
	scale[0][3] = origin[0] - factor[0] * origin[0];
	scale[1][1] = factor[1];
	scale[1][3] = origin[1] - factor[1] * origin[1];
	scale[2][2] = factor[2];
	scale[2][3] = origin[2] - factor[2] * origin[2];
	node->trans = node->trans*scale; 	
	scale[0][0] = 1/factor[0];
	scale[0][3] = origin[0] - 1/factor[0] * origin[0];
	scale[1][1] = 1/factor[1];
	scale[1][3] = origin[1] - 1/factor[1] * origin[1];
	scale[2][2] = 1/factor[2];
	scale[2][3] = origin[2] - 1/factor[2] * origin[2];
	node->invtrans = scale*node->invtrans; 
}

Matrix4x4 Raytracer::initInvViewMatrix( Point3D eye, Vector3D view, 
		Vector3D up ) {
	Matrix4x4 mat; 
	Vector3D w;
	view.normalize();
	up = up - up.dot(view)*view;
	up.normalize();
	w = view.cross(up);

	mat[0][0] = w[0];
	mat[1][0] = w[1];
	mat[2][0] = w[2];
	mat[0][1] = up[0];
	mat[1][1] = up[1];
	mat[2][1] = up[2];
	mat[0][2] = -view[0];
	mat[1][2] = -view[1];
	mat[2][2] = -view[2];
	mat[0][3] = eye[0];
	mat[1][3] = eye[1];
	mat[2][3] = eye[2];

	return mat; 
}

void Raytracer::traverseScene( SceneDagNode* node, Ray3D& ray, 
							  Matrix4x4* modelToWorld, Matrix4x4* worldToModel,
							  bool b_shadowRay) {
	SceneDagNode *childPtr;

	if (b_shadowRay)
	{
		if (!ray.intersection.none)
		{	// For shadow rays, can stop once intersection is detected

					//TODO : Need to add cases if the intersection is a refractive object

			return;
		}
		else if (ray.intersection.t_value > 0.999)
		{
			// if intersection is beyond light, there is also no shadow
			// Assume shadowRay dir is NOT normalize, so light is at t=0.999
			ray.intersection.none = true;
			return;
		}
	}
	// Applies transformation of the current node to the global
	// transformation matrices.
	*modelToWorld = *modelToWorld*node->trans;
	*worldToModel = node->invtrans * *worldToModel; 

	if (node->obj) 
	{
		// Perform intersection.
		if (node->obj->intersect(ray, *worldToModel, *modelToWorld, b_shadowRay)) 
		{
			ray.intersection.mat = node->mat;
		}
	}
	// Traverse the children.
	childPtr = node->child;
	while (childPtr != NULL) 
	{
		traverseScene(childPtr, ray, modelToWorld, worldToModel, b_shadowRay);
		childPtr = childPtr->next;
	}

	// Removes transformation of the current node from the global
	// transformation matrices.
	*worldToModel = node->trans * *worldToModel;
	*modelToWorld = *modelToWorld * node->invtrans;
}

// Return a randomized point near the given center
// on a plane perpendicular to the incident direction
// radius is the max linear distance allowed from center point
Point3D getRandomPt(Vector3D inDir, Point3D center, double radius)
{
	Point3D retPt;
	// Randomized radius and angle
	double ranAngle, ranRadius;
	ranRadius = (double)rand() / RAND_MAX;
	ranRadius = ranRadius * radius;
	ranAngle = (double)rand() / RAND_MAX;
	ranAngle = ranAngle * 2.0 * M_PI;

	double u_coor = ( ranRadius * cos(ranAngle) );
	double v_coor = ( ranRadius * sin(ranAngle) );

	// Get orthogonal vectors pair to incident vector
	Vector3D vect01;
	if (inDir[1] != 0 || inDir[2] != 0)
	{
		vect01 = Vector3D(1.0, 0.0, 0.0);
	}
	else
	{
		vect01 = Vector3D(0.0, 1.0, 0.0);
	}
	
	Vector3D u_vect = vect01.cross(inDir);
	Vector3D v_vect = u_vect.cross(inDir);
	u_vect.normalize();
	v_vect.normalize();

	retPt = center + (u_coor * u_vect) + (v_coor * v_vect);

	return retPt;
}

void Raytracer::computeShading( Ray3D& ray, 
							   Matrix4x4* modelToWorld, Matrix4x4* worldToModel ) {
	LightListNode* curLight = _lightSource;
	for (;;) {
		if (curLight == NULL) break;
		// Each lightSource provides its own shading function.

		// Implement shadows.
		// ray should contain information about intersection now
		if (_render_mode &
			(MODE_SHADOW | MODE_SOFTSHADOW_HIGH | MODE_SOFTSHADOW_LOW | MODE_SOFTSHADOW_EXTREME))
		{
			// Vector from intersect pt to light source
			// DO NOT normalize this vect, so light pos is at ~ t=0.999
			Vector3D vec_toLight = (curLight->light->get_position() - ray.intersection.point);
			// Advance the intersection point by a small delta
			Point3D deltaPt = ray.intersection.point + (0.001 * vec_toLight);
			double percentLight = 0;

			if (_render_mode & 
				(MODE_SOFTSHADOW_HIGH | MODE_SOFTSHADOW_LOW | MODE_SOFTSHADOW_EXTREME) )
			{
				// Implement area light, which always faces the incident ray
				// General randomized light positions based on light->pos()
				// 10 or 50 or 200 rays for shawdow checking
				int numRays;
				if (_render_mode & MODE_SOFTSHADOW_EXTREME)
					numRays = 500;
				else if (_render_mode & MODE_SOFTSHADOW_HIGH)
					numRays = 50;
				else 
					numRays = 10;

				double radius = 0.3;

				// If the FIRST consecutive shadowRays not in shadow
				// assume the point is not in shadow to save computation time
				int countNotShadow = 0;
				const int maxNotShadow = 15;
				int countInShadow = 0;
				const int maxInShadow = 40;

				for (int i = 0; i < numRays; i++, 
					countNotShadow < maxNotShadow, countInShadow < maxInShadow)
				{
					Point3D ranLightPOS 
						= getRandomPt(vec_toLight, curLight->light->get_position(), radius);
					Vector3D shadowRayDir = (ranLightPOS - deltaPt);
					Ray3D shadowRay(deltaPt, shadowRayDir);
					SceneDagNode *local_root = _root;
					traverseScene(local_root, shadowRay, modelToWorld, worldToModel, true);
					if (shadowRay.intersection.none)
					{
						percentLight += 1.0;
						if (i < maxNotShadow) {
							countNotShadow++; 
						}
					}
					else
					{
						if (i < maxInShadow) {
							countInShadow++;
						}
					}
				}
				
				// Calculate average light percentage
				if (countNotShadow >= maxNotShadow)
				{
					percentLight = 1.0;
				}
				else if (countInShadow >= maxInShadow)
				{
					percentLight = 0.0;
				}
				else
				{
					percentLight = percentLight / numRays;
				}
			}
			else // (_render_mode & MODE_SHADOW)
			{
				// Create a new ray to light source
				Ray3D shadowRay(deltaPt, vec_toLight);
				// Check for any intersection in scene
				SceneDagNode *local_root = _root;
				traverseScene(local_root, shadowRay, modelToWorld, worldToModel, true);
				// Shade original ray, depending if in shadow
				// If shadowRay hits something, then in shadow.
				percentLight = (shadowRay.intersection.none)? 1.0: 0.0;
			}

			curLight->light->shade(ray, percentLight);
		}
		else
		{
			// not in MODE_SHADOW
			curLight->light->shade(ray, 1.0);
		}

		curLight = curLight->next;
	}
	ray.col.clamp();
}

void Raytracer::initPixelBuffer() {
	int numbytes = _scrWidth * _scrHeight * sizeof(unsigned char);
	_rbuffer = new unsigned char[numbytes];
	_gbuffer = new unsigned char[numbytes];
	_bbuffer = new unsigned char[numbytes];
	for (int i = 0; i < _scrHeight; i++) {
		for (int j = 0; j < _scrWidth; j++) {
			_rbuffer[i*_scrWidth+j] = 0;
			_gbuffer[i*_scrWidth+j] = 0;
			_bbuffer[i*_scrWidth+j] = 0;
		}
	}
}

void Raytracer::flushPixelBuffer( char *file_name ) {
	bmp_write( file_name, _scrWidth, _scrHeight, _rbuffer, _gbuffer, _bbuffer );
	delete _rbuffer;
	delete _gbuffer;
	delete _bbuffer;
}

Colour Raytracer::shadeRay( Ray3D& ray, int depth,
						   Matrix4x4* modelToWorld, Matrix4x4* worldToModel) {
	Colour col(0.0, 0.0, 0.0); 
	SceneDagNode *local_root = _root;
	bool b_internalReflection = false;
	// All input ray should have ray.refractIndex set already
	traverseScene(local_root, ray, modelToWorld, worldToModel, false); 
	
	// Don't bother shading if the ray didn't hit 
	// anything.
	if (!ray.intersection.none && depth >= 0) 
	{
		computeShading(ray, modelToWorld, worldToModel); 

		// Add effects from reflection and refraction
		// Follow the general scheme on p306-7
		if (depth > 0 && _render_mode & MODE_REFRACT)
		{
			double reflectance;
			double rindex = ray.intersection.mat->refract_index;
			double dotProd;
			Ray3D refractRay;
			Ray3D reflectRay;
			double cosTheta;
			// A bit of a hack, to allow adjusting global refract & reflect effect
			double refractBeerCoeff;
			double reflectBeerCoeff;

			if ( rindex > 1.0)
			{
				// ray intersect transparent object
				reflectRay = getReflectRay(ray);
				
				dotProd = ray.intersection.normal.dot(ray.dir);
				if (dotProd < 0.0)
				{
					// Ray hits OUTSIDE of object
					// do refract, get refract.dir
					getRefractRay(ray, &refractRay, true);
					cosTheta = -dotProd;
					// TODO: why = 1.0?? Assume no absorbance in air?
					refractBeerCoeff = 1.0;
					reflectBeerCoeff = 1.0;
				}
				else
				{
					// Ray hits INSIDE of object
					//refractBeerCoeff = exp(-ray.intersection.t_value) / _refract_global_factor;
					//reflectBeerCoeff = exp(-ray.intersection.t_value) / _reflect_fudge_factor;
					refractBeerCoeff = 1.0;
					reflectBeerCoeff = 1.0;

					if (getRefractRay(ray, &refractRay, false))
					{
						// Have refraction
						cosTheta = refractRay.dir.dot(ray.intersection.normal);
					}
					else
					{
						// Total interal reflection, no refraction
						ray.col = ray.col + reflectBeerCoeff * 
							shadeRay(reflectRay, depth - 1, modelToWorld, worldToModel);
					}
				}

				double r0 = ray.intersection.mat->reflectance;
				reflectance = r0 + ((1 - r0) * pow( (1 - cosTheta), 5 ));

				ray.col = ray.col + 
					refractBeerCoeff * (1 - reflectance) *
					shadeRay(refractRay, depth - 1, modelToWorld, worldToModel);

				if (_render_mode & MODE_REFLECT)
				{
					ray.col = ray.col +
						reflectBeerCoeff * reflectance *
						shadeRay(reflectRay, depth - 1, modelToWorld, worldToModel);
				}
			}
			else
			{
				// intersect opaque object, reflect only
				ray.col = ray.col + 
					getReflectRayColor(ray, depth - 1, modelToWorld, worldToModel);
			}

		}
		else if (_render_mode & MODE_REFLECT)
		{
			ray.col = ray.col + 
					getReflectRayColor(ray, depth - 1, modelToWorld, worldToModel);
		}

		ray.col.clamp();
		col = ray.col;
	}
	return col; 
}	

Ray3D Raytracer::getReflectRay( Ray3D& ray)
{
	Vector3D rayDir = ray.dir;
	rayDir.normalize();
	Vector3D reflectRay_dir = rayDir - 2 * rayDir.dot(ray.intersection.normal) * ray.intersection.normal;
	reflectRay_dir.normalize();
	Point3D intersectPtDelta = ray.intersection.point + ( 0.001 * reflectRay_dir);
	Ray3D reflectRay = Ray3D(intersectPtDelta, reflectRay_dir);

	return reflectRay;
}


bool Raytracer::getRefractRay( Ray3D& ray, Ray3D *refractRay, bool b_hitOutside)
{
	bool b_haveRefract = true;

	Vector3D rDir = ray.dir;
	rDir.normalize();

	Vector3D rNormal = b_hitOutside? ray.intersection.normal: -ray.intersection.normal;
	rNormal.normalize();
	
	double inRI = ray.refract_index;
	double outRI = ray.intersection.mat->refract_index;
	//ratio of in/out refractive index
	double ratioRI;
	ratioRI = b_hitOutside?  (inRI / outRI): (outRI / inRI);

	// TODO: have to switch in/out RI depending on b_hitOutside
	//			since theta is always angle from air (or just less RI side?).

	double cosTheta = rDir.dot(rNormal);

	double cosPhiSq = 1 - ( (1 - cosTheta * cosTheta) * pow(ratioRI, 2) );
	Vector3D refractDir;

	if (cosPhiSq < 0.0)
	{
		// Total internal reflection
		b_haveRefract = false;
	}
	else
	{
		refractDir = ( ratioRI  * (rDir - (cosTheta * rNormal))) 
			- (sqrt(cosPhiSq) * rNormal);
		refractDir.normalize();
	}

	if (b_haveRefract)
	{
		Point3D intersectDelta = ray.intersection.point + 0.001 * refractDir;
		*refractRay = Ray3D(intersectDelta, refractDir);
	}

	return b_haveRefract;
}

Colour Raytracer::getReflectRayColor( 
	Ray3D& ray, 
	int depth,
	Matrix4x4* modelToWorld, 
	Matrix4x4* worldToModel)
{
	Colour totalReflectColor(0.0, 0.0, 0.0);

	if (depth > 0)
	{
		//depth--; // Handle depth decrement in "shadeRay"

		Ray3D reflectRay = getReflectRay(ray);

		Colour reflectColor = shadeRay(reflectRay, depth, modelToWorld, worldToModel);
		if (!reflectRay.intersection.none)
		{
			//Dampening + reflectivity affects factor
					
			if (ray.intersection.mat->reflect_factor >= 1.0)
			{
				// Assume reflect_factor of >= 1 is glass
				// Minus the diffuse/ambient component of glass,
				// so image doesn't get brighter than object
				Colour mirror_image = ((reflectColor - 
					ray.intersection.mat->diffuse) - ray.intersection.mat->ambient);
				mirror_image.clamp();
				totalReflectColor = totalReflectColor + mirror_image;
			}
			else 
			{
				// Beer's Law: I = Io * e^(-a*x)
				// a = 1/reflect_factor, x = t_value / _reflect_fudge_factor
				double factor = exp( 
					-( (ray.intersection.t_value / _reflect_fudge_factor) 
					/ ray.intersection.mat->reflect_factor));
				totalReflectColor = totalReflectColor + (factor * reflectColor);
			}
		}
		else
		{
			// Stop recursion if reflectRay not hitting anything
			depth = 0;
		}

		totalReflectColor.clamp();
	}

	return totalReflectColor;
}

void Raytracer::render( int width, int height, Point3D eye, Vector3D view, 
		Vector3D up, double fov, char* fileName ) {
	Matrix4x4 viewToWorld;
	_scrWidth = width;
	_scrHeight = height;
	double factor = (double(height)/2) / tan(fov*M_PI/360.0);

	initPixelBuffer();
	viewToWorld = initInvViewMatrix(eye, view, up);
	
	if (_render_mode & MODE_MULTITHREAD)
	{
		// Seting up 8 threads... hardcoding not pretty!
		int segment = _scrWidth / 8;
		// Arguments for threads
		renderArgs args1(0, segment, factor, viewToWorld);
		renderArgs args2(segment, 2*segment, factor, viewToWorld);
		renderArgs args3(2*segment, 3*segment, factor, viewToWorld);
		renderArgs args4(3*segment, 4*segment, factor, viewToWorld);
		renderArgs args5(4*segment, 5*segment, factor, viewToWorld);
		renderArgs args6(5*segment, 6*segment, factor, viewToWorld);
		renderArgs args7(6*segment, 7*segment, factor, viewToWorld);
		renderArgs args8(7*segment, _scrWidth, factor, viewToWorld);
		// Starts threads
		std::thread t1 (&Raytracer::render_section, this, &args1);
		std::thread t2 (&Raytracer::render_section, this, &args2);
		std::thread t3 (&Raytracer::render_section, this, &args3);
		std::thread t4 (&Raytracer::render_section, this, &args4);
		std::thread t5 (&Raytracer::render_section, this, &args5);
		std::thread t6 (&Raytracer::render_section, this, &args6);
		std::thread t7 (&Raytracer::render_section, this, &args7);
		std::thread t8 (&Raytracer::render_section, this, &args8);
		// Wait on threads
		t1.join();
		t2.join();
		t3.join();
		t4.join();
		t5.join();
		t6.join();
		t7.join();
		t8.join();
	}
	else
	{
		renderArgs args(0, _scrWidth, factor, viewToWorld);
		render_section(&args);
	}
	
	flushPixelBuffer(fileName);
}

void Raytracer::render_section(renderArgs* args) 
{
	int i_start = args->width_start;
	int i_end = args->width_end>_scrWidth? _scrWidth: args->width_end;
	Matrix4x4 viewToWorld = args->viewToWorld;
	double factor = args->factor;
	int width = _scrWidth;
	int height = _scrHeight;
	Matrix4x4 modelToWorld;
	Matrix4x4 worldToModel; 
	

	// Construct a ray for each pixel.
	for (int i = 0; i < _scrHeight; i++) {
		for (int j = i_start; j < i_end; j++) {

			// Super-sample antialias modes
			int i_ssaa = 1;
			if (_render_mode & MODE_SSAA4)
			{
				i_ssaa = 2;
			}
			else if (_render_mode & MODE_SSAA16)
			{
				i_ssaa = 4;
			}
			else if (_render_mode & MODE_SSAA36)
			{
				i_ssaa = 6;
			}
			else if (_render_mode & MODE_SSAA64)
			{
				i_ssaa = 8;
			}

			Colour col(0.0, 0.0, 0.0);
			// Sets up ray origin 
			Point3D origin(0, 0, 0);
			Point3D imagePlane;

			// Stratified sampling described in text p311
			for (int p = 0; p < i_ssaa; p++)
			{
				for (int q = 0; q < i_ssaa; q++)
				{
					double d_rand;
					if (i_ssaa == 1)
					{
						d_rand = 0.5;
					}
					else
					{
						d_rand = (double)rand() / RAND_MAX;
					}
					
					// Sets up direction in view space, 
					// image plane is at z = -1.
					imagePlane[0] = (-double(width)/2 + (p + d_rand)/i_ssaa + j)/factor;
					imagePlane[1] = (-double(height)/2 + (q + d_rand)/i_ssaa + i)/factor;
					imagePlane[2] = -1;

					// Convert ray to world space and call 
					// shadeRay(ray) to generate pixel colour. 	
			
					Ray3D ray;
					ray.origin = viewToWorld * origin;
					ray.dir = viewToWorld * (imagePlane - origin);
					// Assume all ray starting from viewplane is in air
					ray.refract_index = RI_AIR;

					col = col + shadeRay(ray, _reflect_depth, &modelToWorld, &worldToModel ); 
				}
			}
			// Take average of all sample rays
			col = (1.0 / (double)(i_ssaa * i_ssaa)) * col;

			_rbuffer[i*width+j] = int(col[0]*255);
			_gbuffer[i*width+j] = int(col[1]*255);
			_bbuffer[i*width+j] = int(col[2]*255);
		}
	}
}

int main(int argc, char* argv[])
{	
	// Keep a timer
	std::clock_t start_time;
	double duration;
	start_time = std::clock();
	// seed random
	srand(time(0));
	
	/************************************************
			 Set _render_mode here
	************************************************/
	_render_mode = (mode)(MODE_FULL_PHONG | MODE_MULTITHREAD | MODE_SOFTSHADOW_LOW);

	Raytracer raytracer;

	int width = 300; 
	int height = 200; 

	if (argc == 3) {
		width = atoi(argv[1]);
		height = atoi(argv[2]);
	}

	/********************************************************************
		SCENE 1: 
			"Room" with 2 mirror spheres and 2 solid sphere and 2 lights.
	*********************************************************************/



	raytracer.addLightSource( new PointLight(Point3D(0, 3, 1), Colour(1.0, 1.0, 1.0)));
	//raytracer.addLightSource( new PointLight(Point3D(1, 0, 1), Colour(0.5,0.5,0.5)));
	//raytracer.addLightSource( new PointLight(Point3D(0, 3.9, 1), Colour(0.5,0.5,0.5)));

	//SceneDagNode* sphere1 = raytracer.addObject( new UnitSphere(), &mat_glass);
	SceneDagNode* sphere2 = raytracer.addObject( new Hyperboloid(), &mat_copper );
	SceneDagNode* sphere3 = raytracer.addObject( new UnitSphere(), &texture_neptune );
	SceneDagNode* sphere4 = raytracer.addObject( new UnitSphere(), &mat_yellow );

	SceneDagNode* plane_back = raytracer.addObject( new UnitSquare(), &mat_purple );
	SceneDagNode* plane_bottom = raytracer.addObject( new UnitSquare(), &mat_green);
	SceneDagNode* plane_left = raytracer.addObject( new UnitSquare(), &mat_red);
	//SceneDagNode* plane_top = raytracer.addObject( new UnitSquare(), &mat_gold);
	//SceneDagNode* plane_right = raytracer.addObject( new UnitSquare(), &mat_chrome);

	double factor3[3] = {0.5, 0.5, 0.5};
	//raytracer.translate(sphere1, Vector3D(-0.5, 0, -2.5));
	//raytracer.scale(sphere1, Point3D(0,0,0), factor3);

	double factor1[3] = {0.3, 0.3, 0.3};
	raytracer.translate(sphere2, Vector3D(2, 1, -1));
	raytracer.scale(sphere2, Point3D(0,0,0), factor1);

	raytracer.translate(sphere4, Vector3D(1, 1.5, -3));
	//raytracer.scale(sphere4, Point3D(0,0,0), factor1);

	
	raytracer.translate(sphere3, Vector3D(-0.75, -0.5, -3));
	raytracer.scale(sphere3, Point3D(0,0,0), factor3);

	double factor2[3] = { 8.0, 8.0, 8.0 };
    raytracer.translate(plane_back, Vector3D(0, 0, -7));        
    raytracer.scale(plane_back, Point3D(0, 0, 0), factor2);

	raytracer.translate(plane_bottom, Vector3D(0, -2, -3));        
    raytracer.scale(plane_bottom, Point3D(0, 0, 0), factor2);
	raytracer.rotate(plane_bottom, 'x', -90);

	raytracer.translate(plane_left, Vector3D(-4, 0, -3));        
    raytracer.scale(plane_left, Point3D(0, 0, 0), factor2);
	raytracer.rotate(plane_left, 'y', 90);
	/*
	raytracer.translate(plane_top, Vector3D(0, 4, -3));        
    raytracer.scale(plane_top, Point3D(0, 0, 0), factor2);
	raytracer.rotate(plane_top, 'x', 90);

	raytracer.translate(plane_right, Vector3D(4, 0, -3));        
    raytracer.scale(plane_right, Point3D(0, 0, 0), factor2);
	raytracer.rotate(plane_right, 'y', -90);
	*/
	// Camera parameters.
	Point3D eye(0, 0, 1);
	Vector3D view(0, 0, -1);
	Vector3D up(0, 1, 0);
	double fov = 40;

	// Render the scene, feel free to make the image smaller for
	// testing purposes.	
	//raytracer.render(width, height, eye, view, up, fov, "AA1.bmp");
	
	// Render it from a different point of view.
	Point3D eye2(2, 2, 1);
	Vector3D view2(-2, -2, -6);
	raytracer.render(width, height, eye2, view2, up, fov, "AA_36.bmp");

	// Calculate time taken
	duration = (std::clock() - start_time) / 1000.0;
	std::cout << "Finished rendering in " << duration << " sec." << '\n' << '\n';

	return 0;
}
