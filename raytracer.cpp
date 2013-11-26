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
int _refract_rays;

Raytracer::Raytracer() : _lightSource(NULL) {
	_root = new SceneDagNode();
	if (_render_mode & (mode)(MODE_REFLECT | MODE_REFRACT))
	{
		_reflect_depth = 5;
		_reflect_rays = 1;
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
			return;
		}
		else if (ray.intersection.t_value > 0.99)
		{
			// if intersection is beyond light, there is also no shadow
			// Assume shadowRay dir is NOT normalize, so light is at t=0.99
			ray.intersection.none = true;
			return;
		}
	}
	// Applies transformation of the current node to the global
	// transformation matrices.
	*modelToWorld = *modelToWorld*node->trans;
	*worldToModel = node->invtrans * *worldToModel; 
	if (node->obj) {
		// Perform intersection.
		if (node->obj->intersect(ray, *worldToModel, *modelToWorld, b_shadowRay)) {
			ray.intersection.mat = node->mat;
		}
	}
	// Traverse the children.
	childPtr = node->child;
	while (childPtr != NULL) {
		traverseScene(childPtr, ray, modelToWorld, worldToModel, b_shadowRay);
		childPtr = childPtr->next;
	}

	// Removes transformation of the current node from the global
	// transformation matrices.
	*worldToModel = node->trans * *worldToModel;
	*modelToWorld = *modelToWorld * node->invtrans;
}

void Raytracer::computeShading( Ray3D& ray, 
							   Matrix4x4* modelToWorld, Matrix4x4* worldToModel ) {
	LightListNode* curLight = _lightSource;
	for (;;) {
		if (curLight == NULL) break;
		// Each lightSource provides its own shading function.

		// Implement shadows.
		// ray should contain information about intersection now
		if (_render_mode & MODE_SHADOW)
		{
			// Vector from intersect pt to light source
			// DO NOT normalize this vect, so light pos is at ~ t=0.99
			Vector3D vec_toLight = (curLight->light->get_position() - ray.intersection.point);
			// Advance the intersection point by a small delta
			Point3D intersectPlusPt = ray.intersection.point + (0.01 * vec_toLight);
			// Create a new ray to light source
			Ray3D shadowRay(intersectPlusPt, vec_toLight);
			// Check for any intersection in scene
			SceneDagNode *local_root = _root;
			traverseScene(local_root, shadowRay, modelToWorld, worldToModel, true);
		
			// Shade original ray, depending if in shadow
			// If shadowRay hits something, then in shadow.
			curLight->light->shade(ray, !shadowRay.intersection.none);

		}
		else
		{
			// not in MODE_SHADOW
			curLight->light->shade(ray, false);
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

	traverseScene(local_root, ray, modelToWorld, worldToModel, false); 
	
	// Don't bother shading if the ray didn't hit 
	// anything.
	if (!ray.intersection.none) {
		computeShading(ray, modelToWorld, worldToModel); 
 
	// You'll want to call shadeRay recursively (with a different ray, 
	// of course) here to implement reflection/refraction effects.  

	// Generate secondary ray(s) from ray.intersection.point
	// either follow normal, or more rays
	// shaderay(secondary_ray);
	// add secondary_ray.col to ray.col, with some deminishing factor?? try (1/t-value)
		Colour totalReflectColor(0.0, 0.0, 0.0);
		bool b_done = false;
		if (depth > 0 && !b_done)
		{
			depth--;
			//for (int i = 0; i < _reflect_rays; i++)
			{
				Point3D intersectPtDelta = ray.origin + ((ray.intersection.t_value * 0.999) * ray.dir);
				
				Vector3D reflectRay_dir = ray.dir - 2 * ray.dir.dot(ray.intersection.normal) * ray.intersection.normal;

				reflectRay_dir.normalize();
				Ray3D reflectRay = Ray3D(intersectPtDelta, reflectRay_dir);
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
						// damping factor = 1 / (1 + distance)
						//double distance = reflectRay.intersection.t_value;
						//double factor = distance>=1.0? 1.0/distance: 2.0-distance;
						double factor = 1 / ( 1 + reflectRay.intersection.t_value);
						factor = ray.intersection.mat->reflect_factor * factor;
						totalReflectColor = totalReflectColor + (factor * reflectColor);
					}
				}
				else
				{
					// Stop recursion if reflectway not hitting anything
					b_done = true;
				}
			}
			totalReflectColor.clamp();
		}

		ray.col = ray.col + totalReflectColor;
		ray.col.clamp();
		col = ray.col;
	}
	return col; 
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

	//_render_mode = (mode)(MODE_SIGNATURE | MODE_MULTITHREAD);
	//_render_mode = (mode)(MODE_FULL_PHONG | MODE_MULTITHREAD);// | MODE_SSAA4);
	_render_mode = (mode)(MODE_FULL_PHONG | MODE_MULTITHREAD | MODE_SHADOW | MODE_REFLECT | MODE_SSAA4);
	//_render_mode = (mode) (MODE_MULTITHREAD | MODE_DIFFUSE);
	//_render_mode = (mode) (MODE_MULTITHREAD | MODE_SPECULAR);
	
	Raytracer raytracer;

	int width = 800; 
	int height = 800; 

	if (argc == 3) {
		width = atoi(argv[1]);
		height = atoi(argv[2]);
	}

	// Camera parameters.
	Point3D eye(0, 0, 1);
	Vector3D view(0, 0, -1);
	Vector3D up(0, 1, 0);
	double fov = 60;

	/********************************************************************
		SCENE 1: 
			"Room" with 2 mirror spheres and 2 solid sphere and 2 lights.
	*********************************************************************/
	/*
	_render_mode = (mode)(MODE_FULL_PHONG | MODE_MULTITHREAD | MODE_SHADOW |MODE_REFLECT);

	raytracer.addLightSource( new PointLight(Point3D(0, 3.9, 1), Colour(0.5,0.5,0.5)));
	raytracer.addLightSource( new PointLight(Point3D(0, 0, 1), Colour(0.5,0.5,0.5)));
	raytracer.addLightSource( new PointLight(Point3D(0, 3.9, 1), Colour(0.5,0.5,0.5)));

    SceneDagNode* sphere1 = raytracer.addObject( new UnitSphere(), &mat_mirror );
	SceneDagNode* sphere2 = raytracer.addObject( new UnitSphere(), &mat_copper );
	SceneDagNode* sphere3 = raytracer.addObject( new UnitSphere(), &mat_mirror );
	SceneDagNode* sphere4 = raytracer.addObject( new UnitSphere(), &mat_yellow );

    SceneDagNode* plane_back = raytracer.addObject( new UnitSquare(), &mat_blue );
	SceneDagNode* plane_bottom = raytracer.addObject( new UnitSquare(), &mat_green);
	SceneDagNode* plane_left = raytracer.addObject( new UnitSquare(), &mat_red);
	SceneDagNode* plane_top = raytracer.addObject( new UnitSquare(), &mat_gold);
	SceneDagNode* plane_right = raytracer.addObject( new UnitSquare(), &mat_chrome);
	//SceneDagNode* plane_behind = raytracer.addObject( new UnitSquare(), &mat_jade);

	raytracer.translate(sphere1, Vector3D(-1.5, -0.5, -4.5));

	double factor1[3] = {0.3, 0.3, 0.3};
	raytracer.translate(sphere2, Vector3D(2, 0, -4.5));
	raytracer.scale(sphere2, Point3D(0,0,0), factor1);

	raytracer.translate(sphere4, Vector3D(1.5, 1.5, -1));
	raytracer.scale(sphere4, Point3D(0,0,0), factor1);

	double factor3[3] = {0.5, 0.5, 0.5};
	raytracer.translate(sphere3, Vector3D(0, -0.5, -5));
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

	raytracer.translate(plane_top, Vector3D(0, 4, -3));        
    raytracer.scale(plane_top, Point3D(0, 0, 0), factor2);
	raytracer.rotate(plane_top, 'x', 90);

	raytracer.translate(plane_right, Vector3D(4, 0, -3));        
    raytracer.scale(plane_right, Point3D(0, 0, 0), factor2);
	raytracer.rotate(plane_right, 'y', -90);
	*/

	// SCENE 2:
	double factor0[3] = {0.05, 0.05, 0.05};
	double factor1[3] = {0.5, 0.5, 0.5};
	double factor2[3] = { 8.0, 8.0, 8.0 };
	double factor3[3] = {2, 2, 2};

	raytracer.addLightSource( new PointLight(Point3D(0, 0, 1), Colour(0.5,0.5,0.5)));
	raytracer.addLightSource( new PointLight(Point3D(0, 3.9, 1), Colour(0.5,0.5,0.5)));

	SceneDagNode* plane_back = raytracer.addObject( new UnitSquare(), &mat_blue );

	SceneDagNode* plane_bottom = raytracer.addObject( new UnitSquare(), &mat_green);
	SceneDagNode* plane_left = raytracer.addObject( new UnitSquare(), &mat_red);
	SceneDagNode* plane_top = raytracer.addObject( new UnitSquare(), &mat_gold);
	SceneDagNode* plane_right = raytracer.addObject( new UnitSquare(), &mat_chrome);
	
    raytracer.translate(plane_back, Vector3D(0, 0, -7));        
    raytracer.scale(plane_back, Point3D(0, 0, 0), factor2);

	raytracer.translate(plane_bottom, Vector3D(0, -2, -3));        
    raytracer.scale(plane_bottom, Point3D(0, 0, 0), factor2);
	raytracer.rotate(plane_bottom, 'x', -90);

	raytracer.translate(plane_left, Vector3D(-4, 0, -3));        
    raytracer.scale(plane_left, Point3D(0, 0, 0), factor2);
	raytracer.rotate(plane_left, 'y', 90);

	raytracer.translate(plane_top, Vector3D(0, 4, -3));        
    raytracer.scale(plane_top, Point3D(0, 0, 0), factor2);
	raytracer.rotate(plane_top, 'x', 90);

	raytracer.translate(plane_right, Vector3D(4, 0, -3));        
    raytracer.scale(plane_right, Point3D(0, 0, 0), factor2);
	raytracer.rotate(plane_right, 'y', -90);

	SceneDagNode* sphere2 = raytracer.addObject( new UnitSphere(), &mat_yellow );
	raytracer.translate(sphere2, Vector3D(2, 0, -2.5));
	raytracer.scale(sphere2, Point3D(0,0,0), factor1);
		
	SceneDagNode* sphere3 = raytracer.addObject( new UnitSphere(), &mat_jade );
	raytracer.translate(sphere3, Vector3D(1, 2, -3));
	raytracer.scale(sphere3, Point3D(0,0,0), factor1);

	SceneDagNode* hyper = raytracer.addObject( new Hyperboloid(2.5), &mat_mirror);
	raytracer.translate(hyper, Vector3D(-0.5, 0, -4));
	raytracer.scale(hyper, Point3D(0,0,0), factor1);
	raytracer.rotate(hyper, 'x', 70);
	
	SceneDagNode* sphere1 = raytracer.addObject( new UnitSphere(), &mat_copper);
	raytracer.translate(sphere1, Vector3D(-2, 1, -2.7));
	raytracer.scale(sphere1, Point3D(0,0,0), factor1);


	// Render the scene, feel free to make the image smaller for
	// testing purposes.	
	raytracer.render(width, height, eye, view, up, fov, "view1.bmp");
	
	// Render it from a different point of view.
	Point3D eye2(2, 2, 1);
	Vector3D view2(-2, -2, -6);
	raytracer.render(width, height, eye2, view2, up, fov, "view2.bmp");
	
	// Calculate time taken
	duration = (std::clock() - start_time) / 1000.0;
	std::cout << "Finished rendering in " << duration << " sec." << '\n' << '\n';

	return 0;
}
