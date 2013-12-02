/*********************************************
	Andrew Wong
	2013 11 24

	Materials values from:
	http://devernay.free.fr/cours/opengl/materials.html


**********************************************/

#include "util.h"

// Define refractive index
#define RI_AIR  1.000277f

// Define material properties
// Material mat_name( Color(ambient), Color(diffuse), Color(specular), 
//					spec_exp, reflectivity, [refractive index]);

// reflectivity of 1.0 wil be treated as perfectly reflective

// If material is non-opaque, with refractive index > 1, 
// transparency (factor) of material is taken as (1.0 - (amb+diffuse+spec))
// e.g. water normally has ~80% transparency

Material mat_gold("gold", Colour(0.24725, 0.1995, 0.0745), 
				  Colour(0.75164, 0.60648, 0.22648), 
				  Colour(0.628281, 0.555802, 0.366065), 
				  51.2, 0.6 );

Material mat_jade("jade", Colour(0.135, 0.2225, 0.1575), 
				  Colour(0.54, 0.89, 0.63), 
				  Colour(0.316228, 0.316228, 0.316228), 
				  12.8, 0.3 );

Material mat_copper("copper", Colour(0.19125, 0.0735, 0.0225), 
					Colour(0.7038, 0.27048, 0.0828), 
					Colour(0.256777, 0.137622, 0.086014), 
					12.8, 0.5 );

Material mat_chrome("chrome", Colour(0.25, 0.25, 0.25), 
					Colour(0.4,0.4,0.4), 
					Colour(0.774597,0.774597,0.774597), 
					128*0.6, 0.3 );

Material mat_mirror("mirror", Colour(0,0,0), 
					Colour(0.00, 0.0, 0.0), 
					Colour(1,1,1), 
					128, 1);

Material mat_mirror_regular("mirror_regular", Colour(0,0,0), 
					Colour(0.05, 0.05, 0.05), 
					Colour(1,1,1), 
					128, 1);

Material mat_red("red", Colour(0.3,0.1,0.1), 
					Colour(0.4,0.1,0.1), 
					Colour(0.7,0.2,0.2), 
					128*0.6, 0.4);

Material mat_green( "green", Colour(0.1,0.3,0.1), 
					Colour(0.1,0.4,0.1), 
					Colour(0.2,0.7,0.2), 
					128*0.6, 0.4);

Material mat_blue( "blue", Colour(0.1,0.1,0.3), 
					Colour(0.1,0.1,0.4), 
					Colour(0.2,0.2,0.7), 
					128*0.6, 0.4);

Material mat_yellow( "yellow", Colour(0.3,0.3,0.1), 
					Colour(0.6,0.6,0.1), 
					Colour(0.6,0.6,0.2), 
					128*0.6, 0.4);

Material mat_hotpink( "hotpink", Colour(0.3,0.1,0.2), 
					Colour(0.6,0.1,0.5), 
					Colour(0.6,0.2,0.5), 
					128*0.6, 0.2);

Material mat_purple( "purple", Colour(0.25,0.1,0.4), 
					Colour(0.55,0.1,0.7), 
					Colour(0.55,0.2,0.7), 
					128*0.6, 0.2);

Material mat_glass( "glass", Colour(0.1,0.1,0.1), 
					Colour(0, 0, 0), 
					Colour(0, 0, 0), 
					128*0.6, 0.2, 1.6);

Material mat_diamond( "diamond", Colour(0.1,0.1,0.1), 
					Colour(0, 0, 0), 
					Colour(0, 0, 0), 
					128*0.6, 0.2, 2.419);