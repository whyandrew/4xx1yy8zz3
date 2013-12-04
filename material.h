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

Material mat_gold("Gold", Colour(0.24725, 0.1995, 0.0745), 
				  Colour(0.75164, 0.60648, 0.22648), 
				  Colour(0.628281, 0.555802, 0.366065), 
				  51.2, 0.6 );

Material mat_jade("Jade", Colour(0.135, 0.2225, 0.1575), 
				  Colour(0.54, 0.89, 0.63), 
				  Colour(0.316228, 0.316228, 0.316228), 
				  12.8, 0.3 );

Material mat_copper("Copper", Colour(0.19125, 0.0735, 0.0225), 
					Colour(0.7038, 0.27048, 0.0828), 
					Colour(0.256777, 0.137622, 0.086014), 
					12.8, 0.5 );

Material mat_chrome("Chrome", Colour(0.25, 0.25, 0.25), 
					Colour(0.4,0.4,0.4), 
					Colour(0.774597,0.774597,0.774597), 
					128*0.6, 0.3 );

Material mat_mirror("Mirror", Colour(0,0,0), 
					Colour(0.00, 0.0, 0.0), 
					Colour(1,1,1), 
					128, 1);

Material mat_mirror_regular("Mirror_regular", Colour(0.0,0.0,0.0), 
					Colour(0.05, 0.05, 0.05), 
					Colour(1.0,1.0,1.0), 
					128, 1.0);

Material mat_red("Red", Colour(0.3,0.1,0.1), 
					Colour(0.4,0.1,0.1), 
					Colour(0.7,0.2,0.2), 
					128*0.6, 0.4);

Material mat_green( "Green", Colour(0.1,0.3,0.1), 
					Colour(0.1,0.4,0.1), 
					Colour(0.2,0.7,0.2), 
					128*0.6, 0.4);

Material mat_blue( "Blue", Colour(0.1,0.1,0.3), 
					Colour(0.1,0.1,0.4), 
					Colour(0.2,0.2,0.7), 
					128*0.6, 0.4);

Material mat_yellow( "Yellow", Colour(0.3,0.3,0.1), 
					Colour(0.6,0.6,0.1), 
					Colour(0.6,0.6,0.2), 
					128*0.6, 0.4);

Material mat_hotpink( "Hotpink", Colour(0.3,0.1,0.2), 
					Colour(0.6,0.1,0.5), 
					Colour(0.6,0.2,0.5), 
					128*0.6, 0.2);

Material mat_purple( "Purple", Colour(0.25,0.1,0.4), 
					Colour(0.55,0.1,0.7), 
					Colour(0.55,0.2,0.7), 
					128*0.6, 0.2);

Material mat_glass( "Glass", Colour(0.1,0.1,0.1), 
					Colour(0, 0, 0), 
					Colour(0, 0, 0), 
					128*0.6, 0.2, 1.6);

Material mat_diamond( "Diamond", Colour(0.1,0.1,0.1), 
					Colour(0, 0, 0), 
					Colour(0, 0, 0), 
					128*0.6, 0.2, 2.419);

Material mat_water( "Water", Colour(0.1,0.1,0.1), 
					Colour(0, 0, 0), 
					Colour(0, 0, 0), 
					128*0.6, 0.2, 1.3330);

Material mat_light( "Light", Colour(1.0,1.0,1.0), 
					Colour(1.0,1.0,1.0), 
					Colour(1.0,1.0,1.0), 
					128*0.6, 0.2);

// For texture, the colour components are multipled to texture
// Set colors to 0.0 if not to be affected by that light component

Material texture_earth( "Texture Earth", "earth.bmp",
					   Colour(0.05, 0.05, 0.05),
					   Colour(1.0, 1.0 , 1.0),
					   Colour(0.0, 0.0 , 0.0),
					   128*0.2, 0.0);

Material texture_moon( "Texture Moon", "moon.bmp",
					   Colour(0.15, 0.15, 0.15),
					   Colour(0.8, 0.8 , 0.8),
					   Colour(0.0, 0.0 , 0.0),
					   128*0.2, 0.0);

Material texture_galaxy( "Texture Galaxy", "galaxy.bmp",
					   Colour(1.0, 1.0, 1.0),
					   Colour(0.0, 0.0 , 0.0),
					   Colour(0.0, 0.0 , 0.0),
					   128*0.2, 0.0);
/*
Material texture_sun( "Texture Sun", "sun.bmp",
					   Colour(1.0, 1.0, 1.0),
					   Colour(1.0, 1.0, 1.0),
					   Colour(0.0, 0.0 , 0.0),
					   128*0.2, 0.0);
					   */
Material texture_jupiter( "Texture Jupiter", "jupiter.bmp",
					   Colour(0.05, 0.05, 0.05),
					   Colour(0.8, 0.8 , 0.8),
					   Colour(0.0, 0.0 , 0.0),
					   128*0.2, 0.0);

Material texture_neptune( "Texture neptune", "neptune.bmp",
					   Colour(0.05, 0.05, 0.05),
					   Colour(0.8, 0.8 , 0.8),
					   Colour(0.0, 0.0 , 0.0),
					   128*0.2, 0.0);

Material texture_hardwood( "Texture Hardwood", "hardwood.bmp",
						Colour(0.1,0.1,0.1),
						Colour(1.0, 1.0, 1.0),
						Colour(1.0, 1.0, 1.0),
						128*0.4, 0.7);