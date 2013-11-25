/*********************************************
	Andrew Wong
	2013 11 24

	Materials values from:
	http://devernay.free.fr/cours/opengl/materials.html


**********************************************/

#include "util.h"

// Define material properties
// Material mat_name( Color(ambient), Color(diffuse), Color(specular), 
//					spec_exp, reflectivity);

// reflectivity of 1.0 wil be treated as perfectly reflective

Material mat_gold( Colour(0.24725, 0.1995, 0.0745), 
				  Colour(0.75164, 0.60648, 0.22648), 
				  Colour(0.628281, 0.555802, 0.366065), 
				  51.2, 0.6 );

Material mat_jade( Colour(0.135, 0.2225, 0.1575), 
				  Colour(0.54, 0.89, 0.63), 
				  Colour(0.316228, 0.316228, 0.316228), 
				  12.8, 0.3 );

Material mat_copper( Colour(0.19125, 0.0735, 0.0225), 
					Colour(0.7038, 0.27048, 0.0828), 
					Colour(0.256777, 0.137622, 0.086014), 
					12.8, 0.5 );

Material mat_chrome( Colour(0.25, 0.25, 0.25), 
					Colour(0.4,0.4,0.4), 
					Colour(0.774597,0.774597,0.774597), 
					128*0.6, 0.7 );

Material mat_glass( Colour(0.1,0.1,0.1), 
					Colour(0.1,0.1,0.1), 
					Colour(0.1,0.1,0.1), 
					128, 1);



Material mat_red( Colour(0.3,0.1,0.1), 
					Colour(0.4,0.1,0.1), 
					Colour(0.7,0.2,0.2), 
					128*0.6, 0.3);

Material mat_green( Colour(0.1,0.3,0.1), 
					Colour(0.1,0.4,0.1), 
					Colour(0.2,0.7,0.2), 
					128*0.6, 0.3);

Material mat_blue( Colour(0.1,0.1,0.3), 
					Colour(0.1,0.1,0.4), 
					Colour(0.2,0.2,0.7), 
					128*0.6, 0.3);

Material mat_yellow( Colour(0.3,0.3,0.1), 
					Colour(0.6,0.6,0.1), 
					Colour(0.6,0.6,0.2), 
					128*0.6, 0.3);