#ifndef CROSS_SECTION
#define CROSS_SECTION

#include <cmath>
#include<gsl/gsl_math.h>
#include "phys_const.hpp"
#include "point_like_charge.hpp"

long double fundamental_CS(long double arg)
{
	long double w = arg;
	long double mass = MUON_MASS;

	return (4*PI * FINE_STRUCT_CONST*FINE_STRUCT_CONST / w) * (2* log( (w / (2*mass)) * 
			(1 + sqrt(1 - (4*mass*mass / (w*w))) )) * (1 + (4*mass*mass*w*w - 8*mass*mass)/(w*w*w*w) ) 
			- sqrt(1 - (4*mass*mass) /(w*w) ) * (1 + 4*mass*mass/(w*w)) );
}



#endif
