
#include <cmath>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>

#include "phys_const.hpp"

/*
 * **************************************************************
 * Parte de fatores de forma
 * **************************************************************
 */

double form_factor_wood_saxon(double q, void* params)
{
	double 

}

/*
 * **************************************************************
 * Parte dos fluxos de fótons
 * **************************************************************
 */

struct flux_paramater
{
	double alpha;
	double lambda;
	double atom_radius;
	double (*form_factor_ptr) (double q, void* params);
};
