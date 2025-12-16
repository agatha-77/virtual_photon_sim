#include <cmath>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>

#include "phys_const.hpp"

/*
 * **************************************************************
 * Parte de fatores de forma
 * **************************************************************
 */

struct wood_saxon_par
{
	int mass_num;
	int atomic_num;
};

double wood_saxon_distribution(double radius, void* params)
{
	struct wood_saxon_par* cast_params = (struct wood_saxon_par*) params;
	const double ATOM_RADIUS = 1.07e-15 * METRE_TO_EV *
		pow(cast_params->mass_num,1.0/3.0);
	const double A_PAR = 0.54e-15 * METRE_TO_EV;
	const double RHO_0 = 0.17e15 * (((double)cast_params->mass_num) / 
		((double)cast_params->atomic_num)) * 1.0 / pow(METRE_TO_EV, 3);

	return RHO_0 / (1 + exp( (radius - ATOM_RADIUS) / A_PAR ) );
}

double form_factor_wood_saxon(double q, void* params)
{
	struct wood_saxon_par* cast_params = (struct wood_saxon_par*) params;

	const double ATOM_RADIUS = 1.07e-15 * METRE_TO_EV *
		pow(cast_params->mass_num,1.0/3.0);
	const double A_PAR = 0.54e-15 * METRE_TO_EV;
	const double RHO_0 = 0.17e15 * METRE_TO_EV * ((double)cast_params->mass_num) / 
		((double)cast_params->atomic_num);

	return (4*PI / (cast_params->mass_num * q*q*q)) * RHO_0 
		* (sin(q*ATOM_RADIUS) - q*ATOM_RADIUS* cos(q*ATOM_RADIUS)) 
		/ (1 + q*q * A_PAR*A_PAR);
}

double form_factor_monopole(double q, void* params)
{
	const double LAMBDA = 0.088e9;

	return (LAMBDA*LAMBDA) / (LAMBDA*LAMBDA + q*q);
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
