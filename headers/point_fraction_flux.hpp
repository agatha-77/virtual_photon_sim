#ifndef POINT_DISTRIBUTION
#define POINT_DISTRIBUTION

#include <cmath>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>

#include "phys_const.hpp"
#include "point_like_charge.hpp"



double EPA_fraction_flux(double arg_fraction, void* params)
{
	struct ion_params* cast_params = (struct ion_params*) params;

	const double IMP_PAR_MIN = 2.0* 5.916 * pow(cast_params->mass_num,1.0/3.0) *
		(1.0 / GSL_CONST_NUM_GIGA);
	double mass = cast_params->atomic_num * PROTON_MASS +
		(cast_params->mass_num - cast_params->atomic_num) * NEUTRON_MASS;
	double gamma = cast_params->energy_CMS / mass;

	double beta = sqrt(1.0 - 1.0 / (gamma*gamma));
	double bessel_arg = arg_fraction * IMP_PAR_MIN * cast_params->energy_CMS
		/ (gamma * beta);

	double frontal_mult = 2 * FINE_STRUCT_CONST
		* cast_params->atomic_num *cast_params->atomic_num /
		(PI*(beta*beta) * arg_fraction * cast_params->energy_CMS*cast_params->energy_CMS);

	return frontal_mult * ( bessel_arg * K0(bessel_arg) * K1(bessel_arg) -
			(beta*beta * bessel_arg*bessel_arg) * ( K1(bessel_arg)*K1(bessel_arg) -
			K0(bessel_arg)*K0(bessel_arg) ) / 2.0);
}

#endif
