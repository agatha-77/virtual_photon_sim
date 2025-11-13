#ifndef POINT_DISTRIBUTION
#define POINT_DISTRIBUTION

#include <cmath>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>

#include "phys_const.hpp"

// Funções de Bessel resumidas para deixar mais legível o resto do codigo
double K0(double b_ARG)
{
	return gsl_sf_bessel_K0(b_ARG);
}

double K1(double b_ARG)
{
	return gsl_sf_bessel_K1(b_ARG);
}

// Parâmetros do íon (massa e energia de centro de massa)
struct ion_params
{
	int mass_num;
	int atomic_num;
	double energy_CMS; // Raiz de s sobre 2
};

double epa_fraction_flux(double arg_fraction, void* params)
{
	struct ion_params* cast_params = (struct ion_params*) params;

	const double IMP_PAR_MIN = 2.0* 5.916 * pow(cast_params->mass_num,1.0/3.0) *
		(1.0 / GSL_CONST_NUM_GIGA);
	double mass = cast_params->atomic_num * PROTON_MASS +
		(cast_params->mass_num - cast_params->atomic_num) * NEUTRON_MASS;
	double gamma = cast_params->energy_CMS / mass;

	double beta = sqrt(1.0 - (1.0 / (gamma*gamma)));
	double bessel_arg = mass * arg_fraction * IMP_PAR_MIN / (gamma*beta);

	double frontal_mult = FINE_STRUCT_CONST
		* cast_params->atomic_num *cast_params->atomic_num /
		(PI * arg_fraction);

	return frontal_mult * ( 2.0 * bessel_arg * K0(bessel_arg) * K1(bessel_arg) -
			(bessel_arg*bessel_arg) * ( K1(bessel_arg)*K1(bessel_arg) -
			K0(bessel_arg)*K0(bessel_arg) ));
}

// Parâmetros do elétron (energia de centro de massa e massa do sistema produzido)
struct electron_params
{
	double energy_CMS;
	double produced_system_mass;
};

double electron_fraction_flux(double arg_fraction, void* params)
{
	struct electron_params* cast_params = (struct electron_params*) params;

	double frontal_mult = FINE_STRUCT_CONST / (2*PI);
	double energy = cast_params->energy_CMS / 2.0;
	double max_virt_sqrd = cast_params->produced_system_mass *
		cast_params->produced_system_mass;
	double min_virt_sqrd = ELECTRON_MASS*ELECTRON_MASS / (1-arg_fraction);

	return frontal_mult * ( (1 + (1-arg_fraction)*(1-arg_fraction) )/(arg_fraction) 
		   *log(max_virt_sqrd/min_virt_sqrd) + 2*ELECTRON_MASS*ELECTRON_MASS *
		   arg_fraction * ((1/max_virt_sqrd) - (1/min_virt_sqrd)) );

}

#endif
