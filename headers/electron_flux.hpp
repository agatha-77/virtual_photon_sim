/*
 * Biblioteca para fluxo de fótons do elétron.
 *
 */
#ifndef ELECTRON_FLUX
#define ELECTRON_FLUX

#include<cmath>

#include<gsl/gsl_math.h>

// Parâmetros do elétron (energia de centro de massa e massa do sistema produzido)
struct electron_params
{
	double energy_CMS;
	double produced_system_mass;
};

// Fluxo de fótons do elétron (Budnev, 1974)
double electron_photon_flux(double freq, void* params)
{
	struct electron_params* cast_params = (struct electron_params*) params;

	double frontal_mult = FINE_STRUCT_CONST / PI;
	double energy = cast_params->energy_CMS / 2.0;

	double max_virt= cast_params->produced_system_mass *
		cast_params->produced_system_mass;
	double min_virt= ELECTRON_MASS * ELECTRON_MASS * freq * freq
		/ (energy* (energy - freq));

	return frontal_mult * (1.0 - (freq/energy) + (freq*freq / energy) ) *
		log(max_virt*max_virt/(min_virt*min_virt)) + (1.0 - freq / energy) * 
		( 1.0 - (min_virt*min_virt)/(max_virt*max_virt) );
}

#endif
