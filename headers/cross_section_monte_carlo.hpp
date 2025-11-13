/*
 * Biblioteca para calculo das secoes de choque usando
 * integracao de monte carlo.
 * 
 * Segue a mesma estrutura da outra biblioteca.
 */

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>

#include "phys_const.hpp"
#include "point_like_charge.hpp"

struct dilepton_params
{
	double beam_energy;
	double lepton_mass;
};

double fundamental_CS_dilepton(double arg, double mass, void* params)
{
	long double w = arg;

	return (4.0*PI * FINE_STRUCT_CONST*FINE_STRUCT_CONST / (w*w))
		* (2.0* log( (w / (2*mass)) * 
		(1.0 + sqrt(1.0 - ( (4*mass*mass) / (w*w))) ))
		* (1.0 + (4.0*mass*mass*w*w - 8.0*mass*mass*mass*mass)/(w*w*w*w) ) 
		- sqrt(1.0 - (4.0*mass*mass) /(w*w) ) * (1.0 + (4.0*mass*mass)/(w*w)) );
}


double integrand_EPA_dilepton(double* var, size_t dim, void* params)
{
	struct dilepton_params* cast_params = (struct dilepton_params*) params;
	double cms_energy = cast_params->beam_energy;
	double mass = cast_params->lepton_mass;
	double dummy = 0.0;

	double var1 = var[0];
	double var2 = var[1];

	struct ion_params gold179_params;
	gold179_params.atomic_num = 79;
	gold179_params.mass_num = 179;
	gold179_params.energy_CMS = cms_energy;

	struct ion_params pb208_params;
	pb208_params.atomic_num = 82;
	pb208_params.mass_num = 208;
	pb208_params.energy_CMS = cms_energy;

	return (epa_photon_flux(var1*cms_energy, &pb208_params)/var1) * 
		(epa_photon_flux(var2*cms_energy, &pb208_params)/var2) *
		fundamental_CS_dilepton(sqrt(4*var1*var2*cms_energy*cms_energy), mass, &dummy);
}


double integral_monte_vegas( double (*fptr) (double* arg, size_t dim, void* params),
		double* x_i, double* x_f, const size_t DIMENSION, double* err, void* params)
{
	double res;
	size_t calls = 1000;

	const gsl_rng_type* T;
	gsl_rng* r;

	gsl_monte_function INTEGR_FUNCTION = {fptr, DIMENSION, params};

	gsl_rng_env_setup();

	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	gsl_monte_vegas_state* s = gsl_monte_vegas_alloc(DIMENSION);
	gsl_monte_vegas_integrate(&INTEGR_FUNCTION,
			x_i, x_f,
			DIMENSION, calls,
			r, s,
			&res, err);

	std::cout << "-----------------------------------------\n";
	std::cout << "Calculando integral:\n";
	std::cout << "resultado = " << res << "\nerro = " << *err << "\n";
	std::cout << "-----------------------------------------\n";

	do {
		gsl_monte_vegas_integrate(&INTEGR_FUNCTION,
				x_i, x_f,
				DIMENSION, calls/5,
				r, s,
				&res, err);
		//std::cout << res << "\t" << *err << "\n";

	} while(fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5);

	gsl_monte_vegas_free(s);
	gsl_rng_free(r);

	return res;
}

double dilepton_TCS_EPA_monte_vegas(double beam_energy, double mass, double* err)
{
	double alpha = 0.0;
	double result;
	double x_i[2];
	double x_f[2];

	struct dilepton_params lepton_params;
	lepton_params.beam_energy = beam_energy;
	lepton_params.lepton_mass = mass;

	gsl_set_error_handler_off();

	x_i[0] = mass / beam_energy; x_f[0] = 1.0;
	x_i[1] = mass*mass / (beam_energy*beam_energy * x_i[0]);
	x_f[1] = 1.0;

	result = integral_monte_vegas(&integrand_EPA_dilepton,
			x_i, x_f, 
			2, err, &lepton_params);

	return result;

}
