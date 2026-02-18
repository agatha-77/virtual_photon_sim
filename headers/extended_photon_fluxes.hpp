/*
 * Fluxos de fótons para distribuições de carga extendidas. 
 * *************************************************************************
 *
 * Algumas das integrais não podem ser postas em forma multidimensional tão
 * facilmente, logo a implementação de monte carlo é mais difícil.
 *
 */
#include <cmath>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>

#include "phys_const.hpp"

double heaviside_func(double x)
{
	if (x < 0.0) {
		return 0.0;
	} else if (x >= 0.0) {
		return 1.0;
	}
}

/*
 * **************************************************************
 * Parte de fatores de forma
 * **************************************************************
 */

struct form_factor_par
{
	int mass_num;
	int atomic_num;
};

double wood_saxon_distribution(double radius, void* params)
{
	struct form_factor_par* cast_params = (struct form_factor_par*) params;
	const double ATOM_RADIUS = 1.07e-15 * METRE_TO_EV *
		pow(cast_params->mass_num,1.0/3.0);
	const double A_PAR = 0.54e-15 * METRE_TO_EV;

	const double RHO_0 = (cast_params->atomic_num * FINE_STRUCT_CONST / sqrt(4*PI))
		/ (1.803085 + ATOM_RADIUS * ATOM_RADIUS * 0.69314 / A_PAR + ATOM_RADIUS
		* PI*PI / 6.0);

	return RHO_0 / (1 + exp( (radius - ATOM_RADIUS) / A_PAR ) );
}

double form_factor_wood_saxon(double q, void* params)
{
	struct form_factor_par* cast_params = (struct form_factor_par*) params;

	const double ATOM_RADIUS = 1.07e-15 * METRE_TO_EV *
		pow(cast_params->mass_num,1.0/3.0);
	const double A_PAR = 0.54e-15 * METRE_TO_EV;

	const double RHO_0 = 3 * cast_params->mass_num /
		(4.0*PI*ATOM_RADIUS*ATOM_RADIUS*ATOM_RADIUS);

	return ((4*PI * RHO_0  / (cast_params->mass_num * q*q*q * (1 + q*q*A_PAR*A_PAR))) 
		* (sin(q*ATOM_RADIUS) - q*ATOM_RADIUS* cos(q*ATOM_RADIUS)) );
}

double form_factor_monopole(double q, void* params)
{
	struct form_factor_par* cast_params = (struct form_factor_par*) params;
	const double LAMBDA = 0.088e9;

	return (LAMBDA*LAMBDA) / (LAMBDA*LAMBDA + q*q);
}

/*
 * **************************************************************
 * Parte dos fluxos de fótons
 * **************************************************************
 */

double J1(double arg)
{
	return gsl_sf_bessel_J1(arg);
}

struct flux_parameter
{
	int mass_num;
	int atomic_num;
	double beam_energy;
	double (*form_factor_ptr) (double q, void* params);
	int npontos;
	double freq_0;
	double freq_f;
};

struct photon_integrand_params
{
	int mass_num;
	int atomic_num;
	double freq;
	double ion_mass;
	double beam_energy;
	double impact_par;
	double (*form_factor_ptr) (double q, void* params);
};

// Integrando da funcao abaixo
double photon_flux_integrand(double u, void* params)
{
	struct photon_integrand_params* cast_params
		= (struct photon_integrand_params*) params;

	struct form_factor_par form_par;
	form_par.mass_num = cast_params->mass_num;
	form_par.atomic_num = cast_params->atomic_num;

	double gamma = cast_params->beam_energy / (cast_params->ion_mass);
	double beta = sqrt(1.0 - 1.0 / (gamma*gamma));
	double impact_par = cast_params->impact_par;
	double xi_arg = cast_params->freq * impact_par / (gamma * beta);

	return u*u * J1(u) *
		cast_params->form_factor_ptr( (u*u + xi_arg*xi_arg) /
				impact_par*impact_par , &form_par) / (u*u + xi_arg*xi_arg);
}

// Função para calcular N(w,b).
// Fluxo dependente do parametro de impacto
double extended_photon_flux(double freq, double impact_par, void* params)
{
	struct flux_parameter* cast_params = (struct flux_parameter*) params;
	int mass_num = cast_params->mass_num;
	int atomic_num = cast_params->atomic_num;
	double beam_energy = cast_params->beam_energy;

	double mass = cast_params->atomic_num * PROTON_MASS + 
		(cast_params->mass_num - cast_params->atomic_num)*NEUTRON_MASS;
	double gamma = beam_energy / mass;
	double beta = sqrt(1.0 - 1.0 / (gamma*gamma));

	struct photon_integrand_params photon_params;
	photon_params.mass_num = mass_num;
	photon_params.atomic_num = atomic_num;
	photon_params.freq = freq;
	photon_params.ion_mass = mass;
	photon_params.impact_par = impact_par;
	photon_params.beam_energy = beam_energy;
	photon_params.form_factor_ptr = cast_params->form_factor_ptr;

	int size = 5000;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(size);
	gsl_set_error_handler_off();

	gsl_function integr_func;
	integr_func.function = &photon_flux_integrand;
	integr_func.params = &photon_params;

	double result, err;

	// Integral sobre intervalo semi infinito [0:infty)
	gsl_integration_qagiu(&integr_func, 0.0, 
			1e-6, 1e-3,
			size, w,
			&result, &err);

	gsl_integration_workspace_free(w);

	return ( (atomic_num * atomic_num * FINE_STRUCT_CONST)
			/ (PI*PI * beta*beta * freq * impact_par * impact_par) )
		* fabs(result*result);
}

/* Integrando da funcao abaixo */
double total_photon_flux_integrand(double* arg, size_t dim, void* params)
{
	if (dim != 3) return 0.0;
	struct photon_integrand_params* cast_params 
		= (struct photon_integrand_params*) params;

	struct form_factor_par form_par;
	form_par.mass_num = cast_params->mass_num;
	form_par.atomic_num = cast_params-> atomic_num;

	double gamma = cast_params->beam_energy / (cast_params->ion_mass);
	double beta = sqrt(1.0 - 1.0 / (gamma*gamma));
	double impact_par_min = cast_params->impact_par;

	double xi_arg = cast_params->freq * arg[0] / (gamma*beta);

	return heaviside_func(arg[0] - impact_par_min) * (1.0 / arg[0])
		* arg[1]*arg[1] * arg[2]*arg[2] 
		* cast_params->form_factor_ptr(
				(xi_arg*xi_arg * arg[1]*arg[1]) / (arg[0]*arg[0]),
				&form_par)
		* cast_params->form_factor_ptr(
				(xi_arg*xi_arg * arg[1]*arg[1]) / (arg[0]*arg[0]),
				&form_par)
		* cast_params->form_factor_ptr(
				(xi_arg*xi_arg * arg[2]*arg[2]) / (arg[0]*arg[0]),
				&form_par)
		* cast_params->form_factor_ptr(
				(xi_arg*xi_arg * arg[2]*arg[2]) / (arg[0]*arg[0]),
			&form_par)
		* J1(arg[1]) * J1(arg[2])
		/ (xi_arg*xi_arg * (xi_arg*xi_arg + arg[1]*arg[1] + arg[2]*arg[2]) 
				+ arg[1]*arg[1] * arg[2]*arg[2]);
}

void calc_extended_total_photon_flux(double* freq, double* result,
		double* err, void* params)
{
	struct flux_parameter* cast_params = (struct flux_parameter*) params;
	int mass_num = cast_params->mass_num;
	int atomic_num = cast_params->atomic_num;
	double beam_energy = cast_params->beam_energy;
	int npontos = cast_params->npontos;
	double step = fabs(cast_params->freq_0 - cast_params->freq_f) / (double) npontos;

	double nuclear_radius = 1.21e-15 * pow(atomic_num, 3) * METRE_TO_EV;

	double mass = cast_params->atomic_num * PROTON_MASS + 
		(cast_params->mass_num - cast_params->atomic_num)*NEUTRON_MASS;
	double gamma = beam_energy / mass;
	double beta = sqrt(1.0 - 1.0 / (gamma*gamma));

	/* struct que é passado ao integrando */
	struct photon_integrand_params photon_params;
	photon_params.mass_num = mass_num;
	photon_params.atomic_num = atomic_num;
	photon_params.ion_mass = mass;
	photon_params.impact_par = 2*nuclear_radius;
	photon_params.beam_energy = beam_energy;
	photon_params.form_factor_ptr = cast_params->form_factor_ptr;

	 // chamadas da função de integração, quanto mais alto mais acurado
	int calls = 10000000;
	double integr_result;
	double integr_err;
	const int DIMENSION = 3;
	const double UPPER_BOUND = 100.0;

	gsl_monte_function F;
	F.f = &total_photon_flux_integrand;
	F.dim = DIMENSION;
	F.params = &photon_params;

	gsl_rng_env_setup();

	const gsl_rng_type* T = gsl_rng_default;
	gsl_rng* r = gsl_rng_alloc(T);

	gsl_monte_miser_state* s = gsl_monte_miser_alloc(DIMENSION);

	double x0[3] = {0.0, 0.0, 0.0};
	double xf[3] = {UPPER_BOUND, UPPER_BOUND, UPPER_BOUND};

	double var_freq = cast_params->freq_0;

	for(int i = 0; i <= npontos; i++){
		freq[i] = var_freq;
		photon_params.freq = var_freq;

		gsl_monte_miser_init(s);
		gsl_monte_miser_integrate(&F, x0, xf, DIMENSION, calls, r, s, 
				&integr_result, &integr_err);

		result[i] = atomic_num*atomic_num * FINE_STRUCT_CONST * fabs(integr_result) /
		(PI*PI*beta*beta * var_freq);
		err[i] = integr_err;

		var_freq += step;
	}

	gsl_monte_miser_free(s);
	gsl_rng_free(r);
}

