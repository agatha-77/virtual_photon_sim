
/*
 * Biblioteca para cálculo de seções de choque e quantidades relacionadas 
 * (valor máximo de seção de choque e afins).
 *
 * Unidades naturais são usadas e métodos de integração são baseados em Monte
 * Vegas.
 */

#ifndef CROSS_SECTION
#define CROSS_SECTION

#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_deriv.h>

#include "phys_const.hpp"
#include "point_like_charge.hpp"


double fundamental_CS_muon(double arg, void* params)
{
	(void)(params);
	long double w = arg;
	long double mass = MUON_MASS;

	return (4.0*PI * FINE_STRUCT_CONST*FINE_STRUCT_CONST / (w*w)) * (2.0* log(
				(w / (2*mass)) * 
			(1.0 + sqrt(1.0 - ( (4*mass*mass) / (w*w))) )) * (1.0 + (4.0*mass*mass*w*w -
					8.0*mass*mass*mass*mass)/(w*w*w*w) ) 
			- sqrt(1.0 - (4.0*mass*mass) /(w*w) ) * (1.0 + (4.0*mass*mass)/(w*w)) );
}

double fundamental_CS_electron(double arg, void* params)
{
	(void)(params);
	double w = arg;
	double mass = ELECTRON_MASS;

	return (4*PI * FINE_STRUCT_CONST*FINE_STRUCT_CONST / (w*w)) * (2* log( (w / (2*mass)) * 
			(1 + sqrt(1 - (4*mass*mass / (w*w))) )) * (1 + (4*mass*mass*w*w -
					8*mass*mass*mass*mass)/(w*w*w*w) ) 
			- sqrt(1 - (4*mass*mass) /(w*w) ) * (1 + 4*mass*mass/(w*w)) );
}


// Calcula o máximo da seção de choque
double max_value_CS(double (*fptr) (double x, void* params), double x_i, double x_f)
{
	const int NSTEPS = 1e+6;
	const double DERIV_ERR_TOL = 1e-8;
	const double ABS_ERR_TOL = 1e-6;

	double step = fabs(x_f - x_i) / (double) NSTEPS;
	double x_var = x_i;
	double result_deriv, abs_err;
	double alpha = 0.0;

	gsl_function F;

	F.function = fptr;
	F.params = &alpha;

	for(int i = 0; i < NSTEPS; i++){
		gsl_deriv_central(&F, x_var, step, &result_deriv, &abs_err);

		if(fabs(result_deriv) < DERIV_ERR_TOL && abs_err < ABS_ERR_TOL)
			return x_var;

		x_var += step;
	}

	return 0.0;
}


// ****************************************************************
// Parte de integrais
// ****************************************************************

double integral( double (*fptr) (double* arg, size_t dim, void* params),
		double* x_i, double* x_f,
		const size_t DIMENSION, double* err, double* params)
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

	std::cout << "\nCalculando integral:\n";
	std::cout << "resultado = " << res << "\terro = " << *err << "\n";
	std::cout << "-----------------------------------------";

	do {
		gsl_monte_vegas_integrate(&INTEGR_FUNCTION,
				x_i, x_f,
				DIMENSION, calls/5,
				r, s,
				&res, err);

	} while(fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5);

	gsl_monte_vegas_free(s);
	gsl_rng_free(r);

	return res;
}


double muon_integrand(double* var, size_t dim, void* params)
{
	(void) (dim);
	(void) (params);
	double dummy = 0.0;

	return (ep_num_total(var[0]) / var[0]) * (ep_num_total(var[1]) / var[1])
		* fundamental_CS_muon(var[0] + var[1], &dummy);
}

double total_muon_cross_section(double beam_energy, double* err)
{
	double x_i[2];
	double x_f[2];

	double alpha = 0.0;

	x_i[0] = (MUON_MASS / beam_energy) * (MUON_MASS / beam_energy);
	x_i[1] = 1.0;

	x_f[0] = 1.0;
	x_f[1] = (MUON_MASS / (beam_energy * x_i[0])) * (MUON_MASS / (beam_energy * x_i[0]));

	std::cout << "\n" << x_i[0] << "\t" << x_i[1] << "\n"
		<< x_f[0] << "\t" << x_f[1] << "\t"
		<< muon_integrand(x_i, 2, &alpha) << "\n";

	return integral(&muon_integrand, x_i, x_f, 2, err, &alpha);
}


#endif
