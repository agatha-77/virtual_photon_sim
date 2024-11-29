
/*
 * Biblioteca para cálculo de seções de choque e quantidades relacionadas 
 * (valor máximo de seção de choque e afins).
 *
 * Unidades naturais são usadas
 */

#ifndef CROSS_SECTION
#define CROSS_SECTION

#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>

#include "phys_const.hpp"
#include "point_like_charge.hpp"
#include "electron_flux.hpp"


double fundamental_CS_muon(double arg, void* params)
{
	(void)(params);
	long double w = arg;
	long double mass = MUON_MASS;

	return (4.0*PI * FINE_STRUCT_CONST*FINE_STRUCT_CONST / (w*w)) * (2.0* log(
				(w / (2*mass)) * 
			(1.0 + sqrt(1.0 - ( (4*mass*mass) / (w*w))) )) * (1.0 +
				(4.0*mass*mass*w*w - 8.0*mass*mass*mass*mass)/(w*w*w*w) ) 
			- sqrt(1.0 - (4.0*mass*mass) /(w*w) ) * (1.0 + (4.0*mass*mass)/(w*w)) );
}

double fundamental_CS_electron(double arg, void* params)
{
	(void)(params);
	double w = arg;
	double mass = ELECTRON_MASS;

	return (4*PI * FINE_STRUCT_CONST*FINE_STRUCT_CONST / (w*w)) *
		(2* log( (w / (2*mass)) * 
		(1 + sqrt(1 - (4*mass*mass / (w*w))) )) *
		(1 + (4*mass*mass*w*w - 8*mass*mass*mass*mass)/(w*w*w*w) ) -
		sqrt(1 - (4*mass*mass) /(w*w) ) * (1 + 4*mass*mass/(w*w)) );
}

double fundamental_CS_tau(double arg, void* params)
{
	(void)(params);
	double w = arg;
	double mass = TAU_MASS;

	return (4*PI * FINE_STRUCT_CONST*FINE_STRUCT_CONST / (w*w)) *
		(2* log( (w / (2*mass)) * 
		(1 + sqrt(1 - (4*mass*mass / (w*w))) )) *
		(1 + (4*mass*mass*w*w - 8*mass*mass*mass*mass)/(w*w*w*w) ) -
		sqrt(1 - (4*mass*mass) /(w*w) ) * (1 + 4*mass*mass/(w*w)) );
}


// Calcula o máximo da seção de choque
double max_value_CS(double (*fptr) (double x, void* params), double x_i, double x_f)
{
	const int NSTEPS = 1e+8;
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

double fundamental_CS_dilepton(double cms_energy, double lepton_mass, void* params)
{
	(void)(params);
	double w = cms_energy;
	double mass = lepton_mass;

	return (4*PI * FINE_STRUCT_CONST*FINE_STRUCT_CONST / (w*w)) *
		(2* log( (w / (2*mass)) * 
		(1 + sqrt(1 - (4*mass*mass / (w*w))) )) *
		(1 + (4*mass*mass*w*w - 8*mass*mass*mass*mass)/(w*w*w*w) ) -
		sqrt(1 - (4*mass*mass) /(w*w) ) * (1 + 4*mass*mass/(w*w)) );
}

/*
 * Parte de Integrais
 * */


struct integrand_params
{
	double beam_energy;
	double var1;
	double produced_mass;
};

struct dilepton_params
{
	double beam_energy;
	double lepton_mass;
};


// Arrumar isso aqui
double dilepton_integrand1_electron(double var2, void* params)
{
	struct integrand_params* cast_params = (struct integrand_params*) params;
	double cms_energy = cast_params->beam_energy;
	double var1 = cast_params->var1;
	double mass = cast_params->produced_mass;
	double dummy = 0.0;

	struct electron_params e_params;
	e_params.energy_CMS = cms_energy / 2.0;
	e_params.produced_system_mass = 2*mass;

	return (electron_photon_flux(var1, &e_params)/var1) *
		(electron_photon_flux(var2, &e_params)/var2) *
		fundamental_CS_dilepton(sqrt(4*var1*var2), mass, &dummy);
}

double dilepton_integrand2_electron(double var1, void* params)
{
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(5000);

	struct dilepton_params* cast_params = (struct dilepton_params*) params;
	double beam_energy = cast_params->beam_energy;
	double mass = cast_params->lepton_mass;
	double result, error;

	struct integrand_params lepton_params;
	lepton_params.beam_energy = beam_energy;
	lepton_params.var1 = var1;
	lepton_params.produced_mass = mass;

	gsl_function integr_func;
	integr_func.function = &dilepton_integrand1_electron;
	integr_func.params = &lepton_params;

	gsl_integration_qag(&integr_func, 
			mass*mass / var1, beam_energy,
			1e-10, 1e-3,
			5000, 1,
			w, &result, &error);

	gsl_integration_workspace_free(w);

	return result;
}

double dilepton_TCS_electron(double beam_energy, double lepton_mass, double* err)
{
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(5000);

	double result;
	struct dilepton_params params;
	params.beam_energy = beam_energy;
	params.lepton_mass = lepton_mass;
	
	gsl_function integr_func;
	integr_func.function = &dilepton_integrand2_electron;
	integr_func.params = &params;

	gsl_integration_qag(&integr_func,
			lepton_mass, beam_energy,
			1e-10, 1e-3,
			5000, 1,
			w, &result, err);

	gsl_integration_workspace_free(w);

	return result;
}


double muon_integrand1_EPA(double var2, void* params)
{
	struct integrand_params* cast_params = (struct integrand_params*) params;
	double cms_energy = cast_params->beam_energy;
	double var1 = cast_params->var1;
	double mass = cast_params->produced_mass;
	double dummy = 0.0;

	struct ion_params gold179_params;
	gold179_params.atomic_num = 79;
	gold179_params.mass_num = 179;
	gold179_params.energy_CMS = cms_energy;

	struct ion_params pb208_params;
	pb208_params.atomic_num = 82;
	pb208_params.mass_num = 208;
	pb208_params.energy_CMS = cms_energy;

	return (ep_num_total(var1, &gold179_params)/var1) * 
		(ep_num_total(var2, &gold179_params)/var2) *
		fundamental_CS_dilepton(sqrt(4*var1*var2), mass, &dummy);
}

double muon_integrand2_EPA(double var1, void* params)
{
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(500);

	struct dilepton_params* cast_params = (struct dilepton_params*) params;
	double beam_energy = cast_params->beam_energy;
	double mass = cast_params->lepton_mass;
	double result, error;

	struct integrand_params lepton_params;
	lepton_params.beam_energy = beam_energy;
	lepton_params.var1 = var1;

	gsl_function integr_func;
	integr_func.function = &muon_integrand1_EPA;
	integr_func.params = &lepton_params;

	gsl_integration_qag(&integr_func,
			mass*mass / var1, beam_energy*beam_energy / var1,
			1e-10, 1e-3,
			500, 3,
			w, &result, &error);

	gsl_integration_workspace_free(w);

	return result;
}

double dilepton_TCS_EPA(double beam_energy, double lepton_mass, double* err)
{
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(500);

	double result;
	struct dilepton_params params;
	params.beam_energy = beam_energy;
	params.lepton_mass = lepton_mass;

	gsl_function integr_func;
	integr_func.function = &muon_integrand2_EPA;
	integr_func.params = &params;

	gsl_integration_qag(&integr_func,
			lepton_mass, beam_energy,
			1e-10, 1e-3,
			500, 3,
			w, &result, err);

	gsl_integration_workspace_free(w);

	return result;
}


#endif
