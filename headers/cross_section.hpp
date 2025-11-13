
/*
 * Biblioteca para calculo de secoes de choque e quantidades relacionadas 
 * (valor maximo de secao de choque e afins).
 *
 * Unidades naturais sao usadas
 */

#ifndef CROSS_SECTION
#define CROSS_SECTION

#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>

#include "phys_const.hpp"
#include "point_like_charge.hpp"
#include "electron_flux.hpp"


// Seções de choque fundamentais dos léptons
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

double fundamental_CS_dilepton(double w, double lepton_mass, void* params)
{
	(void)(params);
	double mass = lepton_mass;

	return (4*PI * FINE_STRUCT_CONST*FINE_STRUCT_CONST / (w*w)) *
		(2* log( (w / (2*mass)) * 
		(1 + sqrt(1 - (4*mass*mass / (w*w))) )) *
		(1 + (4*mass*mass*w*w - 8*mass*mass*mass*mass)/(w*w*w*w) ) -
		sqrt(1 - (4*mass*mass) /(w*w) ) * (1 + (4*mass*mass)/(w*w)) );
}

/* 
 * ****************************************
 * Parte de Integrais
 * ****************************************
 * */



// Estruturas de dados utilizadas nas rotinas
struct integrand_params
{
	double beam_energy; // Energia de feixe
	double var1; // Outra variavel
	double produced_mass; // Massa do sistema produzido
	double (*flux_ptr) (double freq, void* params); // Ponteiro para fluxo
};

struct dilepton_params
{
	double beam_energy; // Energia de feixe
	double lepton_mass; // Massa do sistema produzido
	double (*flux_ptr) (double freq, void* params); // Ponteiro para fluxo
};


// Integrando para EPA
double dilepton_integrand1_EPA(double var2, void* params)
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

	return (cast_params->flux_ptr(var1, &pb208_params) / var1) * 
		(cast_params->flux_ptr(var2, &pb208_params) / var2) *
		fundamental_CS_dilepton(sqrt(4*var1*var2), mass, &dummy);
}

// Integrando "externo"
double dilepton_integrand2_EPA(double var1, void* params)
{
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(5000);

	struct dilepton_params* cast_params = (struct dilepton_params*) params;
	double beam_energy = cast_params->beam_energy;
	double mass = cast_params->lepton_mass;
	double result, error;

	struct ion_params pb208_params;
	pb208_params.atomic_num = 82;
	pb208_params.mass_num = 208;
	pb208_params.energy_CMS = beam_energy;

	const double IMP_PAR_MIN = 5.916 * pow(pb208_params.mass_num,1.0/3.0) *
		(1.0 / GSL_CONST_NUM_GIGA);
	double ion_mass = pb208_params.atomic_num * PROTON_MASS +
		(pb208_params.mass_num - pb208_params.atomic_num) * NEUTRON_MASS;

	double gamma = beam_energy / ion_mass;
	double beta = sqrt(1.0 - 1.0 / (gamma*gamma));

	double max_photon_virt = 4.0*(gamma*beta / IMP_PAR_MIN);

	struct integrand_params lepton_params;
	lepton_params.beam_energy = beam_energy;
	lepton_params.var1 = var1;
	lepton_params.produced_mass = mass;
	lepton_params.flux_ptr = cast_params->flux_ptr;

	gsl_function integr_func;
	integr_func.function = &dilepton_integrand1_EPA;
	integr_func.params = &lepton_params;

	gsl_integration_qags(&integr_func,
			mass*mass / var1, beam_energy,
			1e-10, 1e-6,
			5000,
			w, &result, &error);

	gsl_integration_workspace_free(w);

	return result;
}

// Rotina externa de integracao para secao de choque total
double dilepton_TCS_EPA(double beam_energy, double lepton_mass, 
		double (*flux_ptr) (double freq, void* params), double* err)
{
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(5000);
	gsl_set_error_handler_off();

	struct ion_params pb208_params;
	pb208_params.atomic_num = 82;
	pb208_params.mass_num = 208;
	pb208_params.energy_CMS = beam_energy;

	const double IMP_PAR_MIN = 5.916 * pow(pb208_params.mass_num,1.0/3.0) *
		(1.0 / GSL_CONST_NUM_GIGA);
	double ion_mass = pb208_params.atomic_num * PROTON_MASS +
		(pb208_params.mass_num - pb208_params.atomic_num) * NEUTRON_MASS;

	double gamma = beam_energy / ion_mass;
	double beta = sqrt(1.0 - 1.0 / (gamma*gamma));

	double max_photon_virt = 4.0*(gamma*beta / IMP_PAR_MIN);

	double result;
	struct dilepton_params params;
	params.beam_energy = beam_energy;
	params.lepton_mass = lepton_mass;
	params.flux_ptr = flux_ptr;

	gsl_function integr_func;
	integr_func.function = &dilepton_integrand2_EPA;
	integr_func.params = &params;

	gsl_integration_qags(&integr_func,
			2.0*lepton_mass, beam_energy,
			1e-10, 1e-6,
			5000,
			w, &result, err);

	gsl_integration_workspace_free(w);

	return result;
}

/* 
 * Seções de choque diferenciais à rapidez
 */

// Integrando da função abaixo
double dilepton_DCS_rapidity_integrand_EPA(double w_gamma, void* params)
{
	struct integrand_params* cast_params = (struct integrand_params*) params;
	double beam_energy = cast_params->beam_energy;
	double rapidity = cast_params->var1;
	double lepton_mass = cast_params->produced_mass;
	double dummy = 0.0;

	double omega_1 = w_gamma * exp(rapidity) / 2.0 ;
	double omega_2 = w_gamma * exp(-1.0*rapidity) / 2.0 ;

	struct ion_params pb208_params;
	pb208_params.atomic_num = 82;
	pb208_params.mass_num = 208;
	pb208_params.energy_CMS = beam_energy;
	
	return (cast_params->flux_ptr(omega_1, &pb208_params) / omega_1)
		* (cast_params->flux_ptr(omega_2, &pb208_params) / omega_2)
		* fundamental_CS_dilepton(w_gamma, lepton_mass, &dummy) * (w_gamma / 2.0);
}

// Seção de choque diferencial a rapidez, integral de uma dimensao somente
double dilepton_DCS_rapidity_EPA(double rapidity, void* params, double* err)
{
	int size = 5000;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(size);
	gsl_set_error_handler_off();

	double result;
	struct dilepton_params* cast_params = (struct dilepton_params*) params;
	double lepton_mass = cast_params->lepton_mass;
	double beam_energy = cast_params->beam_energy;

	struct ion_params pb208_params;
	pb208_params.atomic_num = 82;
	pb208_params.mass_num = 208;
	pb208_params.energy_CMS = beam_energy;

	const double IMP_PAR_MIN = 5.916 * pow(pb208_params.mass_num,1.0/3.0) *
		(1.0 / GSL_CONST_NUM_GIGA);
	double ion_mass = pb208_params.atomic_num * PROTON_MASS +
		(pb208_params.mass_num - pb208_params.atomic_num) * NEUTRON_MASS;

	double gamma = beam_energy / ion_mass;
	double beta = sqrt(1.0 - 1.0 / (gamma*gamma));

	double max_photon_virt = 2.0*(gamma*beta / IMP_PAR_MIN);

	struct integrand_params rap_params;
	rap_params.beam_energy = beam_energy;
	rap_params.var1 = rapidity;
	rap_params.produced_mass = 2.0*lepton_mass;
	rap_params.flux_ptr = cast_params->flux_ptr;
	
	gsl_function integr_func;
	integr_func.function = &dilepton_DCS_rapidity_integrand_EPA;
	integr_func.params = &rap_params;

	gsl_integration_qags(&integr_func,
			2*lepton_mass, beam_energy,
			1e-6, 1e-2,
			size, 
			w, &result, err);

	gsl_integration_workspace_free(w);
	
	return result;
}

#endif
