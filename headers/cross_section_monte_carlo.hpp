/*
 * *************************************************************************
 * Biblioteca para calculo das secoes de choque usando integração de monte
 * carlo.
 * *************************************************************************
 * 
 */

#ifndef CROSS_SECTION_MONTE_CARLO
#define CROSS_SECTION_MONTE_CARLO

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>

#include "phys_const.hpp"
#include "point_like_charge.hpp"
#include "extended_photon_fluxes.hpp"

struct cross_integrand_params
{
	double ion_energy;
	double xi;
	double (*form_factor_ptr) (double q, void* params); /* Vai mais coisas aqui acho. */
};


double cross_section_integrand(double* arg, void* params)
{
	// ... TODO
	// escrever o integrando com bases nas notas
}

/* Struct para a função abaixo */
struct cross_params
{
	double ion_energy;
	double lepton_produced_mass;
	double (*form_factor_ptr) (double q, void* params); /* Ponteiro para o fator de forma */
	int atomic_num;
	int mass_num;
	int initial_calls;
	int routine_calls;
	int npontos;
};


/* 
 * Função de integração para a integral de seção de choque total com densidade
 * de carga do íon incidente extendida
 */
void calc_cross_section_to_energy(double* energy, double* result, void* params)
{
	struct cross_params* cast_params = (struct cross_params*) params;
	double ion_mass = cast_params->atomic_num * PROTON_MASS + 
		(cast_params->mass_num - cast_params->atomic_num)*NEUTRON_MASS;

	double gamma = cast_params->ion_energy / cast_params->ion_mass;
	double beta = sqrt(1.0 - 1.0 / (gamma*gamma));
	double nuclear_radius = 1.21e-15 * pow(atomic_num, 3) * METRE_TO_EV;

	const int DIMENSION = 8;
	const double UPPER_BOUND = 800.0;

	double x0[DIMENSION] = {
		cast_params->lepton_produced_mass, 0.0,
		2*nuclear_radius, 2*nuclear_radius,
		0.0, 0.0,
		0.0, 0.0
	};

	double xf[DIMENSION] = {
		cast_params->ion_energy, cast_params->ion_energy,
		UPPER_BOUND, UPPER_BOUND,
		UPPER_BOUND, UPPER_BOUND,
		UPPER_BOUND, UPPER_BOUND
	};


}

#endif
