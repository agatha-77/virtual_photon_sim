/*
 * *************************************************************************
 * Fluxos de fótons para distribuições de carga extendidas. 
 * *************************************************************************
 *
 * Algumas das integrais são calculadas diretamente. O fluxo total tem de ser
 * feito por integração de Monte Carlo.
 *
 */

#ifndef EXTENDED_PHOTON_FLUXES
#define EXTENDED_PHOTON_FLUXES

#include <cmath>
#include <functional>
#include <thread>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

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
	struct photon_integrand_params* cast_params = (struct photon_integrand_params*) params;

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

struct flux_parameter
{
	int mass_num;
	int atomic_num;
	double beam_energy;
	double (*form_factor_ptr) (double q, void* params);
	double b_0;
	double b_f;
	int npontos;
	int integration_size;
};

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

	int size = 500000;
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

/*
 * Função para calcular os pontos do fluxo de carga extendida.
 * Usa os ponteiros para os arrays e preenche iterativamente.
 * */
void calc_extended_photon_fluxes(double* arg, double *result, double* err,
		double freq, void* params)
{
	struct flux_parameter* cast_params = (struct flux_parameter*) params;

	double mass = cast_params->atomic_num * PROTON_MASS + 
		(cast_params->mass_num - cast_params->atomic_num)*NEUTRON_MASS;
	double gamma = cast_params->beam_energy / mass;
	double beta = sqrt(1.0 - 1.0 / (gamma*gamma));

	struct photon_integrand_params photon_params;
	photon_params.mass_num = cast_params->mass_num;
	photon_params.atomic_num = cast_params->atomic_num;
	photon_params.freq = freq;
	photon_params.ion_mass = mass;
	photon_params.impact_par = cast_params->b_0;
	photon_params.beam_energy = cast_params->beam_energy;
	photon_params.form_factor_ptr = cast_params->form_factor_ptr;

	gsl_integration_workspace* w
		= gsl_integration_workspace_alloc(cast_params->integration_size);
//	gsl_set_error_handler_off();

	gsl_function integr_func;
	integr_func.function = &photon_flux_integrand;
	integr_func.params = &photon_params;

	double var_b = cast_params->b_0;
	double integr_result, error;
	double step = fabs(cast_params->b_f - cast_params->b_0)
		/ (double) cast_params->npontos;

	std::cout << "\n-> Calculo de fluxo com energia E = "
		<< cast_params->beam_energy << "\n";

	for(int i = 0; i <= cast_params->npontos; i++)
	{
		arg[i] = var_b;
		gsl_integration_qagiu(&integr_func, 0.0,
				1e-10, 1e-8,
				cast_params->integration_size, w,
				&integr_result, &error);

		err[i] = error;
		result[i] = ( (cast_params->atomic_num
					* cast_params->atomic_num * FINE_STRUCT_CONST)
			/ (PI*PI * beta*beta * freq * var_b * var_b) ) * fabs(integr_result
			* integr_result);

		std::cout << "var_b=" << arg[i]
			<< "\t flux=" << result[i]
			<< "\t integr_result=" << integr_result 
			<< "\t error=" << err[i] << "\n";
		
		var_b += step;
		photon_params.impact_par = var_b;
	}

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

struct flux_total_parameter
{
	int mass_num;
	int atomic_num;
	double beam_energy;
	double (*form_factor_ptr) (double q, void* params);
	int npontos;
	double freq_0;
	double freq_f;
	int init_call;
	int routine_call;
};

/* 
 * Função para calculo dos pontos.
 *
 * Depende de dois vetores definidos na função principal e que são passados com
 * a dimensão dada dentro do struct.
 * */
void calc_extended_total_photon_flux(double* freq, double* result, double* err,
		void* params)
{
	struct flux_total_parameter* cast_params = (struct flux_total_parameter*) params;
	int mass_num = cast_params->mass_num;
	int atomic_num = cast_params->atomic_num;
	double beam_energy = cast_params->beam_energy;
	int npontos = cast_params->npontos;

	double nuclear_radius = 1.21e-15 * pow(atomic_num, 3) * METRE_TO_EV;

	double mass = cast_params->atomic_num * PROTON_MASS + 
		(cast_params->mass_num - cast_params->atomic_num)*NEUTRON_MASS;
	double gamma = beam_energy / mass;
	double beta = sqrt(1.0 - 1.0 / (gamma*gamma));


	 // chamadas da função de integração, quanto mais alto mais acurado
	const int DIMENSION = 3;
	const double UPPER_BOUND = 1000.0;

	double x0[3] = {0.0, 0.0, 0.0};
	double xf[3] = {UPPER_BOUND, UPPER_BOUND, UPPER_BOUND};

	/*
	 * Expressão lambda para usar multiplas threads no programa.  No momento,
	 * estou usando três, mas auxilia para calcular os pontos mais rapidamente
	 * no meu computador.
	 */
	auto thread_func = [&] (int n_0, int n_final, double freq_0, double freq_f) -> void
	{

		/* struct que é passado ao integrando */
		struct photon_integrand_params photon_params;
		photon_params.mass_num = mass_num;
		photon_params.atomic_num = atomic_num;
		photon_params.ion_mass = mass;
		photon_params.impact_par = 2*nuclear_radius;
		photon_params.beam_energy = beam_energy;
		photon_params.form_factor_ptr = cast_params->form_factor_ptr;

		gsl_rng_env_setup();
		double step = fabs(freq_f - freq_0) / fabs((double) n_final - (double) n_0);

		const gsl_rng_type* T = gsl_rng_default;
		gsl_rng* r = gsl_rng_alloc(T);

		gsl_monte_vegas_state* s = gsl_monte_vegas_alloc(DIMENSION);

		gsl_monte_function F;
		F.f = &total_photon_flux_integrand;
		F.dim = DIMENSION;
		F.params = &photon_params;

		double integr_result;
		double integr_err;
		double var_freq = freq_0;

		for(int i = n_0; i <= n_final; i++)
		{

			freq[i] = var_freq;
			photon_params.freq = var_freq;

			gsl_monte_vegas_init(s);
			gsl_monte_vegas_integrate(&F, x0, xf,
					DIMENSION, cast_params->init_call,
					r, s,
					&integr_result, &integr_err);

			std::cout << "[" << i << "] warm up: " << " var_freq = "
				<< var_freq << "\n"
				<< "[" << i << "] warm up: " << " integr_result = "
				<< integr_result << "\n"
				<< "[" << i << "] warm up: " << " integr_err = "
				<< integr_err << "\n"
				<< "[" << i << "] warm up: " << " chisq = "
				<< gsl_monte_vegas_chisq(s) << "\n";

			int convergence_cycle = 0;
			int integrand_call = cast_params->routine_call;

			do {
				/* Depois de 10 ciclos de integração, aumenta as chamadas de
				 * integração em dez mil a cada ciclo */
				if (convergence_cycle >= 10) integrand_call += 10000;

				std::cout << "[" << i << "] converging: " << " var_freq = " <<
					var_freq << "\n"
					<< "[" << i << "] converging: " << " integr_result = " <<
					integr_result << "\n"
					<< "[" << i << "] converging: " << " integr_err = " <<
					integr_err << "\n"
					<< "[" << i << "] converging: " << " chisq = " <<
					gsl_monte_vegas_chisq(s) << "\n"
					<< "[" << i << "] converging: " << " integration_calls = "
					<< integrand_call << "\n"
					<< "[" << i << "] converging: " << " convergence_cycle = "
					<< convergence_cycle << "\n"; 

				gsl_monte_vegas_integrate(&F, x0, xf,
						DIMENSION, integrand_call,
						r, s,
						&integr_result, &integr_err);

				convergence_cycle += 1;

			} while (fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5);

			result[i] = atomic_num*atomic_num * FINE_STRUCT_CONST
				* fabs(integr_result) / (PI*PI*beta*beta * var_freq);
			err[i] = integr_err;
			var_freq += step;
		}

		gsl_monte_vegas_free(s);
		gsl_rng_free(r);
		
	};

	// Definição dos pontos intermediários, tanto dos vetores quanto do argumento da
	// curva.
	int n_middle1 = npontos / 3;
	int n_middle2 = 2*npontos / 3;
	double freq_middle1 = (cast_params->freq_f - cast_params->freq_0) / 3.0;
	double freq_middle2 = 2.0*(cast_params->freq_f - cast_params->freq_0) / 3.0;

	std::cout << "\nCalculando integral por monte carlo (Vegas)\n";
	std::cout << "\t-> npontos = " << npontos << "\n";
	std::cout << "\t-> n_middle1 = " << n_middle1 << "\n";
	std::cout << "\t-> n_middle2 = " << n_middle2 << "\n";
	std::cout << "\t-> Frequency range = [" << cast_params->freq_0 << ":" << cast_params->freq_f << "]\n";
	std::cout << "\t-> freq_middle1 = " << freq_middle1 << "\n";
	std::cout << "\t-> freq_middle2 = " << freq_middle2 << "\n";
	std::cout << "\t-> dimension = " << DIMENSION << "\n";
	std::cout << "\t-> initial_call = " << cast_params->init_call << "\n";
	std::cout << "\t-> call_routine = " << cast_params->routine_call << "\n\n";

	/* Inicializa as threads para calcular os pontos */
	std::thread thread1(thread_func, 0, n_middle1, cast_params->freq_0, freq_middle1);
	std::thread thread2(thread_func, n_middle1 + 1, n_middle2, freq_middle1, freq_middle2);
	std::thread thread3(thread_func, n_middle2 + 1, npontos, freq_middle2, cast_params->freq_f);

	std::cout << "\t-> thread1 id: " << thread1.get_id() << "\n";
	std::cout << "\t-> thread2 id: " << thread2.get_id() << "\n";
	std::cout << "\t-> thread3 id: " << thread3.get_id() << "\n\n";

	/* Espera as threads terminarem antes de concluir a função */
	thread1.join();
	thread2.join();
	thread3.join();
}

#endif
