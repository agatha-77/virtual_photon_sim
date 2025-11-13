
/*
 * Biblioteca para o caso de uma partícula com distribuicao pontual de carga.
 *
 * Aqui, estamos usando unidades SI (MKSA) ao inves das naturais, uma vez que
 * sao estas as disponíveis na biblioteca GSL. As funcoes aqui definidas sao
 * feitas de forma direta, uma vez que os calculos para tais e soluvel
 * analiticamente.
 */

#ifndef PHOTON_SPECTRA
#define PHOTON_SPECTRA

#include <cmath>

// Bibliotecas do GSL que são utilizadas aqui
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>

// Constantes que eu defini separadamente
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


// N(w,b) do pulso paralelo P1
double ep_num_par(double frequency, double imp_par, void* params)
{
	struct ion_params* cast_params = (struct ion_params*) params;

	double mass = cast_params->atomic_num * PROTON_MASS + 
		(cast_params->mass_num - cast_params->atomic_num) * NEUTRON_MASS;
	double gamma = cast_params->energy_CMS / mass;

	double beta = sqrt(1.0 - 1.0 / (gamma*gamma));
	double bessel_arg = frequency * imp_par / gamma * beta;

	return FINE_STRUCT_CONST * cast_params->atomic_num * bessel_arg*bessel_arg
		* K0(bessel_arg)*K0(bessel_arg) /
		( PI*PI * beta*beta * frequency * imp_par*imp_par * gamma*gamma);
}


// N(w,b) do pulso perpendicular P2
double ep_num_perp(double frequency, double imp_par, void* params)
{
	struct ion_params* cast_params = (struct ion_params*) params;
	double mass = cast_params->atomic_num * PROTON_MASS + 
		(cast_params->mass_num - cast_params->atomic_num) * NEUTRON_MASS;
	double gamma = (cast_params->energy_CMS + mass) / mass;

	double beta = sqrt(1.0 - 1.0 / (gamma*gamma));
	double bessel_arg = frequency * imp_par / gamma * beta;

	return FINE_STRUCT_CONST * cast_params->atomic_num * bessel_arg*bessel_arg
		* K1(bessel_arg)*K1(bessel_arg) /
		( PI*PI * beta*beta * frequency * imp_par*imp_par * frequency );
}


// Numero total de fotons equivalentes integrado sobre os parametros de impacto
double epa_photon_flux(double frequency, void* params)
{
	struct ion_params* cast_params = (struct ion_params*) params;

	const double IMP_PAR_MIN = 5.916 * pow(cast_params->mass_num,1.0/3.0) *
		(1.0 / GSL_CONST_NUM_GIGA);
	double mass = cast_params->atomic_num * PROTON_MASS +
		(cast_params->mass_num - cast_params->atomic_num) * NEUTRON_MASS;
	double gamma = cast_params->energy_CMS / mass;

	double beta = sqrt(1.0 - 1.0 / (gamma*gamma));
	double bessel_arg = frequency * IMP_PAR_MIN / (gamma * beta);

	double frontal_mult = 2 * FINE_STRUCT_CONST
		* cast_params->atomic_num *cast_params->atomic_num /
		(PI*(beta*beta));

	return frontal_mult * ( bessel_arg * K0(bessel_arg) * K1(bessel_arg) -
			(beta*beta * bessel_arg*bessel_arg) * ( K1(bessel_arg)*K1(bessel_arg) -
			K0(bessel_arg)*K0(bessel_arg) ) / 2.0);
}

#endif
