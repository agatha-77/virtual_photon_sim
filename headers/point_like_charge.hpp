
/*
 * Biblioteca para o caso de uma partícula com distribuição pontual de carga.
 *
 * Aqui, estamos usando unidades SI (MKSA) ao invés das naturais, uma vez que
 * são estas as disponíveis na biblioteca GSL. As funções aqui definidas são
 * feitas de forma direta, uma vez que os cálculos para tais é solúvel
 * analiticamente.
 */

#ifndef POINT_CHARGE_DISTRIBUTION
#define POINT_CHARGE_DISTRIBUTION

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

// Função de velocidade em termos da frequencia assumindo comportamento
// ondulatório da partícula. 
double vel_f(double freq, double mass){

	double vel = sqrt(1 - (mass*mass) / (freq*freq)) * LIGHT_VEL;

	return vel;
}

// Argumento da função de Bessel 
double bess_arg(double frequency, double imp_par)
{
	double vel = vel_f(frequency, ELECTRON_MASS);
	double beta = vel / LIGHT_VEL;
	double gamma = 1.0/sqrt( 1.0 - (beta*beta) );

	return ( frequency * imp_par )/(gamma * vel); 
}


// N(w,b) do pulso paralelo P1
double ep_num_par(double frequency, double imp_par) {
	double vel = vel_f(frequency, ELECTRON_MASS);
	double beta = vel / LIGHT_VEL;
	double gamma = 1.0/sqrt( 1.0 - beta*beta );
	double bessel_arg = bess_arg(frequency, imp_par);

	return ION_CHARGE * bessel_arg*bessel_arg * K0(bessel_arg)*K0(bessel_arg) /
		( PI*PI * beta*beta * frequency * imp_par*imp_par * gamma*gamma *
		  PLANCK_REDU);
}


// N(w,b) do pulso perpendicular P2
double ep_num_perp(double frequency, double imp_par)
{
	double vel = vel_f(frequency, ELECTRON_MASS);
	double beta = vel / LIGHT_VEL;
	double bessel_arg = bess_arg(frequency, imp_par);

	return ION_CHARGE * bessel_arg*bessel_arg * K1(bessel_arg)*K1(bessel_arg) /
		( PI*PI * beta*beta * frequency * imp_par*imp_par * PLANCK_REDU );
}


// Numero total de fotons equivalentes integrado sobre os parametros de impacto
double ep_num_total(double frequency)
{
	const double IMP_PAR_MIN = 2.0 * GSL_CONST_MKSA_BOHR_RADIUS * METRE_TO_EV;
//	const double IMP_PAR_MIN = ATOMIC_RADIUS;

	double vel = vel_f(frequency, ELECTRON_MASS);
	double beta = vel / LIGHT_VEL;
	double bessel_arg = bess_arg(frequency, IMP_PAR_MIN);

	double frontal_mult = 2 * ION_CHARGE * ION_CHARGE / (PI*(beta*beta) * frequency);

	return frontal_mult * ( bessel_arg * K0(bessel_arg) * K1(bessel_arg) -
			(beta*beta * bessel_arg*bessel_arg) * ( K1(bessel_arg)*K1(bessel_arg) -
				K0(bessel_arg)*K0(bessel_arg) ) / 2.0);
}

#endif
