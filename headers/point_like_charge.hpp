#ifndef POINT_CHARGE_DISTRIBUTION
#define POINT_CHARGE_DISTRIBUTION

/*
 * Biblioteca para o caso de uma partícula com distribuição pontual
 * de carga.
 *
 * Aqui, estamos usando unidades SI (MKSA) ao invés das naturais,
 * uma vez que são estas as disponíveis na biblioteca GSL. As funções aqui
 * definidas são feitas de forma direta, uma vez que os cálculos para tais é
 * solúvel analiticamente.
 */

#include<cmath>
#include<fstream>
#include<iostream>

// Bibliotecas do GSL que são utilizadas aqui
#include<gsl/gsl_sf.h>
#include<gsl/gsl_math.h>

// Constantes que eu defini separadamente
#include"phys_const.hpp" 

// Funções de Bessel resumidas para deixar mais legível o resto do codigo
long double K0(long double b_ARG)
{
	return gsl_sf_bessel_K0(b_ARG);
}

long double K1(long double b_ARG)
{
	return gsl_sf_bessel_K1(b_ARG);
}

// Função de velocidade em termos da frequencia assumindo comportamento ondulatório
// da partícula. 
long double vel_f(long double freq){
	// Comprimento de onda de DE BROGLIE
	long double wavelength = LIGHT_VEL / (2 * PI * sqrt(freq*freq + 2 * freq *
				ELECTRON_MASS * LIGHT_VEL * LIGHT_VEL));
	long double vel = wavelength*freq;
	
	// Velocidade relativística
	/*
	long double vel = (1/LIGHT_VEL * LIGHT_VEL) - ELECTRON_MASS * ELECTRON_MASS *
		LIGHT_VEL*LIGHT_VEL * (freq*freq + (ELECTRON_MASS*LIGHT_VEL*LIGHT_VEL) *
				(ELECTRON_MASS*LIGHT_VEL*LIGHT_VEL));
	*/

	return vel;
}

// Argumento da função de Bessel 
long double bess_arg(long double frequency, long double imp_par)
{
	long double vel = vel_f(frequency);
	long double beta = vel / LIGHT_VEL;
	long double gamma = 1.0/sqrt( 1.0 - (beta*beta) );

	return ( frequency * imp_par )/(gamma * vel); 
}


// N(w,b) do pulso paralelo P1
long double ep_num_par(long double frequency, long double imp_par)
{
	long double vel = vel_f(frequency);
	long double beta = vel / LIGHT_VEL;
	long double gamma = 1.0/sqrt( 1.0 - beta*beta );
	long double bessel_arg = bess_arg(frequency, imp_par);

	return ION_CHARGE * bessel_arg*bessel_arg * K0(bessel_arg)*K0(bessel_arg) / ( PI*PI
			* beta*beta * frequency * imp_par*imp_par * gamma*gamma * PLANCK_REDU);
}


// N(w,b) do pulso perpendicular P2
long double ep_num_perp(long double frequency, long double imp_par)
{
	long double vel = vel_f(frequency);
	long double beta = vel / LIGHT_VEL;
	long double gamma = 1.0/sqrt( 1.0 - pow(beta,2) );
	long double bessel_arg = bess_arg(frequency, imp_par);

	return ION_CHARGE * bessel_arg*bessel_arg * K1(bessel_arg)*K1(bessel_arg) / ( PI*PI
			* beta*beta * frequency * imp_par*imp_par * PLANCK_REDU );
}


// Numero total de fotons equivalentes integrado sobre os parametros de impacto
long double ep_num_total(long double frequency)
{
	const long double IMP_PAR_MIN = GSL_CONST_MKSA_BOHR_RADIUS;

	long double vel = vel_f(frequency);
	long double beta = vel / LIGHT_VEL;
	long double gamma = 1.0/sqrt( 1.0 - pow(beta,2) );
	long double bessel_arg = bess_arg(frequency, IMP_PAR_MIN);

	long double frontal_mult = 2 * ION_CHARGE / (PI*LIGHT_VEL*(beta*beta) *
			PLANCK_REDU);

	return frontal_mult * ( bessel_arg * K0(bessel_arg) * K1(bessel_arg) -
			(beta*beta) * ( K1(bessel_arg)*K1(bessel_arg) - K0(bessel_arg)*K0(bessel_arg) ) / 2);
}

#endif
