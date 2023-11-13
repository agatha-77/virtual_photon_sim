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
double K0(double b_ARG)
{
	return gsl_sf_bessel_K0(b_ARG);
}

double K1(double b_ARG)
{
	return gsl_sf_bessel_K1(b_ARG);
}

// Função de velocidade em termos da frequencia assumindo comportamento ondulatório
// da partícula. 
double vel_f(double freq){
	// Comprimento de onda de DE BROGLIE
	double wavelength = LIGHT_VEL / (2 * PI * sqrt(freq*freq + 2 * freq *
				ELECTRON_MASS * LIGHT_VEL * LIGHT_VEL));
	double vel = wavelength*freq;
	
	// Velocidade relativística
	/*
	double vel = (1/LIGHT_VEL * LIGHT_VEL) - ELECTRON_MASS * ELECTRON_MASS *
		LIGHT_VEL*LIGHT_VEL * (freq*freq + (ELECTRON_MASS*LIGHT_VEL*LIGHT_VEL) *
				(ELECTRON_MASS*LIGHT_VEL*LIGHT_VEL));
	*/

	return vel;
}

// Argumento da função de Bessel 
double bess_arg(double frequency, double imp_par)
{
	double vel = vel_f(frequency);
	double beta = vel / LIGHT_VEL;
	double gamma = 1.0/sqrt( 1.0 - (beta*beta) );

	return ( frequency * imp_par )/(gamma * vel); 
}


// N(w,b) do pulso paralelo P1
double ep_num_par(double frequency, double imp_par)
{
	double vel = vel_f(frequency);
	double beta = vel / LIGHT_VEL;
	double gamma = 1.0/sqrt( 1.0 - beta*beta );
	double bessel_arg = bess_arg(frequency, imp_par);

	return ION_CHARGE * bessel_arg*bessel_arg * K0(bessel_arg)*K0(bessel_arg) / ( PI*PI
			* beta*beta * frequency * imp_par*imp_par * gamma*gamma * PLANCK_REDU);
}


// N(w,b) do pulso perpendicular P2
double ep_num_perp(double frequency, double imp_par)
{
	double vel = vel_f(frequency);
	double beta = vel / LIGHT_VEL;
	double gamma = 1.0/sqrt( 1.0 - pow(beta,2) );
	double bessel_arg = bess_arg(frequency, imp_par);

	return ION_CHARGE * bessel_arg*bessel_arg * K1(bessel_arg)*K1(bessel_arg) / ( PI*PI
			* beta*beta * frequency * imp_par*imp_par * PLANCK_REDU );
}


// Numero total de fotons equivalentes integrado sobre os parametros de impacto
double ep_num_total(double frequency)
{
	const double IMP_PAR_MIN = GSL_CONST_MKSA_BOHR_RADIUS;

	double vel = vel_f(frequency);
	double beta = vel / LIGHT_VEL;
	double gamma = 1.0/sqrt( 1.0 - pow(beta,2) );
	double bessel_arg = bess_arg(frequency, IMP_PAR_MIN);

	double frontal_mult = 2 * ION_CHARGE / (PI*LIGHT_VEL*pow(beta,2) *
			PLANCK_REDU);

	return frontal_mult * ( bessel_arg * K0(bessel_arg) * K1(bessel_arg) -
			pow(beta,2) * ( pow(K1(bessel_arg),2 ) - pow( K0(bessel_arg),2) ) / 2);
}

#endif
