
/* **************************************************************
 * Programa para teste das variaveis no código
 * **************************************************************
 */

#include <iostream>
#include <cmath>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>

#include "headers/phys_const.hpp"
#include "headers/point_like_charge.hpp"

double fundamental_CS_muon(double arg, void* params)
{
	(void) (params);
	long double w = arg;
	long double mass = MUON_MASS;

	return (4.0*PI * FINE_STRUCT_CONST*FINE_STRUCT_CONST / (w*w)) * (2.0* log(
				(w / (2*mass)) * 
			(1.0 + sqrt(1.0 - ( (4*mass*mass) / (w*w))) )) * (1.0 + (4.0*mass*mass*w*w -
					8.0*mass*mass*mass*mass)/(w*w*w*w) ) 
			- sqrt(1.0 - (4.0*mass*mass) /(w*w) ) * (1.0 + (4.0*mass*mass)/(w*w)) );
}

double muon_integrand(double* var, size_t dim, void* params)
{
	(void) (dim);
	(void) (params);
	double dummy = 0.0;

	return (ep_num_total(var[0]) / var[0]) * (ep_num_total(var[1]) / var[1])
		* fundamental_CS_muon(var[0] + var[1], &dummy);
}

int main(int argc, char *argv[]){
	const double FREQ_MIN = 1.0  * GSL_CONST_NUM_GIGA;
	const double FREQ_MAX = 10.0 * GSL_CONST_NUM_GIGA; 

	const double B_MIN = ATOMIC_RADIUS;
	const double VEL = 0.001 * LIGHT_VEL;

	const int NUM_PONTOS = 5;
	
	double freq_step = (FREQ_MAX - FREQ_MIN) / NUM_PONTOS;
	double freq[2] = {FREQ_MIN, FREQ_MIN};

	std::cout.setf(std::ios::scientific);
	std::cout.setf(std::ios::showpos);
	std::cout.precision(13);

	std::cout << "#\n# Resultados com "; 
	std::cout << "[freq_step = " << freq_step << ", B_MIN = " << B_MIN << "]\n#\n";
	std::cout << "###\t\tINICIO DOS DADOS\t\t###\n";

	double alpha = 0.0;

	for(int i = 0; i <= NUM_PONTOS; i++){
		freq[0] = FREQ_MIN;

		for(int j = 0; j <= NUM_PONTOS; j++){
			std::cout << freq[0] << "\t" << freq[1] << "\t"
				<< ep_num_total(freq[0]) << "\t" << ep_num_total(freq[1]) << "\t"
				<< fundamental_CS_muon(freq[0] + freq[1], &alpha) << "\t"
				<< muon_integrand(freq, 2, &alpha) << "\n";
			freq[0] += freq_step;
		}
		std::cout << "\n";

		freq[1] += freq_step;
	}

	return 0;
}
