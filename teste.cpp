// Arquivo de teste das variaveis
#include<iostream>
#include<cmath>

#include"headers/phys_const.hpp"

double vel_f(double freq){
	// Comprimento de onda de COMPTON
	//double wavelength = 2 * PI * PLANCK_REDU * freq / (ELECTRON_MASS *
		//	LIGHT_VEL);
	// Comprimento de onda de DE BROGLIE
	//double wavelength = LIGHT_VEL / (2 * PI * sqrt(freq*freq +
	//			2*freq*ELECTRON_MASS*LIGHT_VEL*LIGHT_VEL));
	
	// Puramente relativístico
	double vel = 1/(LIGHT_VEL*LIGHT_VEL) - ELECTRON_MASS*ELECTRON_MASS*LIGHT_VEL*LIGHT_VEL*
		(freq*freq + ELECTRON_MASS*ELECTRON_MASS * pow(4,LIGHT_VEL));

	return vel;
}

double bess_arg(double frequency, double imp_par){
	double vel = vel_f(frequency);
	double beta = vel / LIGHT_VEL;
	double gamma = 1.0/sqrt( 1.0 - (beta*beta) );

	return ( frequency * imp_par )/(gamma * vel); 
}

int main(int argc, char *argv[]){
	const double FREQ_MIN = 1.0;
	const double FREQ_MAX = 10.0; 

	const double B_MIN = 1.0;
	const double B_MAX = 10.0;

	const int NUM_PONTOS = 4;
	
	double freq_step = (FREQ_MAX - FREQ_MIN) / NUM_PONTOS;
	double b_step = (B_MAX - B_MIN) / NUM_PONTOS;
	double freq = FREQ_MIN;

	std::cout << "#\n# Resultados com "; 
	std::cout << "[b_step = " << b_step << "] ; [freq_step = " << freq_step << "]\n#\n";
	std::cout << "###\t\tINICIO DOS DADOS\t\t###\n";
	std::cout << "#freq\tb\tvel\t\tbeta\t\tgamma\tbess_arg\n";

	for(int i = 0; i <= NUM_PONTOS; i++){
		double b = B_MIN;
		double beta = vel_f(freq) / LIGHT_VEL;
		double gamma = 1.0/sqrt( 1.0 - (beta*beta) );

		for(int j = 0; j <= NUM_PONTOS; j++){
			std::cout << freq << "\t" << b << "\t" << vel_f(freq) << "\t" << beta
				<< "\t" << gamma << "\t" << bess_arg(freq,b) << "\n";
			b += b_step;
		}

		freq += freq_step;
		std::cout << "\n";
	}

	return 0;
}
