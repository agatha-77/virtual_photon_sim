
#include <iostream>
#include <fstream>

#include "headers/cross_section_monte_carlo.hpp"
#include "headers/phys_const.hpp"

int main()
{
	/*
	double err;
	double result = dilepton_TCS_EPA_monte_vegas(5e6, ELECTRON_MASS, &err);
	std::cout << "Valor da seção de choque: " << result << "\n";
	*/
	
	struct dilepton_params params;
	params.beam_energy = 200.0e6;
	params.lepton_mass = ELECTRON_MASS;

	double var[2];
	var[0] = 0.5;
	var[1] = 0.6;

	double result = integrand_EPA_dilepton(var, 2, &params);
	std::cout << "valor do integrando: " << result << "\n";

	return 0;
}
