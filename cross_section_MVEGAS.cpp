
#include <iostream>
#include <fstream>

#include "headers/cross_section_monte_carlo.hpp"
#include "headers/phys_const.hpp"

int main()
{
	struct dilepton_params params;
	params.beam_energy = 500e9;
	params.lepton_mass = MUON_MASS;
	double cms_energy = params.beam_energy;

	struct ion_params pb208_params;
	pb208_params.atomic_num = 82;
	pb208_params.mass_num = 208;
	pb208_params.energy_CMS = cms_energy;

	double var[2];
	var[0] = params.lepton_mass / cms_energy;
	var[1] = params.lepton_mass * params.lepton_mass / (cms_energy*cms_energy * var[0]);

	double dummy = 0.0;

	std::cout << "--------------------------------------\n";
	std::cout << "Parâmetros iniciais do integrando \n";
	std::cout << "--------------------------------------\n";
	std::cout << "var[0] =\t" << var[0] << "\tomega_1 =\t" << var[0]*cms_energy << "\n";
	std::cout << "var[1] =\t" << var[1] << "\tomega_2 =\t" << var[1]*cms_energy << "\n";

	std::cout << "epa(omega_1) =\t"
		<< ep_num_total(var[0]*cms_energy, &pb208_params) << "\t";
	std::cout << "epa(omega_2) =\t"
		<< ep_num_total(var[1]*cms_energy, &pb208_params) << "\n";
	std::cout << "f(omega_1) =\t"
		<< ep_num_total(var[0]*cms_energy, &pb208_params) / var[0] << "\t";
	std::cout << "f(omega_2) =\t"
		<< ep_num_total(var[1]*cms_energy, &pb208_params) / var[1] << "\n";

	std::cout << "fundamental_CS =\t"
		<< fundamental_CS_dilepton(sqrt(4*var[0]*var[1]*cms_energy*cms_energy),
				params.lepton_mass, &dummy) << "\n\n";

	double result = integrand_EPA_dilepton(var, 2, &params);
	std::cout << "valor do integrando:\t" << result << "\n\n";

	double err;
	double result2 = dilepton_TCS_EPA_monte_vegas(cms_energy, MUON_MASS, &err);
	std::cout << "Valor da seção de choque: " << result2 << "\n";

	return 0;
}
