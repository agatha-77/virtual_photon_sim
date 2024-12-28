#include <iostream>
#include <fstream>
#include <cmath>

#include <gsl/gsl_errno.h>

#include "headers/cross_section.hpp"
#include "headers/phys_const.hpp"
#include "headers/point_fraction_flux.hpp"

int main()
{
	std::cout.setf(std::ios::scientific);
	std::cout.setf(std::ios::showpos);
	std::cout.precision(8);

	gsl_set_error_handler_off();

	double cms_energy = 1500e9;
	double produced_mass = MUON_MASS;

	struct dilepton_params lepton_params;
	lepton_params.beam_energy = cms_energy;
	lepton_params.lepton_mass = produced_mass;

	struct ion_params pb208_params;
	pb208_params.atomic_num = 82;
	pb208_params.mass_num = 208;
	pb208_params.energy_CMS = cms_energy;

	double var[2];
	var[0] = produced_mass / cms_energy;
	var[1] = produced_mass * produced_mass / (cms_energy*cms_energy * var[0]);

	double dummy = 0.0;

	const int NPONTOS = 5;
	double passo1 = fabs(var[0] - 1.0) / (double) NPONTOS;

	for(int i = 0; i <= NPONTOS; i++){

		var[1] = produced_mass * produced_mass / (cms_energy*cms_energy * var[0]);
		double passo2 = fabs(var[1] - 1.0) / (double) NPONTOS;

		for(int j = 0; j <= NPONTOS; j++){
			std::cout << var[0] << "\t" << var[1] << "\t"
				<< var[0]*cms_energy << "\t"
				<< var[1]*cms_energy << "\t";
			std::cout << ep_num_total(var[0]*cms_energy, &pb208_params) << "\t"
				<< ep_num_total(var[1]*cms_energy, &pb208_params) << "\t"
				<< ep_num_total(var[0]*cms_energy, &pb208_params) / var[0] << "\t"
				<< ep_num_total(var[1]*cms_energy, &pb208_params) / var[1] << "\n";
			var[1] += passo2;
		}
		std::cout << "\n";
		var[0] += passo1;
	}

	struct integrand_params params;
	params.beam_energy = cms_energy;
	params.var1 = var[0];
	params.produced_mass = produced_mass;

	std::cout << "EPA_integrand1 =\t" << muon_integrand1_EPA(var[1], &params) << "\n";
	std::cout << "EPA_integrand2 =\t" << muon_integrand2_EPA(var[0], &lepton_params)
		<< "\n";

	double error;

	std::cout << "VALOR INTEGRAL =\t" << dilepton_TCS_EPA(cms_energy, produced_mass, &error)
		<< "\n";

	return 0;
}
