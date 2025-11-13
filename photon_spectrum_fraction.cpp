#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>

#include "headers/point_fraction_flux.hpp"

void output_flux(const char FNAME[50]);

int main()
{
	std::cout << "Iniciando plotagem da curva.\n";
	gsl_set_error_handler_off();

	output_flux("data/fraction_spectrum.dat");

	std::cout << "\n Plotagem concluída.\n" ;

	return 0;
}

void output_flux(const char FNAME[50])
{
	const int NUM_PONTOS = 5000;

	// Intervalos de frequencia em eletron-Volt
	const double LOWER_FRACTION = 1e-4;
	const double UPPER_FRACTION = 1.0;

	const double CENTER_OF_MASS_ENERGY = 5.02e12;

	struct ion_params gold179_params;
	gold179_params.atomic_num = 79;
	gold179_params.mass_num = 179;
	gold179_params.energy_CMS = CENTER_OF_MASS_ENERGY;

	struct ion_params pb208_params;
	pb208_params.atomic_num = 82;
	pb208_params.mass_num = 208;
	pb208_params.energy_CMS = CENTER_OF_MASS_ENERGY;

	struct electron_params e_params;
	e_params.energy_CMS = CENTER_OF_MASS_ENERGY;
	e_params.produced_system_mass = 5.0e12;

	double freq = LOWER_FRACTION;
	double step = fabs(UPPER_FRACTION - LOWER_FRACTION)/NUM_PONTOS;

	std::ofstream dados_out(FNAME);
	assert(dados_out.is_open());
	dados_out.setf(std::ios::scientific);
	dados_out.setf(std::ios::showpos);
	dados_out.precision(13);

	double pb208_mass = pb208_params.atomic_num * PROTON_MASS +
		(pb208_params.mass_num - pb208_params.atomic_num) * NEUTRON_MASS;
	double gamma_pb = CENTER_OF_MASS_ENERGY / pb208_mass;
	double beta_pb = sqrt(1.0 - 1.0 / (gamma_pb*gamma_pb));
	double gold179_mass = gold179_params.atomic_num * PROTON_MASS +
		(gold179_params.mass_num - gold179_params.atomic_num) * NEUTRON_MASS;
	double gamma_gold = CENTER_OF_MASS_ENERGY / gold179_mass;
	double beta_gold = sqrt(1.0 - 1.0 / (gamma_gold*gamma_gold));

	std::cout << "\n------------------------------------------------------";
	std::cout << "\n * Parâmetros da curva de espectro por fração de energia \n";
	std::cout << "\tIntervalo de fração de energia: [" << LOWER_FRACTION << ":" <<
		UPPER_FRACTION << "]\n";
	std::cout << "\tNúmero de pontos: " << NUM_PONTOS << "\n";
	std::cout << "\tEnergia de centro de massa: " << CENTER_OF_MASS_ENERGY << " [eV]\n\n";
	std::cout << "\tParâmetro de impacto (Pb208): " << 2.0 * 5.196 *
		pow(pb208_params.mass_num, 1.0/3.0) * (1.0/GSL_CONST_NUM_GIGA) << "\n";
	std::cout << "\tmass_Pb208 = " << pb208_mass << "\n"; 
	std::cout << "\tgamma_Pb208 = " << gamma_pb << "\n";
	std::cout << "\tbeta_Pb208 = " << beta_pb << "\n\n";
	std::cout << "\tParâmetro de impacto (Au179): " << 2.0 * 5.196 * 
		pow(gold179_params.mass_num, 1.0/3.0) * (1.0/GSL_CONST_NUM_GIGA) << "\n";
	std::cout << "\tmass_gold179 = " << gold179_mass << "\n";
	std::cout << "\tgamma_gold179 = " << gamma_gold << "\n";
	std::cout << "\tbeta_gold179 = " << beta_gold << "\n\n";

	for(int i = 0; i <= NUM_PONTOS; i++){
		dados_out << freq << "\t"
			<< freq * epa_fraction_flux(freq, &gold179_params) << "\t"
			<< freq * epa_fraction_flux(freq, &pb208_params) << "\t"
			<< freq * electron_fraction_flux(freq, &e_params) << "\n";
		freq = freq + step;

	}

	std::cout << "\tDados salvos no arquivo \'" << FNAME << "\'\n";

	dados_out.close();
}
