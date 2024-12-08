/*
 * Programa para imprimir as secoes de choque usando aproximacao dos fotons
 * equivalentes
 */

#include <iostream>
#include <fstream>
#include <cmath>

#include "headers/cross_section.hpp"
#include "headers/phys_const.hpp"
#include "headers/point_fraction_flux.hpp"

void output_to_mass(const char FNAME[50]);
void output_to_energy(const char FNAME[50], double lepton_mass);

int main()
{
	std::cout.setf(std::ios::scientific);
	std::cout.setf(std::ios::showpos);
	std::cout.precision(6);

	output_to_mass("data/TCS_mass_pb208_flux.dat");

	//output_to_energy("data/total_electron_lead_CS.dat", ELECTRON_MASS);
	//output_to_energy("data/total_muon_lead_CS.dat", MUON_MASS);
	//output_to_energy("data/total_tau_lead_CS.dat", TAU_MASS);

	return 0;
}

// Imprime para massa
void output_to_mass(const char FNAME[50])
{
	std::ofstream dados_out(FNAME);
	dados_out.setf(std::ios::scientific);
	dados_out.setf(std::ios::showpos);
	dados_out.precision(13);

	const int NPONTOS1 = 300;
	const double LOWER_M = 1e6;
	const double UPPER_M = 20e9;

	const double ENERGY_BEAM1 = 500e9;
	const double ENERGY_BEAM2 = 1000e9;
	const double ENERGY_BEAM3 = 1500e9;

	const double STEP1 = fabs(UPPER_M - LOWER_M) / (double) NPONTOS1;

	double var_M = LOWER_M;
	double err1, err2, err3;

	std::cout << "-----------------------------------------\n";
	std::cout << "Parâmetros de integracao\n";
	std::cout << "-----------------------------------------\n";

	std::cout << "NPONTOS1 =\t" << NPONTOS1 << "\n";
	std::cout << "cms_energy1 =\t" << ENERGY_BEAM1 << "\ncms_energy2 =\t"
		<< ENERGY_BEAM2 << "\ncms_energy3 =\t" << ENERGY_BEAM3 << "\n";

	for (int i = 0; i < NPONTOS1; i++){
		dados_out << var_M / GSL_CONST_NUM_GIGA << "\t"
			<< dilepton_TCS_EPA(ENERGY_BEAM1, var_M, &err1) * EV_TO_BARN << "\t"
			<< dilepton_TCS_EPA(ENERGY_BEAM2, var_M, &err2) * EV_TO_BARN << "\t"
			<< dilepton_TCS_EPA(ENERGY_BEAM3, var_M, &err3) * EV_TO_BARN << "\t"
			<< err1 << "\t" << err2 << "\t" << err3 << "\n";

		var_M += STEP1;
	}

	dados_out.close();
}

// Imprime para energia
void output_to_energy(const char FNAME[50], double lepton_mass)
{
	const int NPONTOS = 500;
	const double LOWER_E = 100e9;
	const double UPPER_E = 1000e9;

	const double STEP = fabs(UPPER_E - LOWER_E) / (double) NPONTOS;

	std::ofstream dados_out(FNAME);
	dados_out.setf(std::ios::scientific);
	dados_out.setf(std::ios::showpos);
	dados_out.precision(13);

	std::cout.setf(std::ios::scientific);
	std::cout.setf(std::ios::showpos);
	std::cout.precision(13);

	double var_E = LOWER_E;
	double err;
	double result = dilepton_TCS_electron(var_E, MUON_MASS, &err);

	std::cout << "* Calculando secao de choque total para fluxo de chumbo.\n";
	std::cout << "\t" << "var_E\t= " << var_E
		<< "\n\tresult\t= " << result << "\n\terr\t= " << err
	   	<< "\n\tSTEP\t= " << STEP << "\n";

	for (int i = 0; i < NPONTOS; i++){
		dados_out << var_E / GSL_CONST_NUM_GIGA << "\t"
			<< dilepton_TCS_EPA(var_E, lepton_mass, &err) * EV_TO_BARN << "\t"
			<< err << "\n";

		var_E += STEP;
	}

	dados_out.close();

}
