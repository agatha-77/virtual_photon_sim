
#include <iostream>
#include <fstream>
#include <cmath>

#include "headers/cross_section.hpp"
#include "headers/phys_const.hpp"

void output_TCS_electron_flux(const char FNAME[50]);

int main()
{
	output_TCS_electron_flux("data/TCS_mass_eflux.dat");

	return 0;
}

void output_TCS_electron_flux(const char FNAME[50])
{
	std::ofstream dados_out(FNAME);
	dados_out.setf(std::ios::scientific);
	dados_out.setf(std::ios::showpos);
	dados_out.precision(13);

	std::cout.setf(std::ios::scientific);
	std::cout.setf(std::ios::showpos);
	std::cout.precision(13);

	const int NPONTOS1 = 300;
	const int NPONTOS2 = 200;
	const double LOWER_M = 1e9;
	const double UPPER_M = 1000e9;

	const double ENERGY_BEAM1 = 500e9;
	const double ENERGY_BEAM2 = 1000e9;
	const double ENERGY_BEAM3 = 1500e9;

	const double STEP1 = fabs(ENERGY_BEAM1 - LOWER_M) / (double) NPONTOS1;
	const double STEP2 = fabs(UPPER_M - ENERGY_BEAM1) / (double) NPONTOS2;

	double var_M = LOWER_M;
	double err;

	for (int i = 0; i < NPONTOS1; i++){
		dados_out << var_M / GSL_CONST_NUM_GIGA << "\t"
			<< dilepton_TCS_electron(ENERGY_BEAM1, var_M, &err) * EV_TO_BARN << "\t"
			<< dilepton_TCS_electron(ENERGY_BEAM2, var_M, &err) * EV_TO_BARN << "\t"
			<< dilepton_TCS_electron(ENERGY_BEAM3, var_M, &err) * EV_TO_BARN << "\n";

		var_M += STEP1;
	}
	
	for (int i = 0; i < NPONTOS2; i++){
		dados_out << var_M / GSL_CONST_NUM_GIGA << "\t"
			<< 0.0 << "\t"
			<< dilepton_TCS_electron(ENERGY_BEAM2, var_M, &err) * EV_TO_BARN << "\t"
			<< dilepton_TCS_electron(ENERGY_BEAM3, var_M, &err) * EV_TO_BARN << "\n";

		var_M += STEP2;
	}


	dados_out.close();
}

