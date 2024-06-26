#include <iostream>
#include <fstream>
#include <cmath>

#include "headers/cross_section.hpp"
#include "headers/phys_const.hpp"

void output_fundamental_electron_CS(const char FNAME[50])
{
	const int NUM_PONTOS = 500;
	const long double LOWER_W = 2*ELECTRON_MASS;
	const long double UPPER_W = 10.0 * GSL_CONST_NUM_GIGA;
	const long double STEP = (UPPER_W - LOWER_W) / NUM_PONTOS;

	std::ofstream dados_out(FNAME);
	
	long double w = LOWER_W;

	for(int i = 0; i < NUM_PONTOS; i++){
		dados_out << w / GSL_CONST_NUM_GIGA << "\t" << fundamental_CS_muon(w) /
			GSL_CONST_NUM_FEMTO << "\n";
		w += STEP;
	}

	dados_out.close();
}


void output_fundamental_muon_CS(const char FNAME[50])
{
	const int NUM_PONTOS = 500;
	const long double LOWER_W = 2*MUON_MASS;
	const long double UPPER_W = 10.0 * GSL_CONST_NUM_GIGA;
	const long double STEP = (UPPER_W - LOWER_W) / NUM_PONTOS;

	std::ofstream dados_out(FNAME);
	
	long double w = LOWER_W;

	for(int i = 0; i < NUM_PONTOS; i++){
		dados_out << w / GSL_CONST_NUM_GIGA << "\t" << fundamental_CS_muon(w) /
			GSL_CONST_NUM_FEMTO << "\n";
		w += STEP;
	}

	dados_out.close();
}

int main()
{
	std::cout << "Realizando cálculo da seção de choque\n";

	output_fundamental_muon_CS("data/fundamental_muon_CS.dat");
	output_fundamental_electron_CS("data/fundamental_electron_CS.dat");

	std::cout << "\n*\tCálculo concluído!\n\n";
	return 0;
}

