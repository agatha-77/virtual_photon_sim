/*
 * ******************************************************************
 * Programa para cálculo das seções de choque fundamental e total dos
 * processos estudados.
 * ******************************************************************
 */

#include <iostream>
#include <fstream>
#include <cmath>

#include "headers/cross_section.hpp"
#include "headers/phys_const.hpp"

void output_fundamental_electron_CS(const char FNAME[50]);
void output_fundamental_muon_CS(const char FNAME[50]);
void output_total_muon_CS(const char FNAME[50]);


int main()
{
	output_total_muon_CS("data/total_muon_CS.dat");

	return 0;
}


void output_fundamental_electron_CS(const char FNAME[50])
{
	const int NUM_PONTOS = 500;
	const long double LOWER_W = 2*ELECTRON_MASS;
	const long double UPPER_W = 10.0 * GSL_CONST_NUM_GIGA;
	const long double STEP = fabs(UPPER_W - LOWER_W) / NUM_PONTOS;

	std::ofstream dados_out(FNAME);
	dados_out.setf(std::ios::scientific);
	dados_out.setf(std::ios::showpos);
	dados_out.precision(13);
	
	long double w = LOWER_W;
	double param = 0.0;

	for(int i = 0; i < NUM_PONTOS; i++){
		dados_out << w / GSL_CONST_NUM_GIGA << "\t" << fundamental_CS_electron(w, &param) *
			EV_TO_BARN / GSL_CONST_NUM_MICRO << "\n";
		w += STEP;
	}

	dados_out.close();
}


void output_fundamental_muon_CS(const char FNAME[50])
{
	const int NUM_PONTOS1 = 10000;
	const int NUM_PONTOS2 = 500;
	const long double LOWER_W = 2.0*MUON_MASS;
	const long double UPPER_W = 10.0 * GSL_CONST_NUM_GIGA;

	const long double MIDDLE_W = max_value_CS(&fundamental_CS_muon, LOWER_W, UPPER_W);
	const long double STEP1 = fabs(MIDDLE_W - LOWER_W) / NUM_PONTOS1;
	const long double STEP2 = fabs(UPPER_W - MIDDLE_W) / NUM_PONTOS2;

	std::ofstream dados_out(FNAME);
	dados_out.setf(std::ios::scientific);
	dados_out.setf(std::ios::showpos);
	dados_out.precision(13);

	std::cout << "* Iniciando plotagem da seção de choque fundamental\n";
	std::cout << "\t- Máximo da seção de choque em torno de: " << MIDDLE_W << " eV\n";
	std::cout << "\t- Número de pontos: " << NUM_PONTOS1 + NUM_PONTOS2 << "\n";
	std::cout << "\t- Intervalo de plotagem: [" << LOWER_W << ":" << UPPER_W << "] [eV]\n\n";
	
	long double w = LOWER_W;
	double param = 0.0;

	for(int i = 0; i < NUM_PONTOS1; i++){
		dados_out << w / GSL_CONST_NUM_GIGA << "\t" << fundamental_CS_muon(w, &param) *
			EV_TO_BARN << "\n";

		w += STEP1;
	}

	w = MIDDLE_W;
	for(int i = 0; i < NUM_PONTOS2; i++){
		dados_out << w / GSL_CONST_NUM_GIGA << "\t" << fundamental_CS_muon(w, &param) *
			EV_TO_BARN << "\n";

		w += STEP2;
	}

	dados_out.close();
}

void output_total_muon_CS(const char FNAME[50])
{
	const int NPONTOS = 10;
	const double LOWER_E = 1.0	* GSL_CONST_NUM_GIGA;
	const double UPPER_E = 10.0 * GSL_CONST_NUM_GIGA;

	const double STEP = fabs(UPPER_E - LOWER_E) / (double) NPONTOS;

	std::ofstream dados_out(FNAME);
	dados_out.setf(std::ios::scientific);
	dados_out.setf(std::ios::showpos);
	dados_out.precision(13);

	double var_E = LOWER_E;
	double err;
	double result = total_muon_cross_section(var_E, &err);

	std::cout << var_E << "\t" << result << "\t" << err << "\n";

	/*
	for (int i = 0; i < NPONTOS; i++){
		dados_out << var_E / GSL_CONST_NUM_GIGA << "\t"
			<< total_muon_cross_section(var_E, &err) << "\t" << err << "\n";

		var_E += STEP;
	}
	*/

	dados_out.close();

}
