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

// Protótipos das subrotinas que fazem a plotagem dos pontos.
// É uma modularização para facilitar a modificação do programa
void output_fundamental_electron_CS(const char FNAME[50]);
void output_fundamental_muon_CS(const char FNAME[50]);
void output_fundamental_tau_CS(const char FNAME[50]);


// Função principal
int main()
{
	output_fundamental_muon_CS("data/fundamental_muon_CS.dat");
	output_fundamental_electron_CS("data/fundamental_electron_CS.dat");
	output_fundamental_tau_CS("data/fundamental_tau_CS.dat");
	//output_total_electron_electron_CS("data/total_electron_electron_CS.dat");
	//output_total_muon_electron_CS("data/total_muon_electron_CS.dat");
	//output_total_tau_electron_CS("data/total_tau_electron_CS.dat");

	return 0;
}


// Imprime a seção de choque fundamental de colisão de fótons em pares de
// elétron-pósitron
void output_fundamental_electron_CS(const char FNAME[50])
{
	const int NUM_PONTOS1 = 10000;
	const int NUM_PONTOS2 = 1000;
	const double LOWER_W = 2*ELECTRON_MASS;
	const double UPPER_W = 10e9;

	const double MIDDLE_W = max_value_CS(&fundamental_CS_electron, LOWER_W, UPPER_W);
	const double STEP1 = fabs(MIDDLE_W - LOWER_W) / NUM_PONTOS1;
	const double STEP2 = fabs(UPPER_W - MIDDLE_W) / NUM_PONTOS2;

	std::ofstream dados_out(FNAME);
	dados_out.setf(std::ios::scientific);
	dados_out.setf(std::ios::showpos);
	dados_out.precision(13);
	
	double w = LOWER_W;
	double param = 0.0;

	for(int i = 0; i < NUM_PONTOS1; i++){
		dados_out << w / GSL_CONST_NUM_GIGA << "\t"
			<< fundamental_CS_electron(w, &param) * EV_TO_BARN << "\n";
		w += STEP1;
	}

	w = MIDDLE_W;
	for(int i = 0; i < NUM_PONTOS2; i++){
		dados_out << w / GSL_CONST_NUM_GIGA << "\t"
			<< fundamental_CS_electron(w, &param) * EV_TO_BARN << "\n";
		w += STEP2;
	}

	dados_out.close();
}


// Imprime a seção de choque fundamental dos pares de muons
void output_fundamental_muon_CS(const char FNAME[50])
{
	const int NUM_PONTOS1 = 10000;
	const int NUM_PONTOS2 = 500;
	const long double LOWER_W = 2.0*MUON_MASS;
	const long double UPPER_W = 10e9;

	// Obtem o maximo da função para separar em dois intervalos de plotagem
	const long double MIDDLE_W = max_value_CS(&fundamental_CS_muon, LOWER_W, UPPER_W);
	const long double STEP1 = fabs(MIDDLE_W - LOWER_W) / NUM_PONTOS1;
	const long double STEP2 = fabs(UPPER_W - MIDDLE_W) / NUM_PONTOS2;

	std::ofstream dados_out(FNAME);
	dados_out.setf(std::ios::scientific);
	dados_out.setf(std::ios::showpos);
	dados_out.precision(13);

	long double w = LOWER_W;
	double param = 0.0;

	for(int i = 0; i < NUM_PONTOS1; i++){
		dados_out << w / GSL_CONST_NUM_GIGA << "\t"
			<< fundamental_CS_muon(w, &param) * EV_TO_BARN << "\n";

		w += STEP1;
	}

	w = MIDDLE_W;
	for(int i = 0; i < NUM_PONTOS2; i++){
		dados_out << w / GSL_CONST_NUM_GIGA << "\t"
			<< fundamental_CS_muon(w, &param) * EV_TO_BARN << "\n";

		w += STEP2;
	}

	dados_out.close();
}


// Imprime a seção de choque do par de tauon
void output_fundamental_tau_CS(const char FNAME[50])
{
	const int NUM_PONTOS1 = 10000;
	const int NUM_PONTOS2 = 500;
	const double LOWER_W = 2*TAU_MASS;
	const double UPPER_W = 10e9;

	const double MIDDLE_W = max_value_CS(&fundamental_CS_tau, LOWER_W, UPPER_W);
	const double STEP1 = fabs(MIDDLE_W - LOWER_W) / NUM_PONTOS1;
	const double STEP2 = fabs(UPPER_W - MIDDLE_W) / NUM_PONTOS2;

	std::ofstream dados_out(FNAME);
	dados_out.setf(std::ios::scientific);
	dados_out.setf(std::ios::showpos);
	dados_out.precision(13);
	
	double w = LOWER_W;
	double param = 0.0;

	for(int i = 0; i < NUM_PONTOS1; i++){
		dados_out << w / GSL_CONST_NUM_GIGA << "\t"
			<< fundamental_CS_tau(w, &param) * EV_TO_BARN << "\n";
		w += STEP1;
	}

	w = MIDDLE_W;
	for(int i = 0; i < NUM_PONTOS2; i++){
		dados_out << w / GSL_CONST_NUM_GIGA << "\t"
			<< fundamental_CS_tau(w, &param) * EV_TO_BARN << "\n";
		w += STEP2;
	}

	dados_out.close();
}


