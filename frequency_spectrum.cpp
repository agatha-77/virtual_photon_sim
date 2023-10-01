// ***********************************************************************
// Programa de cálculo de número de fótons equivalentes usando a biblioteca
// de cálculo numérico científico GNU Scientific Library (GSL)
// ***********************************************************************
//
// O programa faz os cálculos em SI e, após isso, salva o resultado no arquivo
// em unidades naturais. Em nota, achei mais fácil assim.

#include<iostream>
#include<fstream>
#include<cassert>
#include<cmath>

//Inclusão da bilbioteca das funções para partícula pontual.
#include"headers/point_like_charge.hpp"


// Faz o output dos pontos de N1 e N2 em um mapa bidimensional (x,y), levando
// em conta o requisito do GNUPLOT de colocar uma linha em branco para cada
// atualização do valor de x.
void output_N1_and_N2(const char FNAME_N1[50], const char FNAME_N2[50])
{
	const int NUM_PONTOS = 100;

	const double LOWER_B = 0.5e-15;
	const double UPPER_B = 5e-15;

	const double LOWER_F = 0.1e+9;
	const double UPPER_F = 6.0e+9;

	double b = LOWER_B;
	double freq = LOWER_F;
	double step_par = (UPPER_B - LOWER_B)/NUM_PONTOS;
	double step_freq = (UPPER_F - LOWER_F)/NUM_PONTOS;

	std::ofstream dados_out1(FNAME_N1);
	std::ofstream dados_out2(FNAME_N2);
	std::ofstream log_out("data/bess_arg_log.dat");

	assert(dados_out1.is_open());
	assert(dados_out2.is_open());
	assert(log_out.is_open());

	for(int i = 0; i <= NUM_PONTOS; i++){
		b = LOWER_B;

		for(int k = 0; k <= NUM_PONTOS; k++){
			dados_out1 << freq << "\t" << b << "\t" << ep_num_par(freq,b) << "\n";
			dados_out2 << freq << "\t" << b << "\t" << ep_num_perp(freq,b) << "\n";
			log_out << freq << "\t" << b << "\t" << bess_arg(freq, b) << "\n";

			b += step_par;
		}

		dados_out1 << "\n";
		dados_out2 << "\n";
		log_out << "\n";
		freq += step_freq;
	}

	dados_out1.close();
	dados_out2.close();
	log_out.close();
}

void output_ratio(const char FNAME[50])
{
	const int NUM_PONTOS = 500;

	const double LOWER_B = 0.5e-16;
	const double UPPER_B = 50.0e-15;
	const double STEP = (UPPER_B - LOWER_B) / NUM_PONTOS;
	const double FREQ = 10.0e+9;

	std::ofstream dados_out(FNAME);
	assert(dados_out.is_open());

	double b = LOWER_B;
	double ratio;

	for(int i = 0; i <= NUM_PONTOS; i++){
		ratio = ep_num_par(FREQ,b) / ep_num_perp(FREQ,b);
		dados_out << b << "\t" << ratio << "\n";
		b += STEP;
	}

	dados_out.close();
}

void output_n_total(const char FNAME[50])
{
	// Número de pontos do gráfico
	const int NUM_PONTOS = 100;

	// Intervalos de frequencia em eletron-Volt
	const double LOWER_F = 1.0e+9;
	const double UPPER_F = 50.0e+9; 

	double freq = LOWER_F;
	double step = (UPPER_F - LOWER_F)/NUM_PONTOS;

	std::ofstream dados_out(FNAME);
	assert(dados_out.is_open());

	for(int i = 0; i <= NUM_PONTOS; i++){
		dados_out << freq << "\t" << ep_num_total(freq) << "\n";
		freq += step;
	}

	dados_out.close();
}


int main(int argc, char *argv[])
{
	std::cout << "\n";
	std::cout << "-----------------------------------\n";
	std::cout << "-> Iniciando cálculo dos pontos.\n";

	output_N1_and_N2("data/N1_map.dat","data/N2_map.dat");
	output_ratio("data/ratio.dat");
	output_n_total("data/n_total.dat");

	std::cout << "\nCálculo concluído\n\n";
	std::cout << "Dados postos nos arquivos.\n";

	return 0;
}

