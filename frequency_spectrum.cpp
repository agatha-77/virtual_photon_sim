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

//Inclusão da bilbioteca das funções para fluxo de fotons
#include"headers/point_like_charge.hpp"
#include"headers/electron_flux.hpp"


// Faz o output dos pontos de N1 e N2 em um mapa bidimensional (x,y), levando
// em conta o requisito do GNUPLOT de colocar uma linha em branco para cada
// atualização do valor de x.
void output_N1_and_N2(const char FNAME_N1[50], const char FNAME_N2[50])
{
	const int NUM_PONTOS = 50;

	const long double LOWER_B = 0.5e-16;
	const long double UPPER_B = 50e-15;

	const long double LOWER_F = 0.1e+6;
	const long double UPPER_F = 50.0e+6;

	long double b = LOWER_B;
	long double freq = LOWER_F;
	long double step_par = (UPPER_B - LOWER_B)/NUM_PONTOS;
	long double step_freq = (UPPER_F - LOWER_F)/NUM_PONTOS;

	const double CENTRE_OF_MASS_ENERGY = 200 * GSL_CONST_NUM_GIGA;

	struct ion_params gold179_params;
	gold179_params.atomic_num = 79;
	gold179_params.mass_num = 179;
	gold179_params.energy_CMS = CENTRE_OF_MASS_ENERGY;

	std::ofstream dados_out1(FNAME_N1);
	std::ofstream dados_out2(FNAME_N2);

	assert(dados_out1.is_open());
	dados_out1.setf(std::ios::scientific);
	dados_out1.setf(std::ios::showpos);
	dados_out1.precision(13);

	assert(dados_out2.is_open());
	dados_out2.setf(std::ios::scientific);
	dados_out2.setf(std::ios::showpos);
	dados_out2.precision(13);

	std::cout << "\n------------------------------------------------------";
	std::cout << "\n * Parâmetros dos mapas bidimensionais de N1 e N2 \n";
	std::cout << "\tIntervalo de parâmetro de impacto:";
	std::cout << LOWER_B / GSL_CONST_NUM_FEMTO << " fm e "
		<< UPPER_B / GSL_CONST_NUM_FEMTO << " fm\n"; 
	std::cout << "\tIntervalo de frequência: " << LOWER_F / GSL_CONST_NUM_MEGA << " GeV e "
		<< UPPER_F / GSL_CONST_NUM_MEGA << " GeV\n";
	std::cout << "\tNúmero de pontos = " << NUM_PONTOS*NUM_PONTOS << "\n";

	for(int i = 0; i <= NUM_PONTOS; i++){
		b = LOWER_B;

		for(int k = 0; k <= NUM_PONTOS; k++){
			dados_out1 << freq / GSL_CONST_NUM_MEGA << "\t"
				<< b / GSL_CONST_NUM_FEMTO << "\t"
				<< ep_num_par(freq, b * METRE_TO_EV, &gold179_params) << "\n";

			dados_out2 << freq / GSL_CONST_NUM_MEGA << "\t"
				<< b / GSL_CONST_NUM_FEMTO << "\t"
				<< ep_num_perp(freq, b * METRE_TO_EV, &gold179_params) << "\n";

			b += step_par;
		}

		dados_out1 << "\n";
		dados_out2 << "\n";
		freq += step_freq;
	}

	std::cout << "\tDados salvos nos arquivos";
	std::cout << " \'" << FNAME_N1 << "\' e \'" << FNAME_N2 << "\'\n";

	dados_out1.close();
	dados_out2.close();
}

// Printa a razão entre o N1 e o N2 para uma escolha da frequencia
void output_ratio(const char FNAME[50])
{
	const int NUM_PONTOS = 500;
	const long double LOWER_B = 0.5e-16;
	const long double UPPER_B = 100.0e-15;
	const long double STEP = (UPPER_B - LOWER_B) / NUM_PONTOS;
	const long double FREQ1 = 1.0e+6;
	const long double FREQ2 = 10.0e+6;

	const double CENTRE_OF_MASS_ENERGY = 200 * GSL_CONST_NUM_GIGA;

	struct ion_params gold179_params;
	gold179_params.atomic_num = 79;
	gold179_params.mass_num = 179;
	gold179_params.energy_CMS = CENTRE_OF_MASS_ENERGY;

	std::ofstream dados_out(FNAME);
	assert(dados_out.is_open());
	dados_out.setf(std::ios::scientific);
	dados_out.setf(std::ios::showpos);
	dados_out.precision(13);

	long double b = LOWER_B;
	long double ratio1, ratio2;

	std::cout << "\n------------------------------------------------------";
	std::cout << "\n * Parâmetros da curva N1/N2 \n";
	std::cout << "\tIntervalo de parâmetro de impacto: ";
	std::cout << LOWER_B / GSL_CONST_NUM_FEMTO << " fm e " << UPPER_B /
		GSL_CONST_NUM_FEMTO << " fm \n";
	std::cout << "\tFrequências utilizadas: " << FREQ1 / GSL_CONST_NUM_GIGA << " GeV e ";
	std::cout << FREQ2 / GSL_CONST_NUM_GIGA << " GeV\n";
	std::cout << "\tNúmero de pontos = " << NUM_PONTOS << "\n";

	for(int i = 0; i <= NUM_PONTOS; i++){
		ratio1 = ep_num_par(FREQ1, b * METRE_TO_EV, &gold179_params)
			/ ep_num_perp(FREQ1, b * METRE_TO_EV, &gold179_params);
		ratio2 = ep_num_par(FREQ2, b * METRE_TO_EV, &gold179_params)
			/ ep_num_perp(FREQ2, b * METRE_TO_EV, &gold179_params);
		dados_out << b / GSL_CONST_NUM_FEMTO << "\t" << ratio1 << "\t" << ratio2 << "\n";
		b += STEP;
	}
	
	std::cout << "\tDados salvos no arquivo \'" << FNAME << "\'\n";

	dados_out.close();
}

// Faz o output do número de fótons totais integrados sobre o parâmetro de
// impacto
void output_n_total(const char FNAME[50])
{
	const int NUM_PONTOS = 2000;

	// Intervalos de frequencia em eletron-Volt
	const double LOWER_F = 1e-6;
	const double UPPER_F = 1e9;

	const double CENTER_OF_MASS_ENERGY = 200e9;

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
	e_params.produced_system_mass = 2*MUON_MASS;

	double freq = LOWER_F;
	double step = fabs(UPPER_F - LOWER_F)/NUM_PONTOS;

	std::ofstream dados_out(FNAME);
	assert(dados_out.is_open());
	dados_out.setf(std::ios::scientific);
	dados_out.setf(std::ios::showpos);
	dados_out.precision(13);
	
	std::cout << "\n------------------------------------------------------";
	std::cout << "\n * Parâmetros da curva n \n";
	std::cout << "\tIntervalo de frequência: [" << LOWER_F << ":" << UPPER_F << "] [eV]\n";
	std::cout << "\tNúmero de pontos: " << NUM_PONTOS << "\n";
	std::cout << "\tEnergia de centro de massa: " << CENTER_OF_MASS_ENERGY << " [eV]\n";

	for(int i = 0; i <= NUM_PONTOS; i++){
		dados_out << freq << "\t"
			<< ep_num_total(freq, &gold179_params) << "\t"
			<< ep_num_total(freq, &pb208_params) << "\t"
			<< electron_photon_flux(freq, &e_params) << "\n";
		freq = freq + step;
	}

	std::cout << "\tDados salvos no arquivo \'" << FNAME << "\'\n";

	dados_out.close();
}


int main(int argc, char *argv[])
{
	std::cout << "Iniciando cálculo dos pontos.\n";

	output_N1_and_N2("data/N1_map.dat","data/N2_map.dat");
	output_ratio("data/ratio.dat");
	output_n_total("data/n_total.dat");

	std::cout << "\nCálculo concluído\n\n";

	return 0;
}

