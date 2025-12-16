// ***********************************************************************
// Programa de cálculo de número de fótons equivalentes usando a biblioteca
// de cálculo numérico científico GNU Scientific Library (GSL)
// ***********************************************************************

#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>

//Inclusão da bilbioteca das funções para fluxo de fotons
#include "headers/point_like_charge.hpp"
#include "headers/extended_photon_fluxes.hpp"
#include "headers/electron_flux.hpp"

/*
 * Prototipos das funcoes
 */
void outputN_EPA(const char FNAME_N1[50], const char FNAME_N2[50]);
void output_ratio_EPA(const char FNAME[50]);
void output_freq_EPA_ions(const char FNAME[50]);
void output_freq_EPA_energy(const char FNAME[50]);

void output_charge_distribution(const char FNAME[50]);
void output_form_factor(const char FNAME[50]);


int main(int argc, char *argv[])
{
	std::cout << "Iniciando cálculo dos pontos.\n";

//	outputN_EPA("data/N1_map.dat", "data/N2_map.dat");
//	output_ratio_EPA("data/ratio.dat");

	output_freq_EPA_ions("data/n_total.dat");
	output_freq_EPA_energy("data/n_total_energy.dat");

	output_charge_distribution("data/wood_saxon_dist.dat");
	output_form_factor("data/form_factor.dat");

	std::cout << "\nCálculo concluído.\n";

	return 0;
}

// Faz o output dos pontos de N1 e N2 em um mapa bidimensional (x,y), levando
// em conta o requisito do GNUPLOT de colocar uma linha em branco para cada
// atualização do valor de x.
void outputN_EPA(const char FNAME_N1[50], const char FNAME_N2[50])
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

	const double CENTRE_OF_MASS_ENERGY = 500 * GSL_CONST_NUM_GIGA;

	struct ion_params pb208_params;
	pb208_params.atomic_num = 82;
	pb208_params.mass_num = 208;
	pb208_params.energy_CMS = CENTRE_OF_MASS_ENERGY;

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
				<< ep_num_par(freq, b * METRE_TO_EV, &pb208_params) << "\n";

			dados_out2 << freq / GSL_CONST_NUM_MEGA << "\t"
				<< b / GSL_CONST_NUM_FEMTO << "\t"
				<< ep_num_perp(freq, b * METRE_TO_EV, &pb208_params) << "\n";

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
void output_ratio_EPA(const char FNAME[50])
{
	const int NUM_PONTOS = 500;
	const long double LOWER_B = 0.5e-16;
	const long double UPPER_B = 100.0e-15;
	const long double STEP = (UPPER_B - LOWER_B) / NUM_PONTOS;
	const long double FREQ1 = 1.0e+6;
	const long double FREQ2 = 10.0e+6;

	const double CENTRE_OF_MASS_ENERGY = 500 * GSL_CONST_NUM_GIGA;

	struct ion_params pb208_params;
	pb208_params.atomic_num = 82;
	pb208_params.mass_num = 208;
	pb208_params.energy_CMS = CENTRE_OF_MASS_ENERGY;

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
		ratio1 = ep_num_par(FREQ1, b * METRE_TO_EV, &pb208_params)
			/ ep_num_perp(FREQ1, b * METRE_TO_EV, &pb208_params);
		ratio2 = ep_num_par(FREQ2, b * METRE_TO_EV, &pb208_params)
			/ ep_num_perp(FREQ2, b * METRE_TO_EV, &pb208_params);
		dados_out << b / GSL_CONST_NUM_FEMTO << "\t" << ratio1 << "\t" << ratio2 << "\n";
		b += STEP;
	}
	
	std::cout << "\tDados salvos no arquivo \'" << FNAME << "\'\n";

	dados_out.close();
}

// Faz o output do número de fótons totais integrados sobre o parâmetro de
// impacto
void output_freq_EPA_ions(const char FNAME[50])
{
	const int NUM_PONTOS = 6000;

	// Intervalos de frequencia em eletron-Volt
	const double LOWER_F = 1e-6;
	const double UPPER_F = 10e9;

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
	e_params.produced_system_mass = 2*MUON_MASS;

	double freq = LOWER_F;
	double step = fabs(UPPER_F - LOWER_F)/NUM_PONTOS;

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
	std::cout << "\n * Parâmetros da curva n \n";
	std::cout << "\tIntervalo de frequência: [" << LOWER_F << ":" << UPPER_F <<
		"] [eV]\n";
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

	gsl_set_error_handler_off();

	for(int i = 0; i <= NUM_PONTOS; i++){
		dados_out << freq << "\t"
			<< epa_photon_flux(freq, &gold179_params) << "\t"
			<< epa_photon_flux(freq, &pb208_params) << "\t"
			<< electron_photon_flux(freq, &e_params) << "\n";
		freq = freq + step;
	}

	std::cout << "\tDados salvos no arquivo \'" << FNAME << "\'\n";

	dados_out.close();
}

void output_freq_EPA_energy(const char FNAME[50])
{
	const double E1 = 500e9;
	const double E2 = 3e12;
	const double E3 = 5e12;

	struct ion_params pb208_par1;
	pb208_par1.atomic_num = 82;
	pb208_par1.mass_num = 208;
	pb208_par1.energy_CMS = E1;

	struct ion_params pb208_par2;
	pb208_par2.atomic_num = 82;
	pb208_par2.mass_num = 208;
	pb208_par2.energy_CMS = E2;

	struct ion_params pb208_par3;
	pb208_par3.atomic_num = 82;
	pb208_par3.mass_num = 208;
	pb208_par3.energy_CMS = E3;

	const double LOWER_F = 1e-6;
	const double UPPER_F = 10e9;
	const int NUM_PONTOS = 6000;

	double freq = LOWER_F;
	double step = fabs(UPPER_F - LOWER_F)/NUM_PONTOS;

	std::ofstream dados_out(FNAME);
	assert(dados_out.is_open());
	dados_out.setf(std::ios::scientific);
	dados_out.setf(std::ios::showpos);
	dados_out.precision(13);

	for(int i = 0; i <= NUM_PONTOS; i++){
		dados_out << freq << "\t"
			<< epa_photon_flux(freq, &pb208_par1) << "\t"
			<< epa_photon_flux(freq, &pb208_par2) << "\t"
			<< epa_photon_flux(freq, &pb208_par3) << "\n";
		freq = freq + step;
	}


	dados_out.close();
}

void output_charge_distribution(const char FNAME[50])
{
	struct wood_saxon_par pb208_par;
	pb208_par.mass_num = 208;
	pb208_par.atomic_num = 82;

	struct wood_saxon_par au179_par;
	au179_par.mass_num = 179;
	au179_par.atomic_num = 79;

	const double LOWER_RADIUS = 0.0;
	const double UPPER_RADIUS = 10e-15 * METRE_TO_EV; // 100 fm
	const int NUM_PONTOS = 100;

	double radius = LOWER_RADIUS;
	double step = fabs(UPPER_RADIUS - LOWER_RADIUS) / (double) NUM_PONTOS;

	std::ofstream dados_out(FNAME);
	assert(dados_out.is_open());
	dados_out.setf(std::ios::scientific);
	dados_out.setf(std::ios::showpos);
	dados_out.precision(13);

	std::cout << "\n------------------------------------------------------\n";
	std::cout << " * Plotagem da curva de distribuição de wood-saxon:\n"
		<< "\t- LOWER_RADIUS = " << LOWER_RADIUS * 1e15 << " [fm]\n"
		<< "\t- UPPER_RADIUS = " << UPPER_RADIUS * 1e15 / METRE_TO_EV << " [fm]\n"
		<< "\t- NUM_PONTOS = " << NUM_PONTOS << "\n";

	for(int i = 0; i <= NUM_PONTOS; i++) {
		dados_out << radius * 1e15 / METRE_TO_EV << "\t" 
			<< wood_saxon_distribution(radius, &pb208_par) / pow(METRE_TO_EV, 3) << "\t"
			<< wood_saxon_distribution(radius, &au179_par) / pow(METRE_TO_EV, 3) << "\n";
		radius += step;
	}

	std::cout << "\tDados salvos no arquivo \'" << FNAME << "\'\n";

	dados_out.close();
}

void output_form_factor(const char FNAME[50])
{
	struct wood_saxon_par pb208_par;
	pb208_par.mass_num = 208;
	pb208_par.atomic_num = 82;

	struct wood_saxon_par au179_par;
	au179_par.mass_num = 179;
	au179_par.atomic_num = 79;

	const double LOWER_Q = 0.0;
	const double UPPER_Q = 1e9;
	const int NUM_PONTOS = 500;

	double var_q = LOWER_Q;
	double step = fabs(UPPER_Q - LOWER_Q) / (double) NUM_PONTOS;

	std::ofstream dados_out(FNAME);
	assert(dados_out.is_open());
	dados_out.setf(std::ios::scientific);
	dados_out.setf(std::ios::showpos);
	dados_out.precision(13);

	std::cout << "\n------------------------------------------------------\n";
	std::cout << "Plotagem do fator de forma:\n";
	std::cout << "\tLOWER_Q = " << LOWER_Q << "\n";
	std::cout << "\tUPPER_Q = " << UPPER_Q << "\n";
	std::cout << "\tNUM_PONTOS = " << NUM_PONTOS << "\n";

	for(int i = 0; i <= NUM_PONTOS; i++) {
		dados_out << var_q*1e-9 << "\t"
			<< abs(form_factor_wood_saxon(var_q, &pb208_par)) << "\n";
		var_q += step;
	}

	std::cout << "\tDados salvos no arquivo \'" << FNAME << "\'\n";
	dados_out.close();
}
