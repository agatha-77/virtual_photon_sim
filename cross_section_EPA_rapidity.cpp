/*
 * Programa para obter a secao de choque diferencial a rapidez
 */

#include <iostream>
#include <fstream>
#include <cmath>

#include "headers/cross_section.hpp"
#include "headers/phys_const.hpp"
#include "headers/point_like_charge.hpp"

void output_to_rapidity(const char FNAME[50], double lepton_mass);
void test_rapidity(const char FNAME[50]);

int main()
{
	output_to_rapidity("data/DCSR_muon_lead.dat", MUON_MASS);
	output_to_rapidity("data/DCSR_electron_lead.dat", ELECTRON_MASS);
	output_to_rapidity("data/DCSR_tau_lead.dat", TAU_MASS);

	//test_rapidity("data/test.data");

	return 0;
}

void output_to_rapidity(const char FNAME[50], double lepton_mass)
{
	// Numero de pontos, intervalo de rapidez e step
	const int NPONTOS = 150;
	const double LOWER_Y = -1.0;
	const double UPPER_Y = 1.0;

	const double STEP = fabs(UPPER_Y - LOWER_Y) / (double) NPONTOS;

	// Abre o arquivo e coloca a precisao numerica
	std::ofstream dados_out(FNAME);
	dados_out.setf(std::ios::scientific);
	dados_out.setf(std::ios::showpos);
	dados_out.precision(13);

	double var_Y = LOWER_Y;
	double err1, err2, err3;

	const double E1 = 1e12;
	const double E2 = 3e12;
	const double E3 = 5e12;

	// Parametros da colisao (energia do feixe, massa do lepton e fluxo escolhido)
	struct dilepton_params params1 = {E1, lepton_mass, &epa_photon_flux};
	struct dilepton_params params2 = {E2, lepton_mass, &epa_photon_flux};
	struct dilepton_params params3 = {E3, lepton_mass, &epa_photon_flux};

	const double PB_MASS = 82.0 * PROTON_MASS + (208.0 - 82.0) * NEUTRON_MASS;

	std::cout
		<< "Calculando secao de choque diferencial a rapidez para fluxo de chumbo.\n";
	std::cout << "-> massa_lepton = " << lepton_mass << "\n";
	std::cout << "-> massa chumbo = " << PB_MASS << "\n";
	std::cout << "-> gamma1 = " << E1 / PB_MASS << "\n";
	std::cout << "-> beta1 = " << sqrt(1 - PB_MASS*PB_MASS/(E1*E1)) << "\n";
	std::cout << "-> gamma2 = " << E2 / PB_MASS << "\n";
	std::cout << "-> beta2 = " << sqrt(1 - PB_MASS*PB_MASS/(E2*E2)) << "\n";
	std::cout << "-> gamma3 = " << E3 / PB_MASS << "\n";
	std::cout << "-> beta3 = " << sqrt(1 - PB_MASS*PB_MASS/(E3*E3)) << "\n\n";

	for(int i = 0; i < NPONTOS; i++){
		dados_out << var_Y
			<< "\t" << dilepton_DCS_rapidity_EPA(var_Y, &params1,
					&err1)*1e12*EV_TO_BARN
			<< "\t" << err1
			<< "\t" << dilepton_DCS_rapidity_EPA(var_Y, &params2,
					&err2)*1e12*EV_TO_BARN
			<< "\t" << err2
			<< "\t" << dilepton_DCS_rapidity_EPA(var_Y, &params3,
					&err3)*1e12*EV_TO_BARN
			<< "\t" << err2
			<< "\n";

		var_Y += STEP;
	}

	dados_out.close();
}

void test_rapidity(const char FNAME[50])
{
	const int NPONTOS = 70;
	const double LOWER_Y = -2.0;
	const double UPPER_Y = 2.0;
	const double E = 1e12;

	struct ion_params pb208_params;
	pb208_params.atomic_num = 82;
	pb208_params.mass_num = 208;
	pb208_params.energy_CMS = E;

	std::cout << "-------------------------------------\n\n";

	const double PB_MASS = pb208_params.atomic_num * PROTON_MASS + 
		(pb208_params.mass_num - pb208_params.atomic_num) * NEUTRON_MASS;
	const double IMP_PAR_MIN = 5.916 * pow(pb208_params.mass_num,1.0/3.0) *
		(1.0 / GSL_CONST_NUM_GIGA);
	double gamma = E / PB_MASS;
	double beta = sqrt(1.0 - 1.0 / (gamma*gamma));

	std::cout << "\tgamma = " << gamma << "\n";
	std::cout << "\tbeta = " << beta << "\n";
	std::cout << "\tIMP_PAR_MIN = " << IMP_PAR_MIN << "\n\n";

	const double LOWER_W = 2.0*MUON_MASS;
	const double UPPER_W = 4.0*gamma*beta / IMP_PAR_MIN;

	//const double LOWER_W = 0.0;
	//const double UPPER_W = E;

	std::cout << "\tLOWER_W = " << LOWER_W << "\n";
	std::cout << "\tUPPER_W = " << UPPER_W << "\n";

	const double STEP1 = fabs(UPPER_Y - LOWER_Y) / (double) NPONTOS;
	const double STEP2 = fabs(UPPER_W - LOWER_W) / (double) NPONTOS;

	std::cout << "\n\tSTEP1 = " << STEP1 << "\n";
	std::cout << "\tSTEP2 = " << STEP2 << "\n";

	double dummy;

	double var_Y = LOWER_Y;
	double var_W = LOWER_W;

	std::ofstream dados_out(FNAME);
	dados_out.setf(std::ios::scientific);
	dados_out.setf(std::ios::showpos);
	dados_out.precision(13);

	double omega_1 = (var_W / 2.0) * exp(var_Y);
	double omega_2 = (var_W / 2.0) * exp((-1.0)*var_Y);

	struct integrand_params rap_params;
	rap_params.beam_energy = E;
	rap_params.var1 = var_Y;
	rap_params.produced_mass = MUON_MASS;

	gsl_set_error_handler_off();

	for(int i = 0; i <= NPONTOS; i++){
		var_W = LOWER_W;

		for(int j = 0; j <= NPONTOS; j++){
			omega_1 = (var_W / 2.0) * exp(var_Y);
			omega_2 = (var_W / 2.0) * exp((-1.0)*var_Y);

			dados_out << var_Y << "\t" << var_W << "\t"
				<< omega_1 << "\t" << omega_2 << "\t"
				<< epa_photon_flux(omega_1, &pb208_params) / omega_1 << "\t"
				<< epa_photon_flux(omega_2, &pb208_params) / omega_2 << "\t"
				<< fundamental_CS_dilepton(var_W, MUON_MASS, &dummy) * (var_W / 2.0)
				<< "\t"
				<< dilepton_DCS_rapidity_integrand_EPA(var_W, &rap_params) << "\n";
			var_W += STEP2;
		}

		var_Y += STEP1;
		rap_params.var1 = var_Y;
		dados_out << "\n";
	}

	dados_out.close();

}
