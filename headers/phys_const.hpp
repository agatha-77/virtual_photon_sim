#ifndef PHYSICAL_CONST
#define PHYSICAL_CONST

#include <gsl/gsl_const_num.h>
#include <gsl/gsl_const_mksa.h>

// Definicao de constantes fisicas e matematicas que serão usadas mais a frente
// As constantes abaixo são definidas em unidades naturais
const double PI = M_PI;	
const double FINE_STRUCT_CONST = GSL_CONST_NUM_FINE_STRUCTURE;

const double LIGHT_VEL = 1.0;
const double PLANCK_REDU = 1.0;

const double ELECTRON_MASS = 0.51099895000 * GSL_CONST_NUM_MEGA / (LIGHT_VEL * LIGHT_VEL);
const double MUON_MASS = 105.6583755 * GSL_CONST_NUM_MEGA / (LIGHT_VEL * LIGHT_VEL);
const double ION_CHARGE = sqrt(4*PI*PLANCK_REDU*LIGHT_VEL*FINE_STRUCT_CONST); 
const double ATOMIC_RADIUS = 1.0423e-5;


// Definição de constantes de multiplicação para passar para unidades de interesse
const double EV_TO_BARN = 389.4 * GSL_CONST_NUM_EXA; // GEV^{-2}
const double METRE_TO_EV = 1.0 / ( 0.197 * GSL_CONST_NUM_MICRO); // Metro pra eV


#endif
