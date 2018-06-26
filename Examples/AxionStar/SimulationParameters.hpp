/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"

// Problem specific includes:
#include "CCZ4.hpp"
#include "OscillotonPotential.hpp"

class SimulationParameters
{
  public:
    SimulationParameters(GRParmParse &pp) { readParams(pp); }

    void readParams(GRParmParse &pp)
    {
        // The automatically generated read parameters code defined in
        // SimulationParameters.inc
        auto_read_params(pp);

        // Fill in the potential parameters
        potential_params.scalar_mass = scalar_mass;
        potential_params.f_axion = f_axion;

        // Fill in the ccz4Parameters
        ccz4_params.kappa1 = kappa1;
        ccz4_params.kappa2 = kappa2;
        ccz4_params.kappa3 = kappa3;
        ccz4_params.shift_Gamma_coeff = shift_Gamma_coeff;
        ccz4_params.shift_advec_coeff = shift_advec_coeff;
        ccz4_params.eta = eta;
        ccz4_params.lapse_power = lapse_power;
        ccz4_params.lapse_coeff = lapse_coeff;
        ccz4_params.lapse_advec_coeff = lapse_advec_coeff;

    }

// SimulationParameters.inc declares all variables and defines
// auto_read_params(GRParmParse& pp)
#include "SimulationParameters.inc"

    // Collection of parameters necessary for the CCZ4 RHS
    CCZ4::params_t ccz4_params;
    OscillotonPotential::params_t potential_params;

};

#endif /* SIMULATIONPARAMETERS_HPP_ */
