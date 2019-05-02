/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "CCZ4.hpp"
#include "InflationPotential.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // read the problem specific params
        readParams(pp);
    }

    void readParams(GRParmParse &pp)
    {
        // for regridding
        pp.load("regrid_threshold_chi", regrid_threshold_chi);
        pp.load("regrid_threshold_phi", regrid_threshold_phi);

        // Initial and SF data
        pp.load("amplitude_phi", amplitude_phi);
        pp.load("G_Newton", G_Newton, 1.0);
        pp.load("scalar_mass", potential_params.scalar_mass);
    }

    // Problem specific parameters
    Real regrid_threshold_chi, regrid_threshold_phi, amplitude_phi;

    // Initial data for matter and potential
    double G_Newton;
    InflationPotential::params_t potential_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
