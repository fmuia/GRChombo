/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(INFLATIONINITIAL_HPP_)
#error "This file should only be included through InflationInitial.hpp"
#endif

#ifndef INFLATIONINITIAL_IMPL_HPP_
#define INFLATIONINITIAL_IMPL_HPP_

// Compute the value of the initial vars on the grid
template <class data_t>
void InflationInitial::compute(Cell<data_t> current_cell) const
{
    Vars<data_t> vars;
    VarsTools::assign(vars, 0.); // Set only the non-zero components below
    Coordinates<data_t> coords(current_cell, m_dx);

    // set the field vars
    vars.phi = m_amplitude_phi;
    vars.Pi = 0;

    // start with unit lapse and flat metric
    vars.lapse = 1;
    vars.chi = 1;
    FOR1(i) vars.h[i][i] = 1.;

    // K needs to be calculated from the energy density
    vars.K = compute_K(vars);

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

// Compute the value of phi at the current point
template <class data_t>
data_t InflationInitial::compute_K(Vars<data_t> vars) const
{
    // set the potential values
    data_t V_of_phi = 0.0;
    data_t dVdphi = 0.0;
    m_potential.compute_potential(V_of_phi, dVdphi, vars);

    data_t H = sqrt(8.0/3.0 * V_of_phi);

    data_t K = -3*H;

    return K;
}

#endif /* INFLATIONINITIAL_IMPL_HPP_ */
