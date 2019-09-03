/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GEOMETRICINITIALCONDITIONS_HPP_
#define GEOMETRICINITIALCONDITIONS_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "GeometricFields.hpp"
#include "GeometricPotential.hpp"
#include "MatterCCZ4.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which creates a bubble of a ComplexScalar field - re and im out of
//! phase
class GeometricInitialConditions
{
  public:
    //! A structure for the input params for scalar field initial conditions
    struct params_t
    {
        double amplitudeSF; //!< Amplitude of bump in initial SF bubble
        std::array<double, CH_SPACEDIM> centerSF; //!< Centre of perturbation
        double widthSF; //!< Width of bump in initial SF bubble
        double r_zero;  //!< Position of bump relative to centre
    };

    //! The constructor for the class
    GeometricInitialConditions(params_t a_params,
                               GeometricPotential a_potential, double a_dx)
        : m_dx(a_dx), m_params(a_params), m_potential(a_potential)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        MatterCCZ4<GeometricFields>::Vars<data_t> vars;
        VarsTools::assign(vars, 0.); // Set all vars to zero
        Coordinates<data_t> coords(current_cell, m_dx, m_params.centerSF);

        // set the field vars - so real and im parts out of phase
        vars.phi1 = 5.0 / sqrt(8.0 * M_PI);
        vars.phi2 = 0.0;
        vars.Pi1 = 0.0;
        vars.Pi2 = 0.0;
        double dphi2da = 5.0; // in reduced Mpl units

        // start with unit lapse and flat metric (must be relaxed for chi)
        vars.lapse = 1.0;
        vars.chi = 1.0;

        // TODO: fudge this for now - should come from potential
        // but will do eventually with IC solver
        // K = -sqrt (24*pi*rho) when no spatial distortions
        data_t f, df_dphi1;
        GeometricFields::calculate_f(vars.phi1, f, df_dphi1);
        // set the potential values
        data_t V_of_phi = 0.0;
        data_t dVdphi1 = 0.0;
        data_t dVdphi2 = 0.0;

        // compute potential
        m_potential.compute_potential(V_of_phi, dVdphi1, dVdphi2, vars);

        // fix K = - 3H and the value of Pi2 in d/da
        data_t H2 =
            (8.0 * M_PI * V_of_phi) / (3.0 - (0.5 * dphi2da * dphi2da * f * f));
        vars.Pi2 = dphi2da / sqrt(8.0 * M_PI) * sqrt(H2);
        data_t rho = V_of_phi + 0.5 * vars.Pi2 * vars.Pi2 * f * f;
        vars.K = -sqrt(24 * M_PI * rho);

        // conformal metric is flat
        FOR1(i) vars.h[i][i] = 1.;

        // Add small peturbations - should resolve for K and chi
        vars.phi1 += m_params.amplitudeSF * (cos(2 * M_PI * coords.x / 32.0) +
                                             cos(2 * M_PI * coords.y / 32.0) +
                                             cos(2 * M_PI * coords.z / 32.0));
        vars.phi2 += m_params.amplitudeSF * (sin(4 * M_PI * coords.x / 32.0) +
                                             sin(4 * M_PI * coords.y / 32.0) +
                                             sin(4 * M_PI * coords.z / 32.0));

        // Store the initial values of the variables
        current_cell.store_vars(vars);
    }

  protected:
    const double m_dx;
    const params_t m_params; //!< The matter initial condition params
    const GeometricPotential m_potential; //!< The matter potential
};

#endif /* GEOMETRICINITIALCONDITIONS_HPP_ */
