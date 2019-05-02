/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef OSCILLOTONTAGGINGCRITERION_HPP_
#define OSCILLOTONTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "ScalarField.hpp"
#include "SimulationParametersBase.hpp"
#include "Tensor.hpp"

class OscillotonTaggingCriterion
{
  protected:
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    const double m_threshold_chi;
    const double m_threshold_phi;
    const double m_cutoff_phi;
    const int m_level;
    const extraction_params_t m_params;

    template <class data_t>
    using MatterVars = typename ScalarField<>::template Vars<data_t>;

    /// Vars object for chi
    template <class data_t> struct Vars
    {
        data_t chi; //!< Conformal factor

        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            using namespace VarsTools; // define_enum_mapping is part of
                                       // VarsTools
            define_enum_mapping(mapping_function, c_chi, chi);
        }
    };

  public:
    OscillotonTaggingCriterion(
        const double dx, const double threshold_chi, const double threshold_phi,
        const double cutoff_phi,
        const int a_level, const extraction_params_t a_params)
        : m_dx(dx), m_deriv(dx), m_threshold_chi(threshold_chi),
          m_cutoff_phi(cutoff_phi), m_threshold_phi(threshold_phi), m_params(a_params),
          m_level(a_level){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto d2 = m_deriv.template diff2<MatterVars>(current_cell);
        const auto d2chi = m_deriv.template diff2<Vars>(current_cell);

        data_t mod_d2_chi = 0;
        data_t mod_d2_phi = 0;

        FOR2(idir, jdir)
        {
           mod_d2_chi += d2chi.chi[idir][jdir] * d2chi.chi[idir][jdir]; 

	   mod_d2_phi += d2.Pi[idir][jdir] * d2.Pi[idir][jdir]
                      +  d2.phi[idir][jdir] * d2.phi[idir][jdir];
        }

        data_t criterion_chi = m_dx / m_threshold_chi * sqrt(mod_d2_chi);

        data_t criterion_phi = m_dx / m_threshold_phi * sqrt(mod_d2_phi);

        data_t criterion = simd_max(criterion_chi, criterion_phi);

        // regrid if within extraction level and not at required refinement
        double L = 64.0;
        if (m_level == 0)
        {
            const Coordinates<data_t> coords(current_cell, m_dx);
            if (coords.x < 0.5*L && coords.y < 0.5*L && coords.z < 0.5*L)
            {
                criterion = 100;
            }
        }
        if (m_level == 1)
        {
            const Coordinates<data_t> coords(current_cell, m_dx);
            if (coords.x < 0.25*L && coords.y < 0.25*L && coords.z < 0.25*L)
            {
                criterion = 100;
            }
        }
        if (m_level == 2)
        {
            const Coordinates<data_t> coords(current_cell, m_dx);
            if (coords.x < 0.125*L && coords.y < 0.125*L && coords.z < 0.125*L)
            {
                criterion = 100;
            }
        }
        if (m_level == 3)
        {
            const Coordinates<data_t> coords(current_cell, m_dx);
            if (coords.x < 0.0625*L && coords.y < 0.0625*L && coords.z < 0.0625*L)
            {
                criterion = 100;
            }
        }
        if (m_level == 4)
        {
            const Coordinates<data_t> coords(current_cell, m_dx);
            if (coords.x < 0.03125*L && coords.y < 0.03125*L && coords.z < 0.03125*L)
            {
                criterion = 100;
            }
        }
        if (m_level == 5)
        {
            const Coordinates<data_t> coords(current_cell, m_dx);
            if (coords.x < 0.015625*L && coords.y < 0.015625*L && coords.z < 0.015625*L)
            {
                criterion = 100;
            }
        }
        if (m_level == 6)
        {
            const Coordinates<data_t> coords(current_cell, m_dx);
            if (coords.x < 0.0078125*L && coords.y < 0.0078125*L && coords.z < 0.0078125*L)
            {
                criterion = 100;
            }
        }
        if (m_level == 7)
        {
            const Coordinates<data_t> coords(current_cell, m_dx);
            if (coords.x < 0.00390625*L && coords.y < 0.00390625*L && coords.z < 0.00390625*L)
            {
                criterion = 100;
            }
        }

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* OSCILLOTONTAGGINGCRITERION_HPP_ */
