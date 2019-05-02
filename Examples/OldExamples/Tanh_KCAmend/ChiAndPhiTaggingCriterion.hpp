/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CHIANDPHITAGGINGCRITERION_HPP_
#define CHIANDPHITAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "ScalarField.hpp"
#include "Tensor.hpp"

class ChiAndPhiTaggingCriterion
{
  protected:
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    const double m_threshold_chi;
    const double m_threshold_phi;
    const int m_level;

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
    ChiAndPhiTaggingCriterion(const double dx, const double threshold_chi,
                              const double threshold_phi, const int a_level)
        : m_dx(dx), m_deriv(dx), m_threshold_chi(threshold_chi),
          m_threshold_phi(threshold_phi), m_level(a_level){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto d2 = m_deriv.template diff2<MatterVars>(current_cell);
        const auto d2chi = m_deriv.template diff2<Vars>(current_cell);

        data_t mod_d2_chi = 0;
        data_t mod_d2_phi = 0;

        FOR2(idir, jdir)
        {
            mod_d2_chi += d2chi.chi[idir][jdir] * d2chi.chi[idir][jdir];

            mod_d2_phi += d2.Pi[idir][jdir] * d2.Pi[idir][jdir] +
                          d2.phi[idir][jdir] * d2.phi[idir][jdir];
        }

        data_t criterion_chi = m_dx / m_threshold_chi * sqrt(mod_d2_chi);

        data_t criterion_phi = m_dx / m_threshold_phi * sqrt(mod_d2_phi);

        data_t criterion = simd_max(criterion_chi, criterion_phi);

        // make sure the inner part is regridded up to y,z = 128 from the BH
        // path
        if (m_level == 0)
        {
            const Coordinates<data_t> coords(current_cell, m_dx);
            if (coords.x < 128.0 && coords.y < 128.0 && coords.z < 128.0)
            {
                criterion = 100;
            }
        }
        if (m_level == 1)
        {
            const Coordinates<data_t> coords(current_cell, m_dx);
            if (coords.x < 64.0 && coords.y < 64.0 && coords.z < 64.0)
            {
                criterion = 100;
            }
        }
        if (m_level == 2)
        {
            const Coordinates<data_t> coords(current_cell, m_dx);
            if (coords.x < 32.0 && coords.y < 32.0 && coords.z < 32.0)
            {
                criterion = 100;
            }
        }
        if (m_level == 3)
        {
            const Coordinates<data_t> coords(current_cell, m_dx);
            if (coords.x < 16.0 && coords.y < 16.0 && coords.z < 16.0)
            {
                criterion = 100;
            }
        }
        if (m_level == 4)
        {
            const Coordinates<data_t> coords(current_cell, m_dx);
            if (coords.x < 8.0 && coords.y < 8.0 && coords.z < 8.0)
            {
                criterion = 100;
            }
        }
        if (m_level == 5)
        {
            const Coordinates<data_t> coords(current_cell, m_dx);
            if (coords.x < 4.0 && coords.y < 4.0 && coords.z < 4.0)
            {
                criterion = 100;
            }
        }
        if (m_level == 6)
        {
            const Coordinates<data_t> coords(current_cell, m_dx);
            if (coords.x < 2.0 && coords.y < 2.0 && coords.z < 2.0)
            {
                criterion = 100;
            }
        }

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* CHIANDPHITAGGINGCRITERION_HPP_ */
