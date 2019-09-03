/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CHIANDPHITAGGINGCRITERION_HPP_
#define CHIANDPHITAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "ScalarField.hpp"
#include "SimulationParametersBase.hpp"
#include "Tensor.hpp"

class ChiAndPhiTaggingCriterion
{
  protected:
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    const double m_threshold_chi;
    const double m_threshold_phi;
    const int m_level;
    const extraction_params_t m_params;

    template <class data_t>
    using MatterVars = typename GeometricFields::template Vars<data_t>;

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
                              const double threshold_phi, const int a_level,
                              const extraction_params_t a_params)
        : m_dx(dx), m_deriv(dx), m_threshold_chi(threshold_chi),
          m_threshold_phi(threshold_phi), m_params(a_params),
          m_level(a_level){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto d1 = m_deriv.template diff1<MatterVars>(current_cell);
        const auto d1chi = m_deriv.template diff1<Vars>(current_cell);
        const auto d2 = m_deriv.template diff2<MatterVars>(current_cell);
        const auto d2chi = m_deriv.template diff2<Vars>(current_cell);

        data_t mod_d2_chi = 0;
        data_t mod_d2_phi = 0;

        FOR2(idir, jdir)
        {
            mod_d2_chi += d2chi.chi[idir][jdir] * d2chi.chi[idir][jdir] /
                          (1e-2 + abs(d1chi.chi[idir] * d1chi.chi[jdir]));

            mod_d2_phi += d2.Pi1[idir][jdir] * d2.Pi1[idir][jdir] /
                              (abs(d1.Pi1[idir] * d1.Pi1[jdir]) + 3e-6) +
                          d2.phi1[idir][jdir] * d2.phi1[idir][jdir] /
                              (abs(d1.phi1[idir] * d1.phi1[jdir]) + 3e-6) +
                          d2.Pi2[idir][jdir] * d2.Pi2[idir][jdir] /
                              (abs(d1.Pi2[idir] * d1.Pi2[jdir]) + 3e-6) +
                          d2.phi2[idir][jdir] * d2.phi2[idir][jdir] /
                              (abs(d1.phi2[idir] * d1.phi2[jdir]) + 3e-6);
        }

        data_t criterion_chi = m_dx / m_threshold_chi * sqrt(mod_d2_chi);

        data_t criterion_phi = m_dx / m_threshold_phi * sqrt(mod_d2_phi);

        data_t criterion = simd_max(criterion_chi, criterion_phi);

        // regrid if within extraction level and not at required refinement
//        if (m_level < m_params.extraction_level)
//        {
//            const Coordinates<data_t> coords(current_cell, m_dx,
//                                             m_params.extraction_center);
//            const data_t r = coords.get_radius();
            // add a 20% buffer to extraction zone so not too near to boundary
//            data_t regrid =
//                simd_compare_lt(r, 1.2 * m_params.extraction_radius);
//            criterion = simd_conditional(regrid, 1.0, criterion);
//        }

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* CHIANDPHITAGGINGCRITERION_HPP_ */
