/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INFLATIONPOTENTIAL_HPP_
#define INFLATIONPOTENTIAL_HPP_

#include "simd.hpp"

//! The basic idea is that the potential is divided into three sections
//! The reheating part, the flat part to the left and the reheating minimum
//! which is roughly m^2 phi^2 to the right.
class InflationPotential
{
  public:
    struct params_t
    {
        double Lambda;
        double mu;
    };

  private:
    params_t m_params;

  public:
    //! The constructor
    InflationPotential(params_t a_params) : m_params(a_params) {}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi,
                           const vars_t<data_t> &vars) const
    {
        // note the values below work for the params in this example, but you
        // would need to check for others that they give a smooth potential
        // where the sections join...
        data_t phi_here = vars.phi;
        data_t ratio = phi_here / m_params.mu;
        double mpl = 1.0e-9;
        double phi_flat = 0.0;
        double phi_reheat = 0.95 * m_params.mu;
        double phi_reheat_min = 1.05 * m_params.mu;
        double lambda_reheat = 3.75e-18;

        //  inflationary potentials
        auto phi_lt_phi_flat = simd_compare_lt(phi_here, phi_flat);
        auto phi_lt_phi_reheat = simd_compare_lt(phi_here, phi_reheat);

        double V_of_phi_flat = m_params.Lambda;
        data_t V_of_phi_inflation = m_params.Lambda * (1.0 - pow(ratio, 4.0));
        data_t V_of_phi_reheat =
            lambda_reheat * pow((phi_here - phi_reheat_min) / mpl, 2.0);

        double dVdphi_flat = 0.0;
        data_t dVdphi_inflation =
            -4.0 / m_params.mu * m_params.Lambda * pow(ratio, 3.0);
        data_t dVdphi_reheat =
            2.0 * lambda_reheat * (phi_here - phi_reheat_min) * pow(mpl, -2.0);

        // set V_of_phi and dVdphi
        V_of_phi = simd_conditional(phi_lt_phi_flat, V_of_phi_flat,
                                    V_of_phi_inflation);
        V_of_phi =
            simd_conditional(phi_lt_phi_reheat, V_of_phi, V_of_phi_reheat);
        dVdphi =
            simd_conditional(phi_lt_phi_flat, dVdphi_flat, dVdphi_inflation);
        dVdphi = simd_conditional(phi_lt_phi_reheat, dVdphi, dVdphi_reheat);
    }
};

#endif /* POTENTIAL_HPP_ */
