/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include <iostream>
#include "parstream.H"

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include "simd.hpp"

class Potential
{
  public:
    struct params_t
    {
        double overall_normalization;
	double decay_constant;
    };

  private:
    params_t m_params;

  public:
    //! The constructor
    Potential(params_t a_params) : m_params(a_params) {}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi,
                           const vars_t<data_t> &vars) const
    {
        // The potential value at phi
        // V0(1 - Exp(-a phi))^2
        float pi = M_PI;
        V_of_phi = pow(m_params.overall_normalization,2.0)*pow((16.0*pi),-1.0)*pow(tanh(sqrt(8.0*pi)*pow(m_params.overall_normalization,-1.0)*vars.phi),2.0);

        // The potential gradient at phi
        // 2 V0 a Exp[-a phi] (1-Exp[-a phi])
        dVdphi = m_params.overall_normalization*pow(sqrt(8.0*pi),-1.0)*pow(cosh(sqrt(8.0*pi)*pow(m_params.overall_normalization,-1.0)*vars.phi),-2.0)*tanh(sqrt(8.0*pi)*pow(m_params.overall_normalization,-1.0)*vars.phi);
    }
};

#endif /* POTENTIAL_HPP_ */
