/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GEOMETRICPOTENTIAL_HPP_
#define GEOMETRICPOTENTIAL_HPP_

#include "simd.hpp"

class GeometricPotential
{
  public:
    struct params_t
    {
        double scalar_mass1;
        double scalar_mass2;
    };

  private:
    params_t m_params;

  public:
    //! The constructor
    GeometricPotential(params_t a_params) : m_params(a_params) {}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi1, data_t &dVdphi2,
                           const vars_t<data_t> &vars) const
    {
        // The potential value at phi1, phi2
        V_of_phi =
            0.5 * pow(m_params.scalar_mass1 * vars.phi1, 2.0) / (8.0 * M_PI);
        //         + 0.5 * pow(m_params.scalar_mass2 * vars.phi2, 2.0) /
        //         (8.0*M_PI);

        // The potential gradient at phi wrt the two fields
        dVdphi1 = pow(m_params.scalar_mass1, 2.0) * vars.phi1 / (8.0 * M_PI);
        // dVdphi2 = pow(m_params.scalar_mass2, 2.0) * vars.phi2 / (8.0*M_PI);
        dVdphi2 = 0.0;
    }
};

#endif /* GEOMETRICPOTENTIAL_HPP_ */
