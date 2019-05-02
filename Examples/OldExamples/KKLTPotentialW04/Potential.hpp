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
//	double decay_constant;
    };

  private:
    params_t modulus_params;

  public:
    //! The constructor
    Potential(params_t a_params) : modulus_params(a_params) {}

    //! Set the potential function for the scalar field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &V_of_phi, data_t &dVdphi,
                           const vars_t<data_t> &vars) const
    {
        // The potential value at phi
        // W0 = 10-4
        //
	double phi_min = 1.28823196577545;
	double pi = M_PI;
        V_of_phi = modulus_params.overall_normalization * ((5.5699469787e-9*exp(-4.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))
							+ 41.8879020479*exp(-2.0*pi*exp(2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))*exp(-4.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))
							  *(-0.0002999999*exp(pi*exp(2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))+10.0*(3.0+pi*exp(2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))))));

        // The potential gradient at phi
        // 2 V0 a Exp[-a phi] (1-Exp[-a phi])
        dVdphi = modulus_params.overall_normalization * (-1.2863241550e-8*exp(-4.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))
						      + 41.8879020479*exp(-2.0*pi*exp(2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))*exp(-4.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))
							*(36.2759872847*exp(2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))-0.0010882796*exp(pi*exp(2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))*exp(2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))
						      - 96.7359660925*exp(-2.0*pi*exp(2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))*exp(-4.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))
                                                          *(-0.0002999999*exp(pi*exp(2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))+10.0*(3.0+pi*exp(2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))))
						      - 303.9050004140*exp(-2.0*pi*exp(2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))*exp(-2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))
                                                          *(-0.0002999999*exp(pi*exp(2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))+10.0*(3.0+pi*exp(2.0*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))));
	
	//pout() << V_of_phi ;
	}
};

#endif /* POTENTIAL_HPP_ */
