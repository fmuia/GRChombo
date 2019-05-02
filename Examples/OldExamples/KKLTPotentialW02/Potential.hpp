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
        // W0 = 10-4
        //
	double phi_min = 0.883851068963227;
	double pi = M_PI;
        V_of_phi = m_params.overall_normalization * (5.539562263296e-7*exp(-8.0*sqrt(pi)*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))
							+ 0.001666666666666*exp(-0.2*exp(4.0*sqrt(pi)*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))*exp(-8.0*sqrt(pi)*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))
							  *(0.3-0.02999999*exp(0.1*exp(4.0*sqrt(pi)*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))+0.01*exp(4.0*sqrt(pi)*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))));

        // The potential gradient at phi
        // 2 V0 a Exp[-a phi] (1-Exp[-a phi])
        dVdphi = m_params.overall_normalization * (-4.5350256114910e-6*exp(-8.0*sqrt(pi)*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))
						      - 0.01364435610595318*exp(-0.2*exp(4.0*sqrt(pi)*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))*exp(-8.0*sqrt(pi)*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))
							*(0.3-0.03*exp(0.1*exp(4.0*sqrt(pi)*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))+0.01*exp(4.0*sqrt(pi)*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))
						      - 0.001364435610595318*exp(-0.2*exp(4.0*sqrt(pi)*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))*exp(-4.0*sqrt(pi)*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))
                                                        *(0.3-0.03*exp(0.1*exp(4.0*sqrt(pi)*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))+0.01*exp(4.0*sqrt(pi)*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))
						      + 0.0016666666666666666*exp(-0.2*exp(4.0*sqrt(pi)*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))*exp(-8.0*sqrt(pi)*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))
                                                          *(0.0409330683178595*exp(4.0*sqrt(pi)*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))-0.01227992049*exp(0.1*exp(4.0*sqrt(pi)*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi)))*exp(4.0*sqrt(pi)*pow(sqrt(3.0),-1.0)*(phi_min+vars.phi))));
	
	//pout() << V_of_phi ;
	}
};

#endif /* POTENTIAL_HPP_ */
