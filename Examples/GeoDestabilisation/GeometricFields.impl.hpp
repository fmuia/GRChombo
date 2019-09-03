/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(GEOMETRICFIELDS_HPP_)
#error "This file should only be included through GeometricFields.hpp"
#endif

#ifndef GEOMETRICFIELDS_IMPL_HPP_
#define GEOMETRICFIELDS_IMPL_HPP_

// Calculate the stress energy tensor elements
// including both fields
// vars = variables on grid
// d1 = spatial derivatives of the fields
// h_UU is spatial metric with raised indices
// chris is spatial metric Christoffel symbol
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> GeometricFields::compute_emtensor(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const Tensor<2, data_t> &h_UU, const Tensor<3, data_t> &chris_ULL) const
{
    // The total em tensor
    emtensor_t<data_t> out;

    // Geometric couplings
    data_t f, df_dphi1;
    calculate_f(vars.phi1, f, df_dphi1);

    // Useful quantity Vt - this is d_mu phi * d^mu phi
    // calculate separately for each field as different geometric
    // factors
    data_t Vt1 = -vars.Pi1 * vars.Pi1;
    data_t Vt2 = -vars.Pi2 * vars.Pi2;
    FOR2(i, j)
    {
        Vt1 += vars.chi * h_UU[i][j] * d1.phi1[i] * d1.phi1[j];
        Vt2 += vars.chi * h_UU[i][j] * d1.phi2[i] * d1.phi2[j];
    }

    // Calculate components of EM Tensor, in 3+1 decomposition
    // S_ij = T_ij
    FOR2(i, j)
    {
        out.Sij[i][j] = -0.5 * 1.0 * vars.h[i][j] * Vt1 / vars.chi +
                        1.0 * d1.phi1[i] * d1.phi1[j] -
                        0.5 * f * f * vars.h[i][j] * Vt2 / vars.chi +
                        f * f * d1.phi2[i] * d1.phi2[j];
    }

    // S = Tr_S_ij
    out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij, h_UU);

    // S_i (note lower index) = - n^a T_ai
    FOR1(i)
    {
        out.Si[i] =
            -1.0 * d1.phi1[i] * vars.Pi1 - f * f * d1.phi2[i] * vars.Pi2;
    }

    // rho = n^a n^b T_ab
    out.rho = 1.0 * vars.Pi1 * vars.Pi1 + f * f * vars.Pi2 * vars.Pi2 +
              0.5 * 1.0 * Vt1 + 0.5 * f * f * Vt2;

    // set the potential values
    data_t V_of_phi = 0.0;
    data_t dVdphi1 = 0.0;
    data_t dVdphi2 = 0.0;

    // compute potential and derivatives
    // currently V = 0.5 m^2 phi1^2
    my_potential.compute_potential(V_of_phi, dVdphi1, dVdphi2, vars);

    // calculate total emtensor including the potential
    out.rho += V_of_phi;
    out.S += -3.0 * V_of_phi;
    FOR2(i, j) { out.Sij[i][j] += -vars.h[i][j] * V_of_phi / vars.chi; }

    return out;
}

// Calculates dphidt and dPidt for the two fields
// and puts it into the total_rhs
// NB dphidt = lapse * Pi + \beta^i d_i \phi
// d1 and d2 are first and second derivs respectively
// advec is \beta^i d_i (var)
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
void GeometricFields::add_matter_rhs(vars_t<data_t> &total_rhs,
                                     const vars_t<data_t> &vars,
                                     const vars_t<Tensor<1, data_t>> &d1,
                                     const diff2_vars_t<Tensor<2, data_t>> &d2,
                                     const vars_t<data_t> &advec) const
{
    // first get the non matter stuff
    // for normal minimally coupled fields
    // (this is well tested, so we reuse it and just add additional
    // terms below, although it makes it a bit ugly)
    ScalarField<>::SFObject<data_t> rhs_sf1;
    ScalarField<>::SFObject<data_t> rhs_sf2;

    // the advection terms
    ScalarField<>::SFObject<data_t> advec_sf1;
    advec_sf1.phi = advec.phi1;
    advec_sf1.Pi = advec.Pi1;
    ScalarField<>::SFObject<data_t> advec_sf2;
    advec_sf2.phi = advec.phi2;
    advec_sf2.Pi = advec.Pi2;

    // the vars
    ScalarField<>::SFObject<data_t> vars_sf1;
    vars_sf1.phi = vars.phi1;
    vars_sf1.Pi = vars.Pi1;
    ScalarField<>::SFObject<data_t> vars_sf2;
    vars_sf2.phi = vars.phi2;
    vars_sf2.Pi = vars.Pi2;

    // go and get the minimally coupled rhs terms
    ScalarField<>::matter_rhs_excl_potential(rhs_sf1, vars, vars_sf1, d1,
                                             d1.phi1, d2.phi1, advec_sf1);
    ScalarField<>::matter_rhs_excl_potential(rhs_sf2, vars, vars_sf2, d1,
                                             d1.phi2, d2.phi2, advec_sf2);

    // set the potential values
    data_t V_of_phi = 0.0;
    data_t dVdphi1 = 0.0;
    data_t dVdphi2 = 0.0;

    // compute potential
    // currently V = 0.5 m^2 phi1^2
    my_potential.compute_potential(V_of_phi, dVdphi1, dVdphi2, vars);

    // adjust RHS for geometric coupling and potential
    data_t f, df_dphi1;
    calculate_f(vars.phi1, f, df_dphi1);

    // calculate inverse of spatial metric
    const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);

    // calculate right hand side (ie dphidt) inc the potential
    total_rhs.phi1 = rhs_sf1.phi;
    total_rhs.Pi1 = rhs_sf1.Pi - vars.lapse * dVdphi1 * 1.0;

    total_rhs.phi2 = rhs_sf2.phi;
    total_rhs.Pi2 = rhs_sf2.Pi - vars.lapse * dVdphi2 / f / f;

    // now add the additional geometric terms
    data_t deriv_coupling_1 = 0.0;
    data_t deriv_coupling_2 = 0.0;
    FOR2(i, j)
    {
        deriv_coupling_1 += vars.chi * h_UU[i][j] * d1.phi2[i] * d1.phi2[j];
        deriv_coupling_2 += vars.chi * h_UU[i][j] * d1.phi1[i] * d1.phi2[j];
    }

    total_rhs.Pi1 +=
        vars.lapse * df_dphi1 * f * (vars.Pi2 * vars.Pi2 - deriv_coupling_1);

    total_rhs.Pi2 += -vars.lapse * 2.0 * df_dphi1 / f *
                     (vars.Pi1 * vars.Pi2 - deriv_coupling_2);
}

// Calculated geometric terms from field space metric
template <class data_t>
void GeometricFields::calculate_f(const data_t &phi1, data_t &f,
                                  data_t &df_dphi1)
{
    double big_M = 0.1 / sqrt(8.0 * M_PI);

    // zero coupling case
    // data_t f_squared = 1.0;
    // f = sqrt(f_squared);
    // df_dphi1 = 0.0;

    // exponential case
    data_t f_squared = exp(-2.0 * phi1 / big_M);
    f = sqrt(f_squared);
    df_dphi1 = -f / big_M;

    // quadratic case
    // data_t f_squared = 1.0 + 2.0 * phi1*phi1 / (big_M*big_M);
    // f = sqrt(f_squared);
    // df_dphi1 = 2.0 * phi1 / f / big_M / big_M;
}

#endif /* GEOMETRICFIELDS_IMPL_HPP_ */
