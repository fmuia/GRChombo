/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "ScalarFieldLevel.hpp"
#include "BoxLoops.hpp"
#include "NanCheck.hpp"
//#include "PositiveChiAndAlpha.hpp"
#include "TraceARemoval.hpp"

// For RHS update
#include "MatterCCZ4.hpp"

// For constraints calculation
#include "MatterConstraints.hpp"

// For tag cells
#include "ChiAndPhiTaggingCriterion.hpp"

// Problem specific includes
#include "ComputePack.hpp"
#include "Density.hpp"
#include "InflationInitial.hpp"
#include "InflationPotential.hpp"
#include "ScalarField.hpp"
#include "SetValue.hpp"

// Things to do at each advance step, after the RK4 is calculated
void ScalarFieldLevel::specificAdvance()
{
    // Enforce trace free A_ij
//    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
    BoxLoops::loop(TraceARemoval(),
                   m_state_new, m_state_new, FILL_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new, SKIP_GHOST_CELLS,
                       disable_simd());
}

// Initial data for field and metric variables
void ScalarFieldLevel::initialData()
{
    CH_TIME("ScalarFieldLevel::initialData");
    if (m_verbosity)
        pout() << "ScalarFieldLevel::initialData " << m_level << endl;

    // First set everything to zero ... we don't want undefined values in
    // constraints etc, then  initial conditions for scalar field
    InflationPotential potential(m_p.potential_params);
    BoxLoops::loop(
        make_compute_pack(SetValue(0.0),
                          InflationInitial(m_p.amplitude_phi, potential, m_dx, m_p.L)),
        m_state_new, m_state_new, FILL_GHOST_CELLS);
}

// Things to do before outputting a checkpoint file
void ScalarFieldLevel::prePlotLevel()
{
    fillAllGhosts();
    InflationPotential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    Density<ScalarFieldWithPotential> my_density(scalar_field, m_dx);
    MatterConstraints<ScalarFieldWithPotential> my_constraints(
        scalar_field, m_dx, m_p.G_Newton);
    auto compute_pack =
        make_compute_pack(my_density, my_constraints);
    BoxLoops::loop(compute_pack, m_state_new, m_state_new, SKIP_GHOST_CELLS);
}

// Things to do in RHS update, at each RK4 step
void ScalarFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{
    // Enforce trace free A_ij
    BoxLoops::loop(
        TraceARemoval(), a_soln, a_soln, FILL_GHOST_CELLS);

    // Calculate MatterCCZ4 right hand side with matter_t = ScalarField
    // We don't want undefined values floating around in the constraints so
    // zero these, and diagnostics rho and rhoAij
    InflationPotential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    MatterCCZ4<ScalarFieldWithPotential> my_ccz4_matter(
        scalar_field, m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation,
        m_p.G_Newton);
    SetValue set_diagnostics_zero(0.0, Interval(c_rho, c_Mom3));
    auto compute_pack2 = make_compute_pack(my_ccz4_matter,
                                           set_diagnostics_zero);
    BoxLoops::loop(compute_pack2, a_soln, a_rhs, SKIP_GHOST_CELLS);
}

// Things to do at ODE update, after soln + rhs
void ScalarFieldLevel::specificUpdateODE(GRLevelData &a_soln,
                                         const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce trace free A_ij
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, FILL_GHOST_CELLS);
}

// Specify if you want any plot files to be written, with which vars
void ScalarFieldLevel::specificWritePlotHeader(
    std::vector<int> &plot_states) const
{
    plot_states = {c_phi, c_chi, c_K};
}

void ScalarFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
    BoxLoops::loop(ChiAndPhiTaggingCriterion(m_dx, m_p.regrid_threshold_chi,
                                                m_p.regrid_threshold_phi),
                   current_state, tagging_criterion);
}
