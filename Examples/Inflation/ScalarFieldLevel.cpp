/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "ScalarFieldLevel.hpp"
#include "BoxLoops.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "TraceARemoval.hpp"

// For RHS update
#include "MatterCCZ4.hpp"

// For constraints calculation
#include "MatterConstraints.hpp"

// For tag cells
#include "PhiAndRhoAijTaggingCriterion.hpp"

// Problem specific includes
#include "ChiRelaxation.hpp"
#include "ComputePack.hpp"
#include "Density.hpp"
#include "DensityAij.hpp"
#include "InflationInitial.hpp"
#include "InflationPotential.hpp"
#include "ProperTime.hpp"
#include "ScalarField.hpp"
#include "SetValue.hpp"

// Things to do at each advance step, after the RK4 is calculated
void ScalarFieldLevel::specificAdvance()
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
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
    // constraints etc, then  initial conditions for scalar field - here a
    // bubble
    BoxLoops::loop(
        make_compute_pack(SetValue(0.0),
                          InflationInitial(m_p.initial_params, m_dx, m_p.L)),
        m_state_new, m_state_new, FILL_GHOST_CELLS);
}

// Things to do before outputting a checkpoint file
void ScalarFieldLevel::preCheckpointLevel()
{
    fillAllGhosts();
    InflationPotential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    Density<ScalarFieldWithPotential> my_density(scalar_field, m_dx);
    DensityAij my_densityAij(m_p.G_Newton);
    MatterConstraints<ScalarFieldWithPotential> my_constraints(
        scalar_field, m_dx, m_p.G_Newton);
    auto compute_pack =
        make_compute_pack(my_density, my_densityAij, my_constraints);
    BoxLoops::loop(compute_pack, m_state_new, m_state_new, SKIP_GHOST_CELLS);
}

// Things to do in RHS update, at each RK4 step
void ScalarFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{

    // Relaxation function for chi - this will eventually be done separately
    // with hdf5 as input
    if (m_time < m_p.relaxtime)
    {
        // Calculate chi relaxation right hand side
        // Note this assumes conformal chi and Mom constraint trivially
        // satisfied  No evolution in other variables, which are assumed to
        // satisfy constraints per initial conditions
        InflationPotential potential(m_p.potential_params);
        ScalarFieldWithPotential scalar_field(potential);
        ChiRelaxation<ScalarFieldWithPotential> relaxation(
            scalar_field, m_dx, m_p.relaxspeed, m_p.G_Newton);
        SetValue set_other_values_zero(0.0, Interval(c_chi, c_Mom3));
        auto compute_pack1 =
            make_compute_pack(set_other_values_zero, relaxation);
        BoxLoops::loop(compute_pack1, a_soln, a_rhs, SKIP_GHOST_CELLS);
    }
    else
    {
        // Enforce trace free A_ij and positive chi and alpha
        BoxLoops::loop(
            make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()), a_soln,
            a_soln, FILL_GHOST_CELLS);

        // Calculate MatterCCZ4 right hand side with matter_t = ScalarField
        // We don't want undefined values floating around in the constraints so
        // zero these, and diagnostics rho and rhoAij
        InflationPotential potential(m_p.potential_params);
        ScalarFieldWithPotential scalar_field(potential);
        MatterCCZ4<ScalarFieldWithPotential> my_ccz4_matter(
            scalar_field, m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation,
            m_p.G_Newton);
        ProperTime my_proper_time(m_dx);
        SetValue set_diagnostics_zero(0.0, Interval(c_rho, c_Mom3));
        auto compute_pack2 = make_compute_pack(my_ccz4_matter, my_proper_time,
                                               set_diagnostics_zero);
        BoxLoops::loop(compute_pack2, a_soln, a_rhs, SKIP_GHOST_CELLS);
    }
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
    plot_states = {c_phi, c_K};
}

void ScalarFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
    BoxLoops::loop(PhiAndRhoAijTaggingCriterion(m_dx, m_p.regrid_threshold_phi,
                                                m_p.regrid_threshold_rhoAij),
                   current_state, tagging_criterion);
}
