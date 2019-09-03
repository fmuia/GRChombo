/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "ScalarFieldLevel.hpp"
#include "BoxLoops.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha2.hpp"
#include "TraceARemoval.hpp"

// For RHS update
#include "MatterCCZ4.hpp"

// For constraints calculation
#include "MatterConstraints.hpp"

// For tag cells
#include "ChiAndPhiTaggingCriterion.hpp"

// Problem specific includes
#include "ComputePack.hpp"
#include "GeometricFields.hpp"
#include "GeometricInitialConditions.hpp"
#include "GeometricPotential.hpp"
#include "SetValue.hpp"

// Things to do at each advance step, after the RK4 is calculated
void ScalarFieldLevel::specificAdvance()
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha2()),
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
    // constraints Then the physical setup
    GeometricPotential potential(m_p.potential_params);
    BoxLoops::loop(make_compute_pack(SetValue(0.0),
                                     GeometricInitialConditions(
                                         m_p.initial_params, potential, m_dx)),
                   m_state_new, m_state_new, FILL_GHOST_CELLS);
}

// Things to do before outputting a checkpoint file
void ScalarFieldLevel::preCheckpointLevel()
{
    fillAllGhosts();
    GeometricPotential potential(m_p.potential_params);
    GeometricFields geometric_fields(potential);
    BoxLoops::loop(MatterConstraints<GeometricFields>(geometric_fields, m_dx,
                                                      m_p.G_Newton),
                   m_state_new, m_state_new, SKIP_GHOST_CELLS);
}

// Things to do in RHS update, at each RK4 step
void ScalarFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha2()),
                   a_soln, a_soln, FILL_GHOST_CELLS);

    // Calculate MatterCCZ4 right hand side with matter_t = GeometricFields
    // We don't want undefined values floating around in the constraints so
    // zero these
    GeometricPotential potential(m_p.potential_params);
    GeometricFields geometric_fields(potential);
    MatterCCZ4<GeometricFields> my_ccz4_matter(geometric_fields,
                                               m_p.ccz4_params, m_dx, m_p.sigma,
                                               m_p.formulation, m_p.G_Newton);
    SetValue set_constraints_zero(0.0, Interval(c_Ham, c_Mom3));
    auto compute_pack2 =
        make_compute_pack(my_ccz4_matter, set_constraints_zero);
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
    plot_states = {c_phi1, c_phi2, c_chi, c_lapse, c_K};
}

void ScalarFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
    BoxLoops::loop(ChiAndPhiTaggingCriterion(m_dx, m_p.regrid_threshold_chi,
                                             m_p.regrid_threshold_phi, m_level,
                                             m_p.extraction_params),
                   current_state, tagging_criterion);
}
