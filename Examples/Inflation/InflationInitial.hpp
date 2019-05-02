/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INFLATIONINITIAL_HPP_
#define INFLATIONINITIAL_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4.hpp"
#include "ScalarField.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"
#include "InflationPotential.hpp"

//! Class which creates a bubble of a scalar field given params for initial
//! matter config
class InflationInitial
{
  public:
    template <class data_t>
    using Vars = MatterCCZ4<ScalarField<>>::Vars<data_t>;

    //! The constructor
    InflationInitial(double a_amplitude_phi, InflationPotential a_potential, 
                     double a_dx, double a_L)
        : m_amplitude_phi(a_amplitude_phi), m_potential(a_potential),
          m_dx(a_dx), m_L(a_L)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    double m_dx;
    double m_L;
    double m_amplitude_phi; //!< Amplitude of fluctuations in phi
    InflationPotential m_potential;

    //! Function to compute the value of K at each point
    template <class data_t>
    data_t compute_K(Vars<data_t> vars) const;
};

#include "InflationInitial.impl.hpp"

#endif /* INFLATIONINITIAL_HPP_ */
