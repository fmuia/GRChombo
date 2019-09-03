/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// This compute class enforces the positive chi and alpha condition
#ifndef POSITIVECHIANDALPHA2_HPP_
#define POSITIVECHIANDALPHA2_HPP_

#include "Cell.hpp"
#include "UserVariables.hpp"
#include "simd.hpp"

class PositiveChiAndAlpha2
{
  public:
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        auto chi = current_cell.load_vars(c_chi);
        auto lapse = current_cell.load_vars(c_lapse);

        chi = simd_max(chi, 1e-20);
        lapse = simd_max(lapse, 1e-4);

        current_cell.store_vars(chi, c_chi);
        current_cell.store_vars(lapse, c_lapse);
    }
};

#endif /* POSITIVECHIANDALPHA2_HPP_ */
