/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GEOMETRICFIELDS_HPP_
#define GEOMETRICFIELDS_HPP_

#include "CCZ4Geometry.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GeometricPotential.hpp"
#include "ScalarField.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS, total num of components
#include "VarsTools.hpp"

//!  Calculates the matter type specific elements such as
//!  the EMTensor and matter evolution for a complex scalar field
/*!
     This class is an example of a matter_t object which calculates the matter
     type specific elements for the RHS update and the evaluation of the
     constraints. This includes the Energy Momentum Tensor, and the matter
     evolution terms. In this case, a pair of fields related by a geometric
   field space.
*/
class GeometricFields
{
  protected:
    //! The local copy of the potential
    GeometricPotential my_potential;

  public:
    //!  Constructor of class GeometricFields
    GeometricFields(const GeometricPotential a_potential)
        : my_potential(a_potential)
    {
    }

    //! Structure containing the rhs variables for the matter fields
    template <class data_t> struct Vars
    {
        data_t phi1;
        data_t Pi1;
        data_t phi2;
        data_t Pi2;

        /// Defines the mapping between members of Vars and Chombo
        /// grid variables (enum in UserVariables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(mapping_function, c_phi1, phi1);
            VarsTools::define_enum_mapping(mapping_function, c_Pi1, Pi1);
            VarsTools::define_enum_mapping(mapping_function, c_phi2, phi2);
            VarsTools::define_enum_mapping(mapping_function, c_Pi2, Pi2);
        }
    };

    //! Structure containing the rhs variables for fields requiring 2nd derivs
    template <class data_t> struct Diff2Vars
    {
        data_t phi1;
        data_t phi2;

        /// Defines the mapping between members of Vars and
        /// Chombo grid variables (enum in UserVariables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(mapping_function, c_phi1, phi1);
            VarsTools::define_enum_mapping(mapping_function, c_phi2, phi2);
        }
    };

    //! The function which calculates the EM Tensor, given the vars and
    //! derivatives
    template <class data_t, template <typename> class vars_t>
    emtensor_t<data_t> compute_emtensor(
        const vars_t<data_t>
            &vars, //!< the value of the variables at the point.
        const vars_t<Tensor<1, data_t>> &d1, //!< the value of the 1st derivs
        const Tensor<2, data_t> &h_UU, //!< the inverse metric (raised indices)
        const Tensor<3, data_t> &chris_ULL //!< the conformal christoffel symbol
        ) const;

    //! The function which adds in the RHS for the matter field vars
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    void add_matter_rhs(
        vars_t<data_t> &total_rhs,  //!< contains the value of the RHS vars.
        const vars_t<data_t> &vars, //!< the value of the vars at the point.
        const vars_t<Tensor<1, data_t>> &d1, //!< the value of the 1st derivs.
        const diff2_vars_t<Tensor<2, data_t>> &d2, //!< the 2nd derivs of vars
        const vars_t<data_t> &advec                //!< the advection terms
        ) const;

    // Calculated geometric terms from fields space metric
    template <class data_t>
    static void calculate_f(const data_t &phi1, data_t &f, data_t &df_dphi1);
};

#include "GeometricFields.impl.hpp"

#endif /* GEOMETRICFIELDS_HPP_ */
