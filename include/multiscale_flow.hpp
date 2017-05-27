#ifndef _MULTISCALE_FLOW_H
#define _MULTISCALE_FLOW_H

#include "hgflow.hpp"

namespace hgf
{
  /** \brief Contains functions that communicate between scales.
   *
   */
  namespace multiscale
  {
    /** \brief Contains functions that communicate between scales for fluid flow problems.
     * 
     */
    namespace flow
    {
      double compute_permeability_x(const parameters& par, const std::vector< int >& pressure_ib_list, \
                                                           const std::vector< degree_of_freedom >& velocity_u, \
                                                           const std::vector< degree_of_freedom >& velocity_v, \
                                                           const std::vector< degree_of_freedom >& velocity_w, \
                                                           const std::vector< double > solution);
      double compute_permeability_y(const parameters& par, const std::vector< int >& pressure_ib_list, \
                                                           const std::vector< degree_of_freedom >& velocity_u, \
                                                           const std::vector< degree_of_freedom >& velocity_v, \
                                                           const std::vector< degree_of_freedom >& velocity_w, \
                                                           const std::vector< double > solution);
      double compute_permeability_z(const parameters& par, const std::vector< int >& pressure_ib_list, \
                                                           const std::vector< degree_of_freedom >& velocity_u, \
                                                           const std::vector< degree_of_freedom >& velocity_v, \
                                                           const std::vector< degree_of_freedom >& velocity_w, \
                                                           const std::vector< double > solution);
      void compute_permeability_tensor(const parameters& par, const std::vector< int >& pressure_ib_list, \
                                                           const std::vector< degree_of_freedom >& velocity_u, \
                                                           const std::vector< degree_of_freedom >& velocity_v, \
                                                           const std::vector< degree_of_freedom >& velocity_w, \
                                                           const std::vector< double > solution_xflow, \
                                                           const std::vector< double > solution_yflow, \
                                                           const std::vector< double > solution_zflow, \
                                                           std::vector< double >& permeability);
    }
  }
}

#endif
