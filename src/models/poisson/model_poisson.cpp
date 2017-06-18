/* poisson main source */

#include "model_poisson.hpp"

/** \brief hgf::models::poisson::build builds the degrees of freedom, initializes the solution and rhs vectors, and sets up the linear system for the Poisso equation.
 *
 * Immersed boundary and boundary condition information are not set by this function.
 * @param[in] par - parameters struct containing problem information.
 * @param[in] msh - mesh object containing a quadrilateral or hexagonal representation of geometry from problem folder addressed in parameters& par.
 */
 void
 hgf::models::poisson::build(const parameters& par, const hgf::mesh::voxel& msh)

  if (par.dimension == 2) {

    // setup the degrees of freedom
    build_degrees_of_freedom_2d(par, msh);

    // initialize solution and rhs
    solution.resize(cc_dof.size());
    rhs.resize(cc_dof.size());

    // setup the linear system
    build_array_2d(par, msh);

  } else {

    // setup the degrees of freedom
    build_degrees_of_freedom_3d(par, msh);

    // initialize solution and rhs
    solution.resize(cc_dof.size());
    rhs.resize(cc_dof.size());

    // setup the linear system
    build_array_3d(par, msh);

  }

}