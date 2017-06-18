/* poisson main source */

#include "model_poisson.hpp"

// 1d->2d index
#define idx(i, j, ldi) ((i * ldi) + j)

/** \brief hgf::models::poisson::build builds the degrees of freedom, initializes the solution and rhs vectors, and sets up the linear system for the Poisso equation.
 *
 * Immersed boundary and boundary condition information are not set by this function.
 * @param[in] par - parameters struct containing problem information.
 * @param[in] msh - mesh object containing a quadrilateral or hexagonal representation of geometry from problem folder addressed in parameters& par.
 */
 void
 hgf::models::poisson::build(const parameters& par, const hgf::mesh::voxel& msh)

  if (par.dimension == 2) {

  } else {

  }

}