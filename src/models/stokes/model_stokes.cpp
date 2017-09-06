/* stokes main source */

#include <ctime>
// hgf includes
#include "model_stokes.hpp"

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

/** \brief hgf::models::stokes::build builds the degrees of freedom, initializes the solution and rhs vectors, and sets up the linear system for Stokes flow.
 *
 * Immersed boundary and boundary condition information are not set by this function.
 * @param[in] par - parameters struct containing problem information.
 * @param[in] msh - mesh object containing a quadrilateral or hexagonal representation of geometry from problem folder addressed in parameters& par.
 */
void
hgf::models::stokes::build(const parameters& par, const hgf::mesh::voxel& msh)
{
  // seed rand in case of IB generation
  srand(time(NULL));

  // viscosity is initialized to 1.0
  viscosity = 1.0;

  if (par.dimension == 2) {

    // setup the degrees of freedom
    build_degrees_of_freedom_2d(par, msh);

    // initialize solution and rhs
    int nU = std::accumulate(interior_u.begin(), interior_u.end(), 0);
    int nV = std::accumulate(interior_v.begin(), interior_v.end(), 0);
    int nP = (int)pressure.size();
    solution_int.resize(nU + nV + nP);
    rhs.resize(nU + nV + nP);
    pressure_ib_list.assign(nP, 0);

    // setup the linear system
    build_array_2d(par, msh);

  }
  else {

    // setup the degrees of freedom
    build_degrees_of_freedom_3d(par, msh);

    // initialize solution and rhs
    int nU = std::accumulate(interior_u.begin(), interior_u.end(), 0);
    int nV = std::accumulate(interior_v.begin(), interior_v.end(), 0);
    int nW = std::accumulate(interior_w.begin(), interior_w.end(), 0);
    int nP = (int)pressure.size();
    solution_int.resize(nU + nV + nW + nP);
    rhs.resize(nU + nV + nW + nP);
    pressure_ib_list.assign(nP, 0);

    // setup the linear system
    build_array_3d(par, msh);

  }
#ifdef _ARRAY_DEBUG
  std::cout << "\nArray size = " << coo_array.size() << "\n";
  std::cout << "\nnU = " << std::accumulate(interior_u.begin(), interior_u.end(), 0) << ",\tnV = " << std::accumulate(interior_v.begin(), interior_v.end(), 0);
  if (par.dimension == 3) {
    std::cout << ",\tNW = " << std::accumulate(interior_w.begin(), interior_w.end(), 0);
  }
  std::cout << ",\tnP = " << (int)pressure.size() << "\n";
  for (int ii = 0; ii < coo_array.size(); ii++) {
    std::cout << coo_array[ii].i_index << "\t" << coo_array[ii].j_index << "\t" << coo_array[ii].value << "\n";
  }
  std::cout << "\n";
#endif
}

/** \brief hgf::models::stokes::setup_xflow_bc setups up the boundary conditions for a problem with flow in the positive direction along the x-axis.
 *
 * Contributions to the linear system coo_array and the rhs vector are set by this function.
 * @param[in] par - parameters struct containing problem information.
 * @param[in] msh - mesh object containing a quadrilateral or hexagonal representation of geometry from problem folder addressed in parameters& par.
 */
void
hgf::models::stokes::setup_xflow_bc(const parameters& par, const hgf::mesh::voxel& msh)
{

  if (par.dimension == 2) xflow_2d(par, msh);
  else xflow_3d(par, msh);

}

/** \brief hgf::models::stokes::setup_yflow_bc setups up the boundary conditions for a problem with flow in the positive direction along the y-axis.
 *
 * Contributions to the linear system coo_array and the rhs vector are set by this function.
 * @param[in] par - parameters struct containing problem information.
 * @param[in] msh - mesh object containing a quadrilateral or hexagonal representation of geometry from problem folder addressed in parameters& par.
 */
void
hgf::models::stokes::setup_yflow_bc(const parameters& par, const hgf::mesh::voxel& msh)
{

  if (par.dimension == 2) yflow_2d(par, msh);
  else yflow_3d(par, msh);

}

/** \brief hgf::models::stokes::setup_zflow_bc setups up the boundary conditions for a problem with flow in the positive direction along the z-axis.
 *
 * Contributions to the linear system coo_array and the rhs vector are set by this function.
 * @param[in] par - parameters struct containing problem information.
 * @param[in] msh - mesh object containing a quadrilateral or hexagonal representation of geometry from problem folder addressed in parameters& par.
 */
void
hgf::models::stokes::setup_zflow_bc(const parameters& par, const hgf::mesh::voxel& msh)
{

  zflow_3d(par, msh);

}


