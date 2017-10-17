/*  poisson main source */

// hgf includes
#include "model_poisson.hpp"

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

#define midpoint_2d( mid, first, last ) do { \
  mid[0] = 0.5 * (first[0] + last[0]); \
  mid[1] = 0.5 * (first[1] + last[1]); \
} while(0)

#define midpoint_3d( mid, first, second, third, fourth ) do { \
  mid[0] = 0.25 * (first[0] + second[0] + third[0] + fourth[0]); \
  mid[1] = 0.25 * (first[1] + second[1] + third[1] + fourth[1]); \
  mid[2] = 0.25 * (first[2] + second[2] + third[2] + fourth[2]); \
} while(0)

/** \brief hgf::models::poisson::build builds the degrees of freedom, initializes the solution and rhs vectors, and sets up the linear system for Poisson model.
 *
 * @param[in] par - parameters struct containing problem information.
 * @param[in] msh - mesh object containing a quadrilateral or hexagonal representation of geometry from problem folder addressed in parameters& par.
 */
void
hgf::models::poisson::build(const parameters& par, const hgf::mesh::voxel& msh)
{

  phi.resize(msh.els.size());
  NTHREADS = omp_get_max_threads();
  block_size = ((int)phi.size() % NTHREADS) ? (int)((phi.size() / NTHREADS) + 1) : (int)(phi.size() / NTHREADS);

  if (par.dimension == 2) {

    if (alpha.size() != msh.els.size()) {
      alpha.resize(msh.els.size());
      for (int cell = 0; cell < alpha.size(); cell++) {
        alpha[cell].assign(4, 0);
        alpha[cell][0] = 1;
        alpha[cell][3] = 1;
      }
    }

    // setup the degrees of freedom
    build_degrees_of_freedom_2d(par, msh);

    // initialize solution and rhs
    solution.resize(phi.size());
    rhs.resize(phi.size());

    // setup the linear system
    build_array_2d(par, msh);

  }
  else {

    if (alpha.size() != msh.els.size()) {
      alpha.resize(msh.els.size());
      for (int cell = 0; cell < alpha.size(); cell++) {
        alpha[cell].assign(9, 0);
        alpha[cell][0] = 1;
        alpha[cell][4] = 1;
        alpha[cell][8] = 1;
      }
    }

    // setup the degrees of freedom
    build_degrees_of_freedom_3d(par, msh);

    // initialize solution and rhs
    solution.resize(phi.size());
    rhs.resize(phi.size());

    // setup the linear system
    build_array_3d(par, msh);

  }
}

/** \brief hgf::models::poisson::setup_dirichlet_bc setups up the boundary conditions for Dirichlet boundary conditions. If nothing is added to the force,
 * the results is homogeneous Dirichlet BCs on the entire boundary. For nonhomogeneous conditions, add_nonhomogeneous_dirichlet_bc should be called subsequently.
 *
 * Contributions to the linear system coo_array and the rhs vector are set by this function.
 * @param[in] par - parameters struct containing problem information.
 * @param[in] msh - mesh object containing a quadrilateral or hexagonal representation of geometry from problem folder addressed in parameters& par.
 */
void
hgf::models::poisson::setup_dirichlet_bc(const parameters& par, const hgf::mesh::voxel& msh)
{
  if (par.dimension == 2) homogeneous_dirichlet_2d(par, msh);
  else homogeneous_dirichlet_3d(par, msh);
}

/** \brief hgf::models::poisson::setup_mixed_bc setups up the boundary conditions for mixed Dirichlet and Neumann boundary conditions. If nothing is added to the force,
 * the results is homogeneous BCs on the entire boundary. For nonhomogeneous conditions, add_nonhomogeneous_dirichlet_bc and add_nonhomogeneous_neumann_bc 
 * should be called subsequently.
 *
 * Contributions to the linear system coo_array and the rhs vector are set by this function.
 * @param[in] par - parameters struct containing problem information.
 * @param[in] msh - mesh object containing a quadrilateral or hexagonal representation of geometry from problem folder addressed in parameters& par.
 * @param[in] f - pointer to heuristic function. Heuristic should take a cell index and coordinates as inputs, 
 *                and return true if the location has a dirichlet bc or false if the location has a neumann bc.
 * 
 */
void
hgf::models::poisson::setup_mixed_bc(const parameters& par, const hgf::mesh::voxel& msh, bool (*f)( const parameters& par, int dof_num, double coords[3] ))
{
  if (par.dimension == 2) homogeneous_mixed_2d(par, msh, f);
  else homogeneous_mixed_3d(par, msh, f);
}