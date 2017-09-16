/*  poisson main source */

// hgf includes
#include "model_poisson.hpp"

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

/** \brief hgf::models::poisson::build builds the degrees of freedom, initializes the solution and rhs vectors, and sets up the linear system for Poisson model.
 *
 * @param[in] par - parameters struct containing problem information.
 * @param[in] msh - mesh object containing a quadrilateral or hexagonal representation of geometry from problem folder addressed in parameters& par.
 */
void
hgf::models::poisson::build(const parameters& par, const hgf::mesh::voxel& msh)
{

  if (par.dimension == 2) {

    if (alpha.size() != msh.els.size()) {
      std::cout << "\nSetting alpha...\n";
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

/** \brief hgf::models::poisson::setup_homogeneous_dirichlet_bc setups up the boundary conditions for ... 
 *
 * Contributions to the linear system coo_array and the rhs vector are set by this function.
 * @param[in] par - parameters struct containing problem information.
 * @param[in] msh - mesh object containing a quadrilateral or hexagonal representation of geometry from problem folder addressed in parameters& par.
 */
void
hgf::models::poisson::setup_homogeneous_dirichlet_bc(const parameters& par, const hgf::mesh::voxel& msh)
{
  if (par.dimension == 2) homogeneous_dirichlet_2d(par, msh);
  else homogeneous_dirichlet_3d(par, msh);
}

/** \brief hgf::models::poisson::set_constant_force sets a constant value to the force (right hand side) in the Poisson model.
*
* This function replaces existing values in the force vector by the new value.
* @param[in] force_in - double precision floating point value assigned as the force in all cells. 
*/
void
hgf::models::poisson::set_constant_force(const parameters& par, const double& force_in) 
{
  rhs.assign( phi.size(), force_in );  
}

/** \brief hgf::models::poisson::set_constant_scalar_alpha sets a constant scalar value to the alpha coefficient in the Poisson model.
*
* This function replaces existing values in the alpha tensor by the new constant scalar value.
* @param[in] alpha_in - double precision floating point value assigned as constant scalar alpha.
*/
void
hgf::models::poisson::set_constant_scalar_alpha(const parameters& par, const hgf::mesh::voxel& msh, const double& alpha_in)
{
  std::vector< double > temp;
  temp.assign( par.dimension * par.dimension, 0 );
  temp[0] = alpha_in;
  if (par.dimension == 2) temp[3] = alpha_in;
  else {
    temp[4] = alpha_in;
    temp[8] = alpha_in;
  }
  alpha.assign( msh.els.size(), temp );
}

/** \brief hgf::models::poisson::set_constant_tensor_alpha sets a constant tensor value to the alpha coefficient in the Poisson model.
*
* This function replaces existing values in the alpha tensor by the new constant tensor.
* @param[in] alpha_in - 2x2 or 3x3 vector of double precision floating point values giving holding the constant alpha tensor value in row-major format..
*/
void
hgf::models::poisson::set_constant_tensor_alpha(const parameters& par, const hgf::mesh::voxel& msh, const std::vector< double >& alpha_in)
{
  alpha.assign( msh.els.size(), alpha_in );
}
