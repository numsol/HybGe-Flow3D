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

/** \brief hgf::models::poisson::set_nonhomogeneous_dirichlet_bc sets the force vector for a nonhomogeneous Dirichlet condition.
 * 
 * This function uses a user-defined heurstic function taken as a function pointer argument to assign dirichlet values in the force vector.
 * This should be called after setup_homogeneous_dirichlet_bc, which sets the array entries for Dirichlet conditions. 
 * The function adds values to the force vector according to the heuristic. It does not over-write existing values. 
 * @param[in] par - parameters struct containing problem information.
 * @param[in] msh - mesh object containing a quadrilateral or hexagonal representation of geometry from problem folder addressed in parameters& par.
 * @param[in] f - pointer to heuristic function. Heuristic should take a cell index and coordinates as inputs, 
 *                and return a Dirichlet BC value based on one or both of those inputs.
 * 
 */
void 
hgf::models::poisson::set_nonhomogeneous_dirichlet_bc(const parameters& par, const hgf::mesh::voxel& msh, double (*f)( int dof_num, double coords[3] ) )
{

  if (par.dimension == 2) {

#pragma omp parallel for schedule(dynamic) num_threads(NTHREADS)
    for (int kk = 0; kk < NTHREADS; kk++) { 
      // loop over cells
      for (int cell = kk*block_size; cell < std::min((kk + 1)*block_size, (int)phi.size()); cell++) {

        double dx, dy; 
        double value = 0;
        double coords[3];
        int nbrs[4];
        int bc_contributor[4] = { 0, 0, 0, 0 };
        int nnbr = 0;
        for (int jj = 0; jj < 4; jj++) {
          nbrs[jj] = phi[cell].neighbors[jj];
          if (nbrs[jj] != -1) nnbr++;
          else bc_contributor[jj] = 1;
        }
        if (nnbr == 4) goto bcexit2;

        dx = msh.els[cell].vtx[1].coords[0] - msh.els[cell].vtx[0].coords[0];
        dy = msh.els[cell].vtx[3].coords[1] - msh.els[cell].vtx[0].coords[1];

        // S neighbor?
        if (bc_contributor[0]) {
          midpoint_2d(coords, msh.els[cell].vtx[0].coords, msh.els[cell].vtx[1].coords);
          value += 2 * f( cell, coords ) * dx / dy;
        } 
      
        // E neighbor?
        if (bc_contributor[1]) {
          midpoint_2d(coords, msh.els[cell].vtx[1].coords, msh.els[cell].vtx[2].coords);
          value += 2 * f( cell, coords ) * dy / dx;
        }
      
        // N neighbor?
        if (bc_contributor[2]) {
          midpoint_2d(coords, msh.els[cell].vtx[2].coords, msh.els[cell].vtx[3].coords);
          value += 2 * f( cell, coords ) * dx / dy;
        }
      
        // W neighbor?
        if (bc_contributor[3]) {
          midpoint_2d(coords, msh.els[cell].vtx[0].coords, msh.els[cell].vtx[3].coords);
          value += 2 * f( cell, coords ) * dy / dx;
        }

        rhs[cell] += value;

      bcexit2:;
      }
    }

  } else { // 3d

#pragma omp parallel for schedule(dynamic) num_threads(NTHREADS)
    for (int kk = 0; kk < NTHREADS; kk++) { 
      // loop over cells
      for (int cell = kk*block_size; cell < std::min((kk + 1)*block_size, (int)phi.size()); cell++) {
        
        double value = 0;
        double coords[3];
        int nbrs[6];
        int bc_contributor[6] = { 0, 0, 0, 0, 0, 0 };
        int nnbr = 0;
        double dx, dy, dz;
        for (int jj = 0; jj < 6; jj++) {
          nbrs[jj] = phi[cell].neighbors[jj];
          if (nbrs[jj] != -1) nnbr++;
          else bc_contributor[jj] = 1;
        }
        if (nnbr == 6) goto bcexit3;

        dx = msh.els[cell].vtx[1].coords[0] - msh.els[cell].vtx[0].coords[0];
        dy = msh.els[cell].vtx[3].coords[1] - msh.els[cell].vtx[0].coords[1];
        dz = msh.els[cell].vtx[7].coords[2] - msh.els[cell].vtx[0].coords[2];
  
        // y- neighbor?
        if (bc_contributor[0]) {
          midpoint_3d(coords, msh.els[cell].vtx[0].coords, msh.els[cell].vtx[1].coords, \
                              msh.els[cell].vtx[6].coords, msh.els[cell].vtx[7].coords);
          value += 2 * f( cell, coords ) * dx * dz / dy;
        }
  
        // x+ neighbor?
        if (bc_contributor[1]) {
          midpoint_3d(coords, msh.els[cell].vtx[1].coords, msh.els[cell].vtx[2].coords, \
                              msh.els[cell].vtx[5].coords, msh.els[cell].vtx[6].coords);
          value += 2 * f( cell, coords ) * dy * dz / dx;
        }
  
        // y+ neighbor?
        if (bc_contributor[2]) {
          midpoint_3d(coords, msh.els[cell].vtx[2].coords, msh.els[cell].vtx[3].coords, \
                              msh.els[cell].vtx[4].coords, msh.els[cell].vtx[5].coords);
          value += 2 * f( cell, coords ) * dx * dz / dy;
        }
  
        // x- neighbor?
        if (bc_contributor[3]) {
          midpoint_3d(coords, msh.els[cell].vtx[0].coords, msh.els[cell].vtx[3].coords, \
                              msh.els[cell].vtx[4].coords, msh.els[cell].vtx[7].coords);
          value += 2 * f( cell, coords ) * dy * dz / dx;
        }
  
        // z- neighbor?
        if (bc_contributor[4]) {
          midpoint_3d(coords, msh.els[cell].vtx[0].coords, msh.els[cell].vtx[1].coords, \
                              msh.els[cell].vtx[2].coords, msh.els[cell].vtx[3].coords);
          value += 2 * f( cell, coords ) * dx * dy / dz;
        }
  
        // z+ neighbor?
        if (bc_contributor[5]) {
          midpoint_3d(coords, msh.els[cell].vtx[4].coords, msh.els[cell].vtx[5].coords, \
                              msh.els[cell].vtx[6].coords, msh.els[cell].vtx[7].coords);
          value += 2 * f( cell, coords ) * dx * dy / dz;
        }

      bcexit3:;
      }
    }
  }
   
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
