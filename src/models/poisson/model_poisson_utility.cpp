/*  poisson utility source */

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

/** \brief hgf::models::poisson::add_nonhomogeneous_bc adds to the force vector the values corersponding to a nonhomogeneous Dirichlet and Neumann boundary conditions.
 * 
 * This function uses a user-defined heurstic function taken as a function pointer argument to assign dirichlet and neumann values in the force vector.
 * This should be called after setup_homogeneous_bc, which sets the array entries for mixed Dirichlet and Neumann boundary conditions and marks
 * edges (2d) or faces (3d) as Dirichlet or Neumann boundaries. 
 * The function adds values to the force vector according to the heuristic. It does not over-write existing values in the force vector. 
 * Note that the Neumann condition includes the alpha coefficient, i.e. the boundary condition that is imposed is "alpha grad u dot n = bc_value".
 * @param[in] par - parameters struct containing problem information.
 * @param[in] msh - mesh object containing a quadrilateral or hexagonal representation of geometry from problem folder addressed in parameters& par.
 * @param[in] bc_value - pointer to heuristic function. Heuristic must take parameters struct, cell index, coordinates, and a bool which can be used to specific Dirichlet or Neuamnn outputs.
 *                       The function then returns the BC value based on one or more of these inputs.
 * 
 */
void 
hgf::models::poisson::add_nonhomogeneous_bc(const parameters& par, const hgf::mesh::voxel& msh, \
                                            double (*bc_value)( const parameters& par, int dof_num, double coords[3] ) )
{

  if (par.dimension == 2) {

#pragma omp parallel for schedule(dynamic) num_threads(NTHREADS)
    for (int kk = 0; kk < NTHREADS; kk++) { 
      // loop over cells
      for (int cell = kk*block_size; cell < std::min((kk + 1)*block_size, (int)phi.size()); cell++) {

        bool alpha_diag;
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

        alpha_diag = (alpha[cell][1] == 0.0 && alpha[cell][2] == 0.0);

        dx = msh.els[cell].vtx[1].coords[0] - msh.els[cell].vtx[0].coords[0];
        dy = msh.els[cell].vtx[3].coords[1] - msh.els[cell].vtx[0].coords[1];

        if (alpha_diag) {

          // S neighbor?
          if (bc_contributor[0]) {
            midpoint_2d(coords, msh.els[cell].vtx[0].coords, msh.els[cell].vtx[1].coords);
            if (bc_types[cell][0] == 1) value += 2 * alpha[cell][3] * bc_value( par, cell, coords ) * dx / dy;
            else value += bc_value( par, cell, coords );
          } 
      
          // E neighbor?
          if (bc_contributor[1]) {
            midpoint_2d(coords, msh.els[cell].vtx[1].coords, msh.els[cell].vtx[2].coords);
            if (bc_types[cell][1] == 1) value += 2 * alpha[cell][0] * bc_value( par, cell, coords ) * dy / dx;
            else value += bc_value( par, cell, coords );
          }
      
          // N neighbor?
          if (bc_contributor[2]) {
            midpoint_2d(coords, msh.els[cell].vtx[2].coords, msh.els[cell].vtx[3].coords);
            if (bc_types[cell][2] == 1) value += 2 * alpha[cell][3] * bc_value( par, cell, coords ) * dx / dy;
            else value += bc_value( par, cell, coords );
          }
      
          // W neighbor?
          if (bc_contributor[3]) {
            midpoint_2d(coords, msh.els[cell].vtx[0].coords, msh.els[cell].vtx[3].coords);
            if (bc_types[cell][3] == 1) value += 2 * alpha[cell][0] * bc_value( par, cell, coords ) * dy / dx;
            else value += bc_value( par, cell, coords );
          }

        }
        else { // Nondiag alpha, TODO

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
        
        bool alpha_diag;
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

        alpha_diag = (alpha[cell][1] == 0.0 && \
          alpha[cell][2] == 0.0 && \
          alpha[cell][3] == 0.0 && \
          alpha[cell][5] == 0.0 && \
          alpha[cell][6] == 0.0 && \
          alpha[cell][7] == 0.0);

        dx = msh.els[cell].vtx[1].coords[0] - msh.els[cell].vtx[0].coords[0];
        dy = msh.els[cell].vtx[3].coords[1] - msh.els[cell].vtx[0].coords[1];
        dz = msh.els[cell].vtx[7].coords[2] - msh.els[cell].vtx[0].coords[2];

        if (alpha_diag) {
  
          // y- neighbor?
          if (bc_contributor[0]) {
            midpoint_3d(coords, msh.els[cell].vtx[0].coords, msh.els[cell].vtx[1].coords, \
                                msh.els[cell].vtx[6].coords, msh.els[cell].vtx[7].coords);
            if (bc_types[cell][0] == 1) value += 2 * alpha[cell][4] * bc_value( par, cell, coords ) * dx * dz / dy;
            else value += bc_value( par, cell, coords );
          }
  
          // x+ neighbor?
          if (bc_contributor[1]) {
            midpoint_3d(coords, msh.els[cell].vtx[1].coords, msh.els[cell].vtx[2].coords, \
                                msh.els[cell].vtx[5].coords, msh.els[cell].vtx[6].coords);
            if (bc_types[cell][1] == 1) value += 2 * alpha[cell][0] * bc_value( par, cell, coords ) * dy * dz / dx;
            else value += bc_value( par, cell, coords );
          }
  
          // y+ neighbor?
          if (bc_contributor[2]) {
            midpoint_3d(coords, msh.els[cell].vtx[2].coords, msh.els[cell].vtx[3].coords, \
                                msh.els[cell].vtx[4].coords, msh.els[cell].vtx[5].coords);
            if (bc_types[cell][2] == 1) value += 2 * alpha[cell][4] * bc_value( par, cell, coords ) * dx * dz / dy;
            else value += bc_value( par, cell, coords );
          }
  
          // x- neighbor?
          if (bc_contributor[3]) {
            midpoint_3d(coords, msh.els[cell].vtx[0].coords, msh.els[cell].vtx[3].coords, \
                                msh.els[cell].vtx[4].coords, msh.els[cell].vtx[7].coords);
            if (bc_types[cell][3] == 1) value += 2 * alpha[cell][0] * bc_value( par, cell, coords ) * dy * dz / dx;
            else value += bc_value( par, cell, coords );
          }
  
          // z- neighbor?
          if (bc_contributor[4]) {
            midpoint_3d(coords, msh.els[cell].vtx[0].coords, msh.els[cell].vtx[1].coords, \
                                msh.els[cell].vtx[2].coords, msh.els[cell].vtx[3].coords);
            if (bc_types[cell][4] == 1) value += 2 * alpha[cell][8] * bc_value( par, cell, coords ) * dx * dy / dz;
            else value += bc_value( par, cell, coords );
          }
  
          // z+ neighbor?
          if (bc_contributor[5]) {
            midpoint_3d(coords, msh.els[cell].vtx[4].coords, msh.els[cell].vtx[5].coords, \
                                msh.els[cell].vtx[6].coords, msh.els[cell].vtx[7].coords);
            if (bc_types[cell][5] == 1) value += 2 * alpha[cell][8] * bc_value( par, cell, coords ) * dx * dy / dz;
            else value += bc_value( par, cell, coords );
          }

        }
        else { // Nondiag alpha, TODO

        }

        rhs[cell] += value;

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
