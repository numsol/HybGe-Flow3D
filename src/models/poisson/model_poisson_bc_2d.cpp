/* poisson bc 2d source */

// hgf includes
#include "model_poisson.hpp"

// simple distance formula
#define distance(x1,y1,x2,y2) \
  sqrt(pow((x1-x2),2) + pow((y1-y2),2))

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

void
hgf::models::poisson::homogeneous_dirichlet_2d(const parameters& par, const hgf::mesh::voxel& msh)
{
  // define temp coo array to store results in parallel region
  std::vector< std::vector< array_coo > > temp_arrays;
  temp_arrays.resize(NTHREADS);
  int maxp = 2*block_size;
  for (int ii = 0; ii < NTHREADS; ii++) temp_arrays[ii].reserve(maxp);

  bool alpha_diag;

#pragma omp parallel for private(alpha_diag) schedule(dynamic) num_threads(NTHREADS)
  for (int kk = 0; kk < NTHREADS; kk++) {

    int nbrs[4];
    array_coo temp_coo;
    double dx, dy;

    for (int ii = kk*block_size; ii < std::min((kk + 1)*block_size, (int)phi.size()); ii++) {
      double value = 0;
      int bc_contributor[4] = { 0, 0, 0, 0 };
      int nnbr = 0;
      for (int jj = 0; jj < 4; jj++) {
        nbrs[jj] = phi[ii].neighbors[jj];
        if (nbrs[jj] != -1) nnbr++;
        else {
          bc_types[ii][jj] = 1;
          bc_contributor[jj] = 1;
        }
      }
      if (nnbr == 4) goto phiexit;

      alpha_diag = (alpha[ii][1] == 0.0 && alpha[ii][2] == 0.0);

      dx = msh.els[ii].vtx[1].coords[0] - msh.els[ii].vtx[0].coords[0];
      dy = msh.els[ii].vtx[2].coords[1] - msh.els[ii].vtx[1].coords[1]; 

      if (alpha_diag) {
      
        // S neighbor?
        if (bc_contributor[0]) {
          value += alpha[ii][3] * dx / (0.5*dy);
        } 

        // E neighbor?
        if (bc_contributor[1]) {
          value += alpha[ii][0] * dy / (0.5*dx); 
        }

        // N neighbor?
        if (bc_contributor[2]) {
          value += alpha[ii][3] * dx / (0.5*dy); 
        }

        // W neighbor?
        if (bc_contributor[3]) {
          value += alpha[ii][0] * dy / (0.5*dx);
        }

      }
      else { // Nondiag or non constant alpha, TODO

      }

      temp_coo.i_index = ii;
      temp_coo.j_index = ii;
      temp_coo.value = value;

      temp_arrays[kk].push_back(temp_coo);

    phiexit:;
    }
  }
  // paste
  for (int ii = 0; ii < NTHREADS; ii++) {
    for (int jj = 0; jj < temp_arrays[ii].size(); jj++) {
      coo_array.push_back(temp_arrays[ii][jj]);
    }
  }
}

void 
hgf::models::poisson::homogeneous_mixed_2d(const parameters& par, const hgf::mesh::voxel& msh, \
                                           bool (*is_dirichlet)( const parameters& par, int dof_num, double coords[3] ))
{
  // define temp coo array to store results in parallel region
  std::vector< std::vector< array_coo > > temp_arrays;
  temp_arrays.resize(NTHREADS);
  int maxp = 2*block_size;
  for (int ii = 0; ii < NTHREADS; ii++) temp_arrays[ii].reserve(maxp);

  bool alpha_diag;

#pragma omp parallel for private(alpha_diag) schedule(dynamic) num_threads(NTHREADS)
  for (int kk = 0; kk < NTHREADS; kk++) {

    int nbrs[4];
    array_coo temp_coo;
    double dx, dy;

    for (int ii = kk*block_size; ii < std::min((kk + 1)*block_size, (int)phi.size()); ii++) {
      double coords[3];
      double value = 0;
      int bc_contributor[4] = { 0, 0, 0, 0 };
      int nnbr = 0;
      for (int jj = 0; jj < 4; jj++) {
        nbrs[jj] = phi[ii].neighbors[jj];
        if (nbrs[jj] != -1) nnbr++;
        else bc_contributor[jj] = 1;
      }
      if (nnbr == 4) goto phiexit;

      alpha_diag = (alpha[ii][1] == 0.0 && alpha[ii][2] == 0.0);

      dx = msh.els[ii].vtx[1].coords[0] - msh.els[ii].vtx[0].coords[0];
      dy = msh.els[ii].vtx[2].coords[1] - msh.els[ii].vtx[1].coords[1]; 

      if (alpha_diag) {
      
        // S neighbor?
        if (bc_contributor[0]) {
          coords[0] = phi[ii].coords[0];
          coords[1] = phi[ii].coords[1] - 0.5*dy;
          if (is_dirichlet( par, ii, coords )) {
            bc_types[ii][0] = 1;
            value += alpha[ii][3] * dx / (0.5*dy);
          }
          else bc_types[ii][0] = 2;
        } 

        // E neighbor?
        if (bc_contributor[1]) {
          coords[0] = phi[ii].coords[0] + 0.5*dx;
          coords[1] = phi[ii].coords[1]; 
          if (is_dirichlet( par, ii, coords )) {
            bc_types[ii][1] = 1;
            value += alpha[ii][0] * dy / (0.5*dx); 
          }
          else bc_types[ii][1] = 2;
        }

        // N neighbor?
        if (bc_contributor[2]) {
          coords[0] = phi[ii].coords[0];
          coords[1] = phi[ii].coords[1] + 0.5*dy;
          if (is_dirichlet( par, ii, coords )) {
            bc_types[ii][2] = 1;
            value += alpha[ii][3] * dx / (0.5*dy); 
          }
          else bc_types[ii][2] = 2;
        }

        // W neighbor?
        if (bc_contributor[3]) {
          coords[0] = phi[ii].coords[0] - 0.5*dx;
          coords[1] = phi[ii].coords[1];
          if (is_dirichlet( par, ii, coords )) {
            bc_types[ii][3] = 1;
            value += alpha[ii][0] * dy / (0.5*dx);
          }
          else bc_types[ii][3] = 2;
        }

      }
      else { // Nondiag alpha, TODO

      }

      temp_coo.i_index = ii;
      temp_coo.j_index = ii;
      temp_coo.value = value;

      if (temp_coo.value) temp_arrays[kk].push_back(temp_coo);

    phiexit:;
    }
  }
  // paste
  for (int ii = 0; ii < NTHREADS; ii++) {
    for (int jj = 0; jj < temp_arrays[ii].size(); jj++) {
      coo_array.push_back(temp_arrays[ii][jj]);
    }
  }

}
