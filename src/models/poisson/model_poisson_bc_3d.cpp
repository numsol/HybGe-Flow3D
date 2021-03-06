/* poisson bc 2d source */

// hgf includes
#include "model_poisson.hpp"

// simple distance formula
#define distance(x1,y1,z1,x2,y2,z2) sqrt(pow((x1-x2),2) + pow((y1-y2),2) + pow((z1-z2),2))

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

void
hgf::models::poisson::homogeneous_dirichlet_3d(const parameters& par, const hgf::mesh::voxel& msh)
{
  
  // define temp coo arays to store results in parallel region
  std::vector< std::vector< array_coo > > temp_arrays;
  temp_arrays.resize(NTHREADS);

  int maxp = 2*block_size;
  for (int ii = 0; ii < NTHREADS; ii++) { 
    temp_arrays[ii].reserve(maxp); 
  }

  bool alpha_diag;

#pragma omp for private(alpha_diag) schedule(dynamic) 
  for (int kk = 0; kk < NTHREADS; kk++) {
      
    int nbrs[6];
    array_coo temp_coo;
    double dx, dy, dz;

    for (int ii = kk*block_size; ii < std::min((kk + 1)*block_size, (int)phi.size()); ii++) {
      double value = 0;
      int bc_contributor[6] = { 0, 0, 0, 0, 0, 0 };
      int nnbr = 0;
      for (int jj = 0; jj < 6; jj++) {
        nbrs[jj] = phi[ii].neighbors[jj];
        if (nbrs[jj] != -1) nnbr++;
        else {
          bc_types[ii][jj] = 1;
          bc_contributor[jj] = 1;
        }
      }
      if (nnbr == 6) goto phiexit;

      alpha_diag = (alpha[ii][1] == 0.0 && \
                    alpha[ii][2] == 0.0 && \
                    alpha[ii][3] == 0.0 && \
                    alpha[ii][5] == 0.0 && \
                    alpha[ii][6] == 0.0 && \
                    alpha[ii][7] == 0.0);

      dx = msh.els[ii].vtx[1].coords[0] - msh.els[ii].vtx[0].coords[0];
      dy = msh.els[ii].vtx[2].coords[1] - msh.els[ii].vtx[1].coords[1];
      dz = msh.els[ii].vtx[7].coords[2] - msh.els[ii].vtx[0].coords[2];

      if (alpha_diag) {

        // y- neighbor?
        if (bc_contributor[0]) {
          value += alpha[ii][4] * dx * dz / (0.5 * dy);
        }

        // x+ neighbor?
        if (bc_contributor[1]) {
          value += alpha[ii][0] * dz * dy / (0.5 * dx);
        }

        // y+ neighbor?
        if (bc_contributor[2]) {
          value += alpha[ii][4] * dx * dz / (0.5 * dy);
        }

        // x- neighbor?
        if (bc_contributor[3]) {
          value += alpha[ii][0] * dy * dz / (0.5 * dx);
        }

        // z- neighbor?
        if (bc_contributor[4]) {
          value += alpha[ii][8] * dx * dy / (0.5 * dz);
        }

        // z+ neighbor?
        if (bc_contributor[5]) {
          value += alpha[ii][8] * dx * dy / (0.5 * dz);
        }

      }
      else { // Nondiag alpha, TODO

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
hgf::models::poisson::homogeneous_mixed_3d(const parameters& par, const hgf::mesh::voxel& msh, \
                                           bool (*is_dirichlet)( const parameters& par, int dof_num, double coords[3] ))
{
  // define temp coo arays to store results in parallel region
  std::vector< std::vector< array_coo > > temp_arrays;
  temp_arrays.resize(NTHREADS);
  int maxp = 2*block_size;
  for (int ii = 0; ii < NTHREADS; ii++) temp_arrays[ii].reserve(maxp); 

  bool alpha_diag;

#pragma omp for private(alpha_diag) schedule(dynamic) 
  for (int kk = 0; kk < NTHREADS; kk++) {
      
    int nbrs[6];
    array_coo temp_coo;
    double dx, dy, dz;

    for (int ii = kk*block_size; ii < std::min((kk + 1)*block_size, (int)phi.size()); ii++) {
      double coords[3];
      double value = 0;
      int bc_contributor[6] = { 0, 0, 0, 0, 0, 0 };
      int nnbr = 0;
      for (int jj = 0; jj < 6; jj++) {
        nbrs[jj] = phi[ii].neighbors[jj];
        if (nbrs[jj] != -1) nnbr++;
        else bc_contributor[jj] = 1;
      }
      if (nnbr == 6) goto phiexit;

      alpha_diag = (alpha[ii][1] == 0.0 && \
                    alpha[ii][2] == 0.0 && \
                    alpha[ii][3] == 0.0 && \
                    alpha[ii][5] == 0.0 && \
                    alpha[ii][6] == 0.0 && \
                    alpha[ii][7] == 0.0);

      dx = msh.els[ii].vtx[1].coords[0] - msh.els[ii].vtx[0].coords[0];
      dy = msh.els[ii].vtx[2].coords[1] - msh.els[ii].vtx[1].coords[1];
      dz = msh.els[ii].vtx[7].coords[2] - msh.els[ii].vtx[0].coords[2];

      if (alpha_diag) {

        // y- neighbor?
        if (bc_contributor[0]) {
          coords[0] = phi[ii].coords[0];
          coords[1] = phi[ii].coords[1] - 0.5 * dy;
          coords[2] = phi[ii].coords[2];
          if (is_dirichlet( par, ii, coords )) {
            bc_types[ii][0] = 1;
            value += alpha[ii][4] * dx * dz / (0.5 * dy);
          }
          else bc_types[ii][0] = 2;
        }

        // x+ neighbor?
        if (bc_contributor[1]) {
          coords[0] = phi[ii].coords[0] + 0.5 * dx;
          coords[1] = phi[ii].coords[1];
          coords[2] = phi[ii].coords[2];
          if (is_dirichlet( par, ii, coords )) {
            bc_types[ii][1] = 1;
            value += alpha[ii][0] * dz * dy / (0.5 * dx);
          }
          else bc_types[ii][1] = 2;
        }

        // y+ neighbor?
        if (bc_contributor[2]) {
          coords[0] = phi[ii].coords[0];
          coords[1] = phi[ii].coords[1] + 0.5 * dy;
          coords[2] = phi[ii].coords[2];
          if (is_dirichlet( par, ii, coords )) {
            bc_types[ii][2] = 1;
            value += alpha[ii][4] * dx * dz / (0.5 * dy);
          }
          else bc_types[ii][2] = 2;
        }

        // x- neighbor?
        if (bc_contributor[3]) {
          coords[0] = phi[ii].coords[0] - 0.5 * dx;
          coords[1] = phi[ii].coords[1];
          coords[2] = phi[ii].coords[2];
          if (is_dirichlet( par, ii, coords )) {
            bc_types[ii][3] = 1;
            value += alpha[ii][0] * dy * dz / (0.5 * dx);
          }
          else bc_types[ii][3] = 2;
        }

        // z- neighbor?
        if (bc_contributor[4]) {
          coords[0] = phi[ii].coords[0];
          coords[1] = phi[ii].coords[1];
          coords[2] = phi[ii].coords[2] - 0.5 * dz;
          if (is_dirichlet( par, ii, coords )) {
            bc_types[ii][4] = 1;
            value += alpha[ii][8] * dx * dy / (0.5 * dz);
          }
          else bc_types[ii][4] = 2;
        }

        // z+ neighbor?
        if (bc_contributor[5]) {
          coords[0] = phi[ii].coords[0];
          coords[1] = phi[ii].coords[1];
          coords[2] = phi[ii].coords[2] + 0.5 * dz;
          if (is_dirichlet( par, ii, coords )) {
            bc_types[ii][5] = 1;
            value += alpha[ii][8] * dx * dy / (0.5 * dz);
          }
          else bc_types[ii][5] = 2;
        }

      }
      else { // Nondiag alpha, TODO

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
