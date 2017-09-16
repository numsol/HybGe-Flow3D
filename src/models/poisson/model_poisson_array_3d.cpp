/* poisson array 3d source */

// hgf includes
#include "model_poisson.hpp"

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

// simple distance formula
#define distance(x1,y1,z1,x2,y2,z2) \
  sqrt(pow((x1-x2),2) + pow((y1-y2),2) + pow((z1-z2),2))

void
hgf::models::poisson::build_array_3d(const parameters& par, const hgf::mesh::voxel& msh)
{
  // define temp coo arrays to store results in parallel region
  std::vector< std::vector< array_coo > > temp_arrays;
  temp_arrays.resize(NTHREADS);
  int maxp = block_size * 7;
  for (int ii = 0; ii < NTHREADS; ii++) temp_arrays[ii].reserve(maxp);

  bool alpha_diag;

#pragma omp parallel
  {
#pragma omp for schedule(dynamic) nowait num_threads(NTHREADS)
    for (int kk = 0; kk < NTHREADS; kk++) { 
      int entries = 0;
      array_coo temp_coo[9] = { 0 };
      double d_dofs[6], d_faces[6];
      int nbrs[6];
      double length, height, width;

      for (int ii = kk*block_size; ii < std::min((kk + 1)*block_size, (int)phi.size()); ii++) {
        entries = 0;
        for (int jj = 0; jj < 6; jj++) nbrs[jj] = phi[ii].neighbors[jj];

        alpha_diag = (alpha[ii][1] == 0.0 && \
                      alpha[ii][2] == 0.0 && \
                      alpha[ii][3] == 0.0 && \
                      alpha[ii][5] == 0.0 && \
                      alpha[ii][6] == 0.0 && \
                      alpha[ii][7] == 0.0 );

        for (int jj = 0; jj < 6; jj++) {
          if (nbrs[jj] > -1) {
            d_dofs[jj] = distance(phi[ii].coords[0], phi[ii].coords[1], phi[ii].coords[2], \
              phi[nbrs[jj]].coords[0], phi[nbrs[jj]].coords[1], phi[nbrs[jj]].coords[2]);
          }
        }

        // face dimensions
        // temporarily uniform...
        d_faces[0] = distance(msh.els[ii].vtx[0].coords[0], msh.els[ii].vtx[0].coords[1], msh.els[ii].vtx[0].coords[2], \
                              msh.els[ii].vtx[1].coords[0], msh.els[ii].vtx[1].coords[1], msh.els[ii].vtx[1].coords[2]) \
                   * distance(msh.els[ii].vtx[0].coords[0], msh.els[ii].vtx[0].coords[1], msh.els[ii].vtx[0].coords[2], \
                              msh.els[ii].vtx[3].coords[0], msh.els[ii].vtx[3].coords[1], msh.els[ii].vtx[3].coords[2]);
        d_faces[1] = d_faces[0];
        d_faces[2] = d_faces[0];
        d_faces[3] = d_faces[0];
        d_faces[4] = d_faces[0];
        d_faces[5] = d_faces[0];

        if (alpha_diag) {
          double alpha_cst; 
          // off diagonal entries
          for (int jj = 0; jj < 6; jj++) {
            if (nbrs[jj] > -1) {
              entries++;
              alpha_cst = (jj == 0 || jj == 2) ? (mean_perm(alpha[ii][4], alpha[nbrs[jj]][4])) : 
                         ((jj == 1 || jj == 3) ? (mean_perm(alpha[ii][0], alpha[nbrs[jj]][0])) : 
                                                  mean_perm(alpha[ii][8], alpha[nbrs[jj]][8]));
              temp_coo[entries - 1].value = -alpha_cst * d_faces[jj] / d_dofs[jj];
              temp_coo[entries - 1].i_index = ii;
              temp_coo[entries - 1].j_index = nbrs[jj];
            }
          }

          // diagonal entry
          entries++;
          temp_coo[entries - 1].value = 0;
          for (int jj = 0; jj < (entries - 1); jj++) {
            temp_coo[entries - 1].value -= temp_coo[jj].value;
          }
          temp_coo[entries - 1].i_index = ii;
          temp_coo[entries - 1].j_index = ii;
      
        }
        else { // nondiag permeability TODO
          std::cout << "Non-diagonal alpha tensors are a work in progress, and are not yet supported. Exiting.\n";
          exit(0);
        }

        // place values into temporoary coo array
        for (int jj = 0; jj < entries; jj++) {
          temp_arrays[kk].push_back(temp_coo[jj]);
        }
      }
    }
  }
  // paste
  for (int ii = 0; ii < NTHREADS; ii++) {
    for (int jj = 0; jj < temp_arrays[ii].size(); jj++) {
      coo_array.push_back(temp_arrays[ii][jj]);
    }
  }
}
