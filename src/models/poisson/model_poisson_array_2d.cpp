/* poisson array 2d source */

// hgf includes
#include "model_poisson.hpp"

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

// simple distance formula
#define distance(x1,y1,x2,y2) \
  sqrt(pow((x1-x2),2) + pow((y1-y2),2))

void
hgf::models::poisson::build_array_2d(const parameters& par, const hgf::mesh::voxel& msh)
{
  // define temp coo arrays to store results in parallel region
  std::vector< std::vector< array_coo > > temp_arrays;
  temp_arrays.resize(NTHREADS);
  int maxp = block_size * 5;
  for (int ii = 0; ii < NTHREADS; ii++) temp_arrays[ii].reserve(maxp);

  bool alpha_diag;

#pragma omp parallel
  {
#pragma omp for private(alpha_diag) schedule(dynamic)
    for (int kk = 0; kk < NTHREADS; kk++) {
      int entries = 0;
      array_coo temp_coo[5] = { 0 };
      double d_dofs[4], d_edges[4];
      int nbrs[4];
      for (int ii = kk*block_size; ii < std::min((kk + 1)*block_size, (int)phi.size()); ii++) {
        entries = 0;
        // grab neighbor numbers
        for (int jj = 0; jj < 4; jj++) { nbrs[jj] = phi[ii].neighbors[jj]; }

        alpha_diag = (alpha[ii][1] == 0.0 && alpha[ii][2] == 0.0);
        for (int jj = 0; jj < 4; jj++) {
          if (nbrs[jj] > -1) {
            if (alpha[nbrs[jj]][1] != 0.0 || alpha[nbrs[jj]][2] != 0.0) {
              alpha_diag = 0;
            }
          }
        }

        // compute cell center distances
        for (int jj = 0; jj < 4; jj++) {
          if (nbrs[jj] > -1) {
            d_dofs[jj] = distance(phi[ii].coords[0], phi[ii].coords[1], \
              phi[nbrs[jj]].coords[0], phi[nbrs[jj]].coords[1]);
          }
        }
        // compute edge distances
        for (int jj = 0; jj < 4; jj++) {
          int nn = (jj < 3) ? (jj + 1) : 0;
          d_edges[jj] = distance(msh.els[ii].vtx[jj].coords[0], msh.els[ii].vtx[jj].coords[1], \
                                 msh.els[ii].vtx[nn].coords[0], msh.els[ii].vtx[nn].coords[1]);
        }

        if (alpha_diag) {
          double alpha_cst; 
          // off diagonal entries
          for (int jj = 0; jj < 4; jj++) {
            if (nbrs[jj] > -1) {
              entries++;
              alpha_cst = (jj == 0 || jj == 2) ? (mean_perm(alpha[ii][3], alpha[nbrs[jj]][3])) : (mean_perm(alpha[ii][0], alpha[nbrs[jj]][0]));
              temp_coo[entries - 1].value = -alpha_cst * d_edges[jj] / d_dofs[jj]; // no alpha yet
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

        } else { // non-diagonal alpha, TODO
          for (int jj = 0; jj < 4; jj++) {
            if (nbrs[jj] > -1) {
              double alpha_tns[4] = { mean_perm(alpha[ii][0], alpha[nbrs[jj]][0]), \
                                      mean_perm(alpha[ii][1], alpha[nbrs[jj]][1]), \
                                      mean_perm(alpha[ii][2], alpha[nbrs[jj]][2]), \
                                      mean_perm(alpha[ii][3], alpha[nbrs[jj]][3]) };
              entries++;
              //--- WIP ---//
              std::cout << "Non-diagonal alpha tensors are a work in progress, and are not yet supported. Exiting.\n";
              exit(0);
            }
          }
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

