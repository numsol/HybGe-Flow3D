/* utilities source */

// system includes
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/filesystem.hpp>

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

// 1d->3d index
#define idx3(i, j, k, ldi1, ldi2) (k + (ldi2 * (j + ldi1 * i)))

#include "hgflow.hpp"

/** \brief Uniformly refines a voxelated input
 *
 * @param[in] par - parameters file containing mesh information.
 * @param[in] refine_len - integer controlling extent of geometry refinement.
 */
void
hgf::mesh::refine_voxel_uniform(parameters& par, int refine_len)
{

  int nx_old = par.nx;
  int ny_old = par.ny;
  int nz_old = par.nz;
  int size_old = par.voxel_geometry.size();
  std::vector< unsigned long > voxel_geometry_old(par.voxel_geometry);
  if (nz_old) { // 3D

    par.voxel_geometry.resize( refine_len * refine_len * refine_len * size_old ); 

    par.nx *= refine_len;
    par.ny *= refine_len;
    par.nz *= refine_len;

    int xi_n, yi_n, zi_n;

    for (int zi = 0; zi < nz_old; zi++) {
      for (int yi = 0; yi < ny_old; yi++) {
        for (int xi = 0; xi < nx_old; xi++) {
          xi_n = xi * refine_len;
          yi_n = yi * refine_len;
          zi_n = zi * refine_len;
          for (int rliz = 0; rliz < refine_len; rliz++) {
            for (int rliy = 0; rliy < refine_len; rliy++) {
              for (int rlix = 0; rlix < refine_len; rlix++) {
                par.voxel_geometry[idx3((zi_n+rliz), (yi_n+rliy), (xi_n+rlix), par.ny, par.nx)] = voxel_geometry_old[idx3(zi, yi, xi, ny_old, nx_old)];
              }
            }
          }
        }
      }
    }

  }
  else { // 2D

    par.voxel_geometry.resize( refine_len * refine_len * size_old ); 

    par.nx *= refine_len;
    par.ny *= refine_len;

    int xi_n, yi_n;

    for (int yi = 0; yi < ny_old; yi++) {
      for (int xi = 0; xi < nx_old; xi++) {
        for (int rliy = 0; rliy < refine_len; rliy++) {
          for (int rlix = 0; rlix < refine_len; rlix++) {
            xi_n = xi * refine_len;
            yi_n = yi * refine_len; 
            par.voxel_geometry[idx2((yi_n+rliy), (xi_n+rlix), par.nx)] = voxel_geometry_old[idx2(yi, xi, nx_old)];
          }
        }
      }
    }

  }

}

/** \brief Function removes cells from a mesh that are boundaries in opposite directions. Returns number of cells removed.
 *
 * @param[in] par - parameters file containing mesh information.
 */
int 
hgf::mesh::geo_sanity(parameters& par)
{
  int totalChanged = 0;
  int nChanged;

  if (par.dimension == 3) goto sanityCheck3;
  else goto sanityCheck2;

sanityCheck3:
  {
    nChanged = 0;
    // sanity
    for (int zi = 0; zi < par.nz; zi++) {
      for (int yi = 0; yi < par.ny; yi++) {
        for (int xi = 0; xi < par.nx; xi++) {
          if (par.voxel_geometry[idx3(zi, yi, xi, par.ny, par.nx)] != 1) {
            // xi sanity
            if (xi == 0) {
              if (par.voxel_geometry[idx3(zi, yi, (xi + 1), par.ny, par.nx)] == 1) {
                par.voxel_geometry[idx3(zi, yi, xi, par.ny, par.nx)] = 1;
                nChanged++;
              }
            }
            else if (xi == par.nx - 1) {
              if (par.voxel_geometry[idx3(zi, yi, (xi - 1), par.ny, par.nx)] == 1) {
                par.voxel_geometry[idx3(zi, yi, xi, par.ny, par.nx)] = 1;
                nChanged++;
              }
            }
            else {
              if (par.voxel_geometry[idx3(zi, yi, (xi + 1), par.ny, par.nx)] == 1 \
                && par.voxel_geometry[idx3(zi, yi, (xi - 1), par.ny, par.nx)] == 1) {
                par.voxel_geometry[idx3(zi, yi, xi, par.ny, par.nx)] = 1;
                nChanged++;
              }
            }
            // yi sanity
            if (yi == 0) {
              if (par.voxel_geometry[idx3(zi, (yi + 1), xi, par.ny, par.nx)] == 1) {
                par.voxel_geometry[idx3(zi, yi, xi, par.ny, par.nx)] = 1;
                nChanged++;
              }
            }
            else if (yi == par.ny - 1) {
              if (par.voxel_geometry[idx3(zi, (yi - 1), xi, par.ny, par.nx)] == 1) {
                par.voxel_geometry[idx3(zi, yi, xi, par.ny, par.nx)] = 1;
                nChanged++;
              }
            }
            else {
              if (par.voxel_geometry[idx3(zi, (yi + 1), xi, par.ny, par.nx)] == 1 \
                && par.voxel_geometry[idx3(zi, (yi - 1), xi, par.ny, par.nx)] == 1) {
                par.voxel_geometry[idx3(zi, yi, xi, par.ny, par.nx)] = 1;
                nChanged++;
              }
            }
            // zi sanity
            if (zi == 0) {
              if (par.voxel_geometry[idx3((zi + 1), yi, xi, par.ny, par.nx)] == 1) {
                par.voxel_geometry[idx3(zi, yi, xi, par.ny, par.nx)] = 1;
                nChanged++;
              }
            }
            else if (zi == par.nz - 1) {
              if (par.voxel_geometry[idx3((zi - 1), yi, xi, par.ny, par.nx)] == 1) {
                par.voxel_geometry[idx3(zi, yi, xi, par.ny, par.nx)] = 1;
                nChanged++;
              }
            }
            else {
              if (par.voxel_geometry[idx3((zi + 1), yi, xi, par.ny, par.nx)] == 1 \
                && par.voxel_geometry[idx3((zi - 1), yi, xi, par.ny, par.nx)] == 1) {
                par.voxel_geometry[idx3(zi, yi, xi, par.ny, par.nx)] = 1;
                nChanged++;
              }
            }
          }
        }
      }
    }
    totalChanged += nChanged;
    if (nChanged != 0) goto sanityCheck3;
    else goto cleanup;
  }

sanityCheck2:
  {
    nChanged = 0;
    // sanity
    for (int yi = 0; yi < par.ny; yi++) {
      for (int xi = 0; xi < par.nx; xi++) {
        if (par.voxel_geometry[idx2(yi, xi, par.nx)] != 1) {
          // xi sanity
          if (xi == 0) {
            if (par.voxel_geometry[idx2(yi, (xi + 1), par.nx)] == 1) {
              par.voxel_geometry[idx2(yi, xi, par.nx)] = 1;
              nChanged++;
            }
          }
          else if (xi == par.nx - 1) {
            if (par.voxel_geometry[idx2(yi, (xi - 1), par.nx)] == 1) {
              par.voxel_geometry[idx2(yi, xi, par.nx)] = 1;
              nChanged++;
            }
          }
          else {
            if (par.voxel_geometry[idx2(yi, (xi + 1), par.nx)] == 1 \
              && par.voxel_geometry[idx2(yi, (xi - 1), par.nx)] == 1) {
              par.voxel_geometry[idx2(yi, xi, par.nx)] = 1;
              nChanged++;
            }
          }
          // yi sanity
          if (yi == 0) {
            if (par.voxel_geometry[idx2((yi + 1), xi, par.nx)] == 1) {
              par.voxel_geometry[idx2(yi, xi, par.nx)] = 1;
              nChanged++;
            }
          }
          else if (yi == par.ny - 1) {
            if (par.voxel_geometry[idx2((yi - 1), xi, par.nx)] == 1) {
              par.voxel_geometry[idx2(yi, xi, par.nx)] = 1;
              nChanged++;
            }
          }
          else {
            if (par.voxel_geometry[idx2((yi + 1), xi, par.nx)] == 1 \
              && par.voxel_geometry[idx2((yi - 1), xi, par.nx)] == 1) {
              par.voxel_geometry[idx2(yi, xi, par.nx)] = 1;
              nChanged++;
            }
          }
        }
      }
    }
    totalChanged += nChanged;
    if (nChanged != 0) goto sanityCheck2;
    else goto cleanup;
  }

cleanup:
  if (totalChanged) {
    std::cout << "\nWarning, input geometry was incompatible.\n";
    std::cout << totalChanged << " cells, representing ";
    if (par.dimension == 3) std::cout << (double)100 * totalChanged / (par.nx*par.ny*par.nz);
    else std::cout << (double)100 * totalChanged / (par.nx*par.ny);
    std::cout << "% of the input geometry, with boundaries on opposite faces \nwere found and removed from void space.\n";
  }
  return totalChanged;
}

/** \brief Remove dead pores from a complex geometry. Returns number of pores removed.
 *
 * @param[in,out] par - parameters struct containing geometry information.
 */
int
hgf::mesh::remove_dead_pores(parameters& par)
{
  std::vector<unsigned long> voxel_geometry_cpy(par.voxel_geometry);
  std::vector<unsigned long> search_queue;
  std::vector<int> current_component;
  int n_components = 0;
  int pores_removed = 0;
  int ii, jj, kk;
  if (par.dimension == 3) {
    for (int k = 0; k < par.nz; k++) {
      for (int j = 0; j < par.ny; j++) {
        for (int i = 0; i < par.nx; i++) {
          if (voxel_geometry_cpy[idx3(k, j, i, par.ny, par.nx)] == 0 || voxel_geometry_cpy[idx3(k, j, i, par.ny, par.nx)] == 2) {
            ii = i;
            jj = j;
            kk = k;
            n_components++;
            voxel_geometry_cpy[idx3(kk, jj, ii, par.ny, par.nx)] = n_components + 2;
            current_component.resize(3);
            current_component[0] = k;
            current_component[1] = j;
            current_component[2] = i;
          }
          else continue;
          search_3d:
          // look y-
          if (jj) {
            if (voxel_geometry_cpy[idx3(kk, (jj - 1), ii, par.ny, par.nx)] == 0 || voxel_geometry_cpy[idx3(kk, (jj - 1), ii, par.ny, par.nx)] == 2) {
              voxel_geometry_cpy[idx3(kk, (jj - 1), ii, par.ny, par.nx)] = n_components + 2;
              current_component.push_back(kk);
              current_component.push_back((jj - 1));
              current_component.push_back(ii);
              search_queue.push_back(kk);
              search_queue.push_back((jj - 1));
              search_queue.push_back(ii);
            }
          }
          // look x+
          if (ii < par.nx - 1) {
            if (voxel_geometry_cpy[idx3(kk, jj, (ii + 1), par.ny, par.nx)] == 0 || voxel_geometry_cpy[idx3(kk, jj, (ii + 1), par.ny, par.nx)] == 2) {
              voxel_geometry_cpy[idx3(kk, jj, (ii + 1), par.ny, par.nx)] = n_components + 2;
              current_component.push_back(kk);
              current_component.push_back(jj);
              current_component.push_back((ii + 1));
              search_queue.push_back(kk);
              search_queue.push_back(jj);
              search_queue.push_back((ii + 1));
            }
          }
          // look y+
          if (jj < par.ny - 1) {
            if (voxel_geometry_cpy[idx3(kk, (jj + 1), ii, par.ny, par.nx)] == 0 || voxel_geometry_cpy[idx3(kk, (jj + 1), ii, par.ny, par.nx)] == 2) {
              voxel_geometry_cpy[idx3(kk, (jj + 1), ii, par.ny, par.nx)] = n_components + 2;
              current_component.push_back(kk);
              current_component.push_back((jj + 1));
              current_component.push_back(ii);
              search_queue.push_back(kk);
              search_queue.push_back((jj + 1));
              search_queue.push_back(ii);
            }
          }
          // look x-
          if (ii) {
            if (voxel_geometry_cpy[idx3(kk, jj, (ii - 1), par.ny, par.nx)] == 0 || voxel_geometry_cpy[idx3(kk, jj, (ii - 1), par.ny, par.nx)] == 2) {
              voxel_geometry_cpy[idx3(kk, jj, (ii - 1), par.ny, par.nx)] = n_components + 2;
              current_component.push_back(kk);
              current_component.push_back(jj);
              current_component.push_back((ii - 1));
              search_queue.push_back(kk);
              search_queue.push_back(jj);
              search_queue.push_back((ii - 1));
            }
          }
         
          // look z-
          if (kk) {
            if (voxel_geometry_cpy[idx3((kk - 1), jj, ii, par.ny, par.nx)] == 0 || voxel_geometry_cpy[idx3((kk - 1), jj, ii, par.ny, par.nx)] == 2) {
              voxel_geometry_cpy[idx3((kk - 1), jj, ii, par.ny, par.nx)] = n_components + 2;
              current_component.push_back((kk - 1));
              current_component.push_back(jj);
              current_component.push_back(ii);
              search_queue.push_back((kk - 1));
              search_queue.push_back(jj);
              search_queue.push_back(ii);
            }
          }
          // look z+
          if (kk < par.nz - 1) {
            if (voxel_geometry_cpy[idx3((kk + 1), jj, ii, par.ny, par.nx)] == 0 || voxel_geometry_cpy[idx3((kk + 1), jj, ii, par.ny, par.nx)] == 2) {
              voxel_geometry_cpy[idx3((kk + 1), jj, ii, par.ny, par.nx)] = n_components + 2;
              current_component.push_back((kk + 1));
              current_component.push_back(jj);
              current_component.push_back(ii);
              search_queue.push_back((kk + 1));
              search_queue.push_back(jj);
              search_queue.push_back(ii);
            }
          }
          if (search_queue.size()) {
            ii = search_queue[search_queue.size() - 1];
            jj = search_queue[search_queue.size() - 2];
            kk = search_queue[search_queue.size() - 3];
            search_queue.resize(search_queue.size() - 3);
            goto search_3d;
          }
          else { // done filling in this component, check if pore needs removing
            int min_i = par.nx;
            int min_j = par.ny;
            int min_k = par.nz;
            int max_i = 0;
            int max_j = 0;
            int max_k = 0;
            for (int cc = 0; cc < ((int)current_component.size() / 3); cc++) {
              if (current_component[idx2(cc, 0, 3)] > max_k) max_k = current_component[idx2(cc, 0, 3)];
              if (current_component[idx2(cc, 0, 3)] < min_k) min_k = current_component[idx2(cc, 0, 3)];
              if (current_component[idx2(cc, 1, 3)] > max_j) max_j = current_component[idx2(cc, 1, 3)];
              if (current_component[idx2(cc, 1, 3)] < min_j) min_j = current_component[idx2(cc, 1, 3)];
              if (current_component[idx2(cc, 2, 3)] > max_i) max_i = current_component[idx2(cc, 2, 3)];
              if (current_component[idx2(cc, 2, 3)] < min_i) min_i = current_component[idx2(cc, 2, 3)];
            }
            if (min_i > 0 || min_j > 0 || max_i < par.nx - 1 || max_j < par.ny - 1 || min_k > 0 || max_k < par.nz - 1) {
              pores_removed++;
              for (int cc = 0; cc < ((int)current_component.size() / 3); cc++) {
                par.voxel_geometry[idx3(current_component[idx2(cc, 0, 3)], current_component[idx2(cc, 1, 3)], current_component[idx2(cc, 2, 3)], par.ny, par.nx)] = 1;
              }
            }
          }
        }
      }
    }
  }
  else {
    for (int j = 0; j < par.ny; j++) {
      for (int i = 0; i < par.nx; i++) {
        if (voxel_geometry_cpy[idx2(j, i, par.nx)] == 0 || voxel_geometry_cpy[idx2(j, i, par.nx)] == 2) {
          ii = i;
          jj = j;
          n_components++;
          voxel_geometry_cpy[idx2(jj, ii, par.nx)] = n_components + 2;
          current_component.resize(2);
          current_component[0] = j;
          current_component[1] = i;
        }
        else continue;
        search_2d:
        // look down
        if (jj) {
          if (voxel_geometry_cpy[idx2((jj - 1), ii, par.nx)] == 0 || voxel_geometry_cpy[idx2((jj - 1), ii, par.nx)] == 2) {
            voxel_geometry_cpy[idx2((jj - 1), ii, par.nx)] = n_components + 2;
            current_component.push_back((jj - 1));
            current_component.push_back(ii);
            search_queue.push_back((jj - 1));
            search_queue.push_back(ii);
          }
        }
        // look right
        if (ii < par.nx - 1) {
          if (voxel_geometry_cpy[idx2(jj, (ii + 1), par.nx)] == 0 || voxel_geometry_cpy[idx2(jj, (ii + 1), par.nx)] == 2) {
            voxel_geometry_cpy[idx2(jj, (ii + 1), par.nx)] = n_components + 2;
            current_component.push_back(jj);
            current_component.push_back((ii + 1));
            search_queue.push_back(jj);
            search_queue.push_back((ii + 1));
          }
        }
        // look up
        if (jj < par.ny - 1) {
          if (voxel_geometry_cpy[idx2((jj + 1), ii, par.nx)] == 0 || voxel_geometry_cpy[idx2((jj + 1), ii, par.nx)] == 2) {
            voxel_geometry_cpy[idx2((jj + 1), ii, par.nx)] = n_components + 2;
            current_component.push_back((jj + 1));
            current_component.push_back(ii);
            search_queue.push_back((jj + 1));
            search_queue.push_back(ii);
          }
        }
        // look left 
        if (ii) {
          if (voxel_geometry_cpy[idx2(jj, (ii - 1), par.nx)] == 0 || voxel_geometry_cpy[idx2(jj, (ii - 1), par.nx)] == 2) {
            voxel_geometry_cpy[idx2(jj, (ii - 1), par.nx)] = n_components + 2;
            current_component.push_back(jj);
            current_component.push_back((ii - 1));
            search_queue.push_back(jj);
            search_queue.push_back((ii - 1));
          }
        }
        if (search_queue.size()) {
          ii = search_queue[search_queue.size() - 1];
          jj = search_queue[search_queue.size() - 2];
          search_queue.resize(search_queue.size() - 2);
          goto search_2d;
        }
        else { // done filling in this component, check if pore needs removing
          int min_i = par.nx;
          int min_j = par.ny;
          int max_i = 0;
          int max_j = 0;
          for (int cc = 0; cc < ((int)current_component.size() / 2); cc++) {
            if (current_component[idx2(cc, 0, 2)] > max_j) max_j = current_component[idx2(cc, 0, 2)];
            if (current_component[idx2(cc, 0, 2)] < min_j) min_j = current_component[idx2(cc, 0, 2)];
            if (current_component[idx2(cc, 1, 2)] > max_i) max_i = current_component[idx2(cc, 1, 2)];
            if (current_component[idx2(cc, 1, 2)] < min_i) min_i = current_component[idx2(cc, 1, 2)];
          }
          if (min_i > 0 || min_j > 0 || max_i < par.nx - 1 || max_j < par.ny - 1) {
            pores_removed++;
            for (int cc = 0; cc < ((int)current_component.size() / 2); cc++) {
              par.voxel_geometry[idx2(current_component[idx2(cc, 0, 2)], current_component[idx2(cc, 1, 2)], par.nx)] = 1;
            }
          }
        }
      }
    }
  }

  return pores_removed;
}
