#include "multiscale_flow.hpp"

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

static inline void
solve_4x4( double *mat, double *v )
{

  double temp;

  // LU (no pivot)
  for (int j = 0; j < 4; j++) {
    for (int i = j + 1; i < 4; i++) {
      mat[i + 4*j] /= mat[j + 4*j];
    }
    for (int jj = j + 1; jj < 4; jj++) {
      temp = mat[j + 4*jj];
      for (int ii = 0; ii < 4; ii++) {
        mat[ii + jj*4] -= mat[ii + j*4] * temp;
      }
    }
  }

  // Form U * x = L^-1 * v, TRSV(LNU)
  for (int j = 0; j < 4; j++) {
    temp = v[j];
    for (int i = j + 1; i < 4; i++) {
      v[i] -= temp * mat[i + j*4];
    }
  }

  // Form x = U^-1 L^-1 * v, TRSV(UNN)
  for (int j = 3; j >= 0; j--) {
    v[j] /= mat[j + 4*j];
    temp = v[j];
    for (int i = j - 1; i >= 0; i--) {
      v[i] -= temp * mat[i + j*4];
    }
  }

}

static inline void
solve_9x9( double *mat, double *v )
{
 
  double temp;

  // LU (no pivot)
  for (int j = 0; j < 9; j++) {
    for (int i = j + 1; i < 9; i++) {
      mat[i + 9*j] /= mat[j + 9*j];
    }
    for (int jj = j + 1; jj < 9; jj++) {
      temp = mat[j + 9*jj];
      for (int ii = 0; ii < 9; ii++) {
        mat[ii + jj*9] -= mat[ii + j*9] * temp;
      }
    }
  }

  // Form U * x = L^-1 * v, TRSV(LNU)
  for (int j = 0; j < 9; j++) {
    temp = v[j];
    for (int i = j + 1; i < 9; i++) {
      v[i] -= temp * mat[i + j*9];
    }
  }

  // Form x = U^-1 L^-1 * v, TRSV(UNN)
  for (int j = 8; j >= 0; j--) {
    v[j] /= mat[j + 9*j];
    temp = v[j];
    for (int i = j - 1; i >= 0; i--) {
      v[i] -= temp * mat[i + j*9];
    }
  }
}

void
compute_averages_x(const parameters& par, const std::vector< int >& pressure_ib_list, \
                                          const std::vector< degree_of_freedom >& velocity_u, \
                                          const std::vector< degree_of_freedom >& velocity_v, \
                                          const std::vector< degree_of_freedom >& velocity_w, \
                                          const std::vector< double > solution, double& v, double& g)
{
  double min_x, max_x, mid_x, midrange_x, min_y, max_y, min_z, max_z, pressure;
  double p1 = 0;
  double p2 = 0;
  double r = 0.09;
  int p1_idx, p2_idx;
  int v_idx = 0;
  int p1_count = 0;
  int p2_count = 0;
  int v_count = 0;
  int n_velocity;

  g = 0;
  v = 0;

  // compute v, g 3d
  if (par.dimension == 3) {
    n_velocity = (int)velocity_u.size() + (int)velocity_v.size() + (int)velocity_w.size();
    // geometry properties
    min_x = 0;
    max_x = par.length;
    min_y = 0;
    max_y = par.width;
    min_z = 0;
    max_z = par.height;
    // adjust limits to avoid recirculation
    min_x += r*par.length;
    max_x -= r*par.length;
    mid_x = 0.5 * (min_x + max_x);
    min_y += r*par.width;
    max_y -= r*par.width;
    min_z += r*par.height;
    max_z -= r*par.height;
    midrange_x = 0.5 * (max_x + mid_x) - 0.5 * (mid_x + min_x);
    // compute averages
    for (int ii = 0; ii < velocity_u.size(); ii++) {
      p1_idx = velocity_u[ii].cell_numbers[0];
      p2_idx = velocity_u[ii].cell_numbers[1];
      if (p1_idx != -1 && p2_idx != -1 && !pressure_ib_list[p1_idx] && !pressure_ib_list[p2_idx]) {
        if (velocity_u[ii].coords[0] > min_x) {
          if (velocity_u[ii].coords[0] < max_x) {
            if (velocity_u[ii].coords[1] > min_y) {
              if (velocity_u[ii].coords[1] < max_y) {
                if (velocity_u[ii].coords[2] > min_z) {
                  if (velocity_u[ii].coords[2] < max_z) {
                    pressure = 0.5 * (solution[n_velocity + p1_idx] + solution[n_velocity + p2_idx]);
                    if (velocity_u[ii].coords[0] < mid_x) {
                      p1 += pressure;
                      p1_count++;
                      v += solution[ii];
                      v_count++;
                    }
                    else if (velocity_u[ii].coords[0] >= mid_x) {
                      p2 += pressure;
                      p2_count++;
                      v += solution[ii];
                      v_count++;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  // compute v, g 2d
  else {
    n_velocity = (int)velocity_u.size() + (int)velocity_v.size();
    // geometry properties
    min_x = 0;
    max_x = par.length;
    min_y = 0;
    max_y = par.width;
    // adjust limits to avoid recirculation
    min_x += r*par.length;
    max_x -= r*par.length;
    mid_x = 0.5 * (min_x + max_x);
    min_y += r*par.width;
    max_y -= r*par.width;
    midrange_x = 0.5 * (max_x + mid_x) - 0.5 * (mid_x + min_x);
    // compute averages
    for (int ii = 0; ii < velocity_u.size(); ii++) {
      p1_idx = velocity_u[ii].cell_numbers[0];
      p2_idx = velocity_u[ii].cell_numbers[1];
      if (p1_idx != -1 && p2_idx != -1 && !pressure_ib_list[p1_idx] && !pressure_ib_list[p2_idx]) {
        if (velocity_u[ii].coords[0] > min_x) {
          if (velocity_u[ii].coords[0] < max_x) {
            if (velocity_u[ii].coords[1] > min_y) {
              if (velocity_u[ii].coords[1] < max_y) {
                pressure = 0.5 * (solution[n_velocity + p1_idx] + solution[n_velocity + p2_idx]);
                if (velocity_u[ii].coords[0] < mid_x) {
                  p1 += pressure;
                  p1_count++;
                  v += solution[ii];
                  v_count++;
                }
                else if (velocity_u[ii].coords[0] >= mid_x) {
                  p2 += pressure;
                  p2_count++;
                  v += solution[ii];
                  v_count++;
                }
              }
            }
          }
        }
      }
    }
  }
  v /= v_count;
  p1 /= p1_count;
  p2 /= p2_count;
  g = (p1 - p2) / midrange_x;
}

void
compute_averages_y(const parameters& par, const std::vector< int >& pressure_ib_list, \
                                          const std::vector< degree_of_freedom >& velocity_u, \
                                          const std::vector< degree_of_freedom >& velocity_v, \
                                          const std::vector< degree_of_freedom >& velocity_w, \
                                          const std::vector< double > solution, double& v, double& g)
{
  double min_x, max_x, min_y, max_y, mid_y, midrange_y, min_z, max_z, pressure;
  double p1 = 0;
  double p2 = 0;
  double r = 0.09;
  int p1_idx, p2_idx;
  int v_idx = 0;
  int p1_count = 0;
  int p2_count = 0;
  int v_count = 0;
  int n_velocity;

  g = 0;
  v = 0;

  // compute v, g 3d
  if (par.dimension == 3) {
    n_velocity = (int)velocity_u.size() + (int)velocity_v.size() + (int)velocity_w.size();
    // geometry properties
    min_x = 0;
    max_x = par.length;
    min_y = 0;
    max_y = par.width;
    min_z = 0;
    max_z = par.height;
    // adjust limits to avoid recirculation
    min_x += r*par.length;
    max_x -= r*par.length;
    min_y += r*par.width;
    max_y -= r*par.width;
    mid_y = 0.5 * (min_y + max_y);
    min_z += r*par.height;
    max_z -= r*par.height;
    midrange_y = 0.5 * (max_y + mid_y) - 0.5 * (mid_y + min_y);
    // compute averages
    for (int ii = 0; ii < velocity_v.size(); ii++) {
      p1_idx = velocity_v[ii].cell_numbers[0];
      p2_idx = velocity_v[ii].cell_numbers[1];
      if (p1_idx != -1 && p2_idx != -1 && !pressure_ib_list[p1_idx] && !pressure_ib_list[p2_idx]) {
        if (velocity_v[ii].coords[0] > min_x) {
          if (velocity_v[ii].coords[0] < max_x) {
            if (velocity_v[ii].coords[1] > min_y) {
              if (velocity_v[ii].coords[1] < max_y) {
                if (velocity_v[ii].coords[2] > min_z) {
                  if (velocity_v[ii].coords[2] < max_z) {
                    pressure = 0.5 * (solution[n_velocity + p1_idx] + solution[n_velocity + p2_idx]);
                    if (velocity_v[ii].coords[1] < mid_y) {
                      p1 += pressure;
                      p1_count++;
                      v += solution[velocity_u.size() + ii];
                      v_count++;
                    }
                    else if (velocity_v[ii].coords[1] >= mid_y) {
                      p2 += pressure;
                      p2_count++;
                      v += solution[velocity_u.size() + ii];
                      v_count++;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  // compute v, g 2d
  else {
    n_velocity = (int)velocity_u.size() + (int)velocity_v.size();
    // geometry properties
    min_x = 0;
    max_x = par.length;
    min_y = 0;
    max_y = par.width;
    // adjust limits to avoid recirculation
    min_x += r*par.length;
    max_x -= r*par.length;
    min_y += r*par.width;
    max_y -= r*par.width;
    mid_y = 0.5 * (min_y + max_y);
    midrange_y = 0.5 * (max_y + mid_y) - 0.5 * (mid_y + min_y);
    // compute averages
    for (int ii = 0; ii < velocity_v.size(); ii++) {
      p1_idx = velocity_v[ii].cell_numbers[0];
      p2_idx = velocity_v[ii].cell_numbers[1];
      if (p1_idx != -1 && p2_idx != -1 && !pressure_ib_list[p1_idx] && !pressure_ib_list[p2_idx]) {
        if (velocity_v[ii].coords[0] > min_x) {
          if (velocity_v[ii].coords[0] < max_x) {
            if (velocity_v[ii].coords[1] > min_y) {
              if (velocity_v[ii].coords[1] < max_y) {
                pressure = 0.5 * (solution[n_velocity + p1_idx] + solution[n_velocity + p2_idx]);
                if (velocity_v[ii].coords[1] < mid_y) {
                  p1 += pressure;
                  p1_count++;
                  v += solution[velocity_u.size() + ii];
                  v_count++;
                }
                else if (velocity_v[ii].coords[1] >= mid_y) {
                  p2 += pressure;
                  p2_count++;
                  v += solution[velocity_u.size() + ii];
                  v_count++;
                }
              }
            }
          }
        }
      }
    }
  }
  v /= v_count;
  p1 /= p1_count;
  p2 /= p2_count;
  g = (p1 - p2) / midrange_y;
}

void
compute_averages_z(const parameters& par, const std::vector< int >& pressure_ib_list, \
                                          const std::vector< degree_of_freedom >& velocity_u, \
                                          const std::vector< degree_of_freedom >& velocity_v, \
                                          const std::vector< degree_of_freedom >& velocity_w, \
                                          const std::vector< double > solution, double& v, double& g)
{
  double min_x, max_x, min_y, max_y, min_z, max_z, mid_z, midrange_z, pressure;
  double p1 = 0;
  double p2 = 0;
  double r = 0.09;
  int p1_idx, p2_idx;
  int v_idx = 0;
  int p1_count = 0;
  int p2_count = 0;
  int v_count = 0;
  int n_velocity;

  g = 0;
  v = 0;

  // compute v, g 3d
  n_velocity = (int)velocity_u.size() + (int)velocity_v.size() + (int)velocity_w.size();
  // geometry properties
  min_x = 0;
  max_x = par.length;
  min_y = 0;
  max_y = par.width;
  min_z = 0;
  max_z = par.height;
  // adjust limits to avoid recirculation
  min_x += r*par.length;
  max_x -= r*par.length;
  min_y += r*par.width;
  max_y -= r*par.width;
  min_z += r*par.height;
  max_z -= r*par.height;
  mid_z = 0.5 * (min_z + max_z);
  midrange_z = 0.5 * (max_z + mid_z) - 0.5 * (mid_z + min_z);
  // compute averages
  for (int ii = 0; ii < velocity_w.size(); ii++) {
    p1_idx = velocity_w[ii].cell_numbers[0];
    p2_idx = velocity_w[ii].cell_numbers[1];
    if (p1_idx != -1 && p2_idx != -1 && !pressure_ib_list[p1_idx] && !pressure_ib_list[p2_idx]) {
      if (velocity_w[ii].coords[0] > min_x) {
        if (velocity_w[ii].coords[0] < max_x) {
          if (velocity_w[ii].coords[1] > min_y) {
            if (velocity_w[ii].coords[1] < max_y) {
              if (velocity_w[ii].coords[2] > min_z) {
                if (velocity_w[ii].coords[2] < max_z) {
                  pressure = 0.5 * (solution[n_velocity + p1_idx] + solution[n_velocity + p2_idx]);
                  if (velocity_w[ii].coords[2] < mid_z) {
                    p1 += pressure;
                    p1_count++;
                    v += solution[velocity_u.size() + velocity_v.size() + ii];
                    v_count++;
                  }
                  else if (velocity_w[ii].coords[2] >= mid_z) {
                    p2 += pressure;
                    p2_count++;
                    v += solution[velocity_u.size() + velocity_v.size() + ii];
                    v_count++;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  v /= v_count;
  p1 /= p1_count;
  p2 /= p2_count;
  g = (p1 - p2) / midrange_z;
}


/** \brief hgf::multiscale::flow::compute_permeability_x computes upscaled permeability in the x direction from a porescale flow solution.
 *
 * @param[in] par - parameters struct containing problem information.
 * @param[in] pressure_ib_list - Integer array indicating if a cell is an immersed boundary cell.
 * @param[in] velocity_u - Degrees of freedom for the x-component of the porescale velocity solution.
 * @param[in] velocity_v - Degrees of freedom for the y-component of the porescale velocity solution.
 * @param[in] velocity_w - Degrees of freedom for the z-component of the porescale velocity solution.
 * @param[in] solution - Porescale flow solution.
 */
double
hgf::multiscale::flow::compute_permeability_x(const parameters& par, const std::vector< int >& pressure_ib_list, \
                                                                     const std::vector< degree_of_freedom >& velocity_u, \
                                                                     const std::vector< degree_of_freedom >& velocity_v, \
                                                                     const std::vector< degree_of_freedom >& velocity_w, \
                                                                     const std::vector< double > solution)
{
  // quick exit
  if (solution.size() <= 0) {
    std::cout << "\nEmpty flow solution, returning -1 permeability.\n";
    return -1.0;
  }

  double v, g, por;

  compute_averages_x(par, pressure_ib_list, velocity_u, velocity_v, velocity_w, solution, v, g);

  // determine porosity (holder)
  int n_voxels = 0;
  int n_void = 0;
  for (int ii = 0; ii < par.voxel_geometry.size(); ii++) {
    n_voxels++;
    if (par.voxel_geometry[ii] == 0) n_void++;
  }
  por = (double)n_void / n_voxels;

  return (v/g) * por;

}

/** \brief hgf::multiscale::flow::compute_permeability_y computes upscaled permeability in the y direction from a porescale flow solution.
 *
 * @param[in] par - parameters struct containing problem information.
 * @param[in] pressure_ib_list - Integer array indicating if a cell is an immersed boundary cell.
 * @param[in] velocity_u - Degrees of freedom for the x-component of the porescale velocity solution.
 * @param[in] velocity_v - Degrees of freedom for the y-component of the porescale velocity solution.
 * @param[in] velocity_w - Degrees of freedom for the z-component of the porescale velocity solution.
 * @param[in] solution - Porescale flow solution.
 */
double
hgf::multiscale::flow::compute_permeability_y(const parameters& par, const std::vector< int >& pressure_ib_list, \
                                                                     const std::vector< degree_of_freedom >& velocity_u, \
                                                                     const std::vector< degree_of_freedom >& velocity_v, \
                                                                     const std::vector< degree_of_freedom >& velocity_w, \
  const std::vector< double > solution)
{
  // quick exit
  if (solution.size() <= 0) {
    std::cout << "\nEmpty flow solution, returning -1 permeability.\n";
    return -1.0;
  }

  double v, g, por;

  compute_averages_y(par, pressure_ib_list, velocity_u, velocity_v, velocity_w, solution, v, g);

  // determine porosity
  int n_voxels = 0;
  int n_void = 0;
  for (int ii = 0; ii < par.voxel_geometry.size(); ii++) {
    n_voxels++;
    if (par.voxel_geometry[ii] == 0) n_void++;
  }
  por = (double)n_void / n_voxels;

  return (v / g) * por;

}

/** \brief hgf::multiscale::flow::compute_permeability_z computes upscaled permeability in the z direction from a porescale flow solution.
 *
 * @param[in] par - parameters struct containing problem information.
 * @param[in] pressure_ib_list - Integer array indicating if a cell is an immersed boundary cell.
 * @param[in] velocity_u - Degrees of freedom for the x-component of the porescale velocity solution.
 * @param[in] velocity_v - Degrees of freedom for the y-component of the porescale velocity solution.
 * @param[in] velocity_w - Degrees of freedom for the z-component of the porescale velocity solution.
 * @param[in] solution - Porescale flow solution.
 */
double
hgf::multiscale::flow::compute_permeability_z(const parameters& par, const std::vector< int >& pressure_ib_list, \
                                                                     const std::vector< degree_of_freedom >& velocity_u, \
                                                                     const std::vector< degree_of_freedom >& velocity_v, \
                                                                     const std::vector< degree_of_freedom >& velocity_w, \
                                                                     const std::vector< double > solution)
{
  // quick exit
  if (solution.size() <= 0) {
    std::cout << "\nEmpty flow solution, returning -1 permeability.\n";
    return -1.0;
  }

  double v, g, por;

  compute_averages_z(par, pressure_ib_list, velocity_u, velocity_v, velocity_w, solution, v, g);

  // determine porosity (holder)
  int n_voxels = 0;
  int n_void = 0;
  for (int ii = 0; ii < par.voxel_geometry.size(); ii++) {
    n_voxels++;
    if (par.voxel_geometry[ii] == 0) n_void++;
  }
  por = (double)n_void / n_voxels;

  return (v / g) * por;
}

/** \brief hgf::multiscale::flow::compute_permeability_tensor computes upscaled permeability tensor given flow solutions for flows in each principal axis direction.
 *
 * @param[in] par - parameters struct containing problem information.
 * @param[in] pressure_ib_list - Integer array indicating if a cell is an immersed boundary cell.
 * @param[in] velocity_u - Degrees of freedom for the x-component of the porescale velocity solution.
 * @param[in] velocity_v - Degrees of freedom for the y-component of the porescale velocity solution.
 * @param[in] velocity_w - Degrees of freedom for the z-component of the porescale velocity solution. Can be empty if par indicates a 2d flow.
 * @param[in] solution_xflow - Solution of flow problem in x-axis direction.
 * @param[in] solution_yflow - Solution of flow problem in y-axis direction.
 * @param[in] solution_zflow - Solution of flow problem in z-axis direction. Can be empty if par indicates a 2d flow.
 * @param[out] permeability - resulting 4x4 (2d flow) or 9x9 (3d flow) permeability tensor. 
 */
void
hgf::multiscale::flow::compute_permeability_tensor(const parameters& par, const std::vector< int >& pressure_ib_list, \
                                                                          const std::vector< degree_of_freedom >& velocity_u, \
                                                                          const std::vector< degree_of_freedom >& velocity_v, \
                                                                          const std::vector< degree_of_freedom >& velocity_w, \
                                                                          const std::vector< double > solution_xflow, \
                                                                          const std::vector< double > solution_yflow, \
                                                                          const std::vector< double > solution_zflow, \
                                                                          std::vector< double >& permeability)
{

  if (par.dimension == 3) {
    double g_val[27], vel[9];
    int g_idx[27], g_jdx[27];

    // Assign Indices
    for (int i = 0; i < 3; i++) {
      // First block
      g_idx[i] = 0;
      g_jdx[i] = i;
      g_idx[3 + i] = 1;
      g_jdx[3 + i] = i;
      g_idx[6 + i] = 2;
      g_jdx[6 + i] = i;
      // Second block
      g_idx[9 + i] = 3;
      g_jdx[9 + i] = 3 + i;
      g_idx[12 + i] = 4;
      g_jdx[12 + i] = 3 + i;
      g_idx[15 + i] = 5;
      g_jdx[15 + i] = 3 + i;
      // Third block
      g_idx[18 + i] = 6;
      g_jdx[18 + i] = 6 + i;
      g_idx[21 + i] = 7;
      g_jdx[21 + i] = 6 + i;
      g_idx[24 + i] = 8;
      g_jdx[24 + i] = 6 + i;
    }

    // Compute averages
    compute_averages_x(par, pressure_ib_list, velocity_u, velocity_v, velocity_w, solution_xflow, vel[0], g_val[0]);
    compute_averages_y(par, pressure_ib_list, velocity_u, velocity_v, velocity_w, solution_xflow, vel[1], g_val[1]);
    compute_averages_z(par, pressure_ib_list, velocity_u, velocity_v, velocity_w, solution_xflow, vel[2], g_val[2]);
    compute_averages_x(par, pressure_ib_list, velocity_u, velocity_v, velocity_w, solution_yflow, vel[3], g_val[3]);
    compute_averages_y(par, pressure_ib_list, velocity_u, velocity_v, velocity_w, solution_yflow, vel[4], g_val[4]);
    compute_averages_z(par, pressure_ib_list, velocity_u, velocity_v, velocity_w, solution_yflow, vel[5], g_val[5]);
    compute_averages_x(par, pressure_ib_list, velocity_u, velocity_v, velocity_w, solution_zflow, vel[6], g_val[6]);
    compute_averages_y(par, pressure_ib_list, velocity_u, velocity_v, velocity_w, solution_zflow, vel[7], g_val[7]);
    compute_averages_z(par, pressure_ib_list, velocity_u, velocity_v, velocity_w, solution_zflow, vel[8], g_val[8]);

    for (int i = 0; i < 9; i++) {
      g_val[i + 9] = g_val[i];
      g_val[i + 18] = g_val[i];
    }

    // column major dense matrix for lapack solve
    double *mat = (double *)calloc(9 * 9, sizeof(double));
    for (int ii = 0; ii < 27; ii++) {
      mat[(g_jdx[ii] * 9 + g_idx[ii])] = g_val[ii];
    }

    // solve linear system for K tensor
    solve_9x9( mat, &vel[0] );

    // determine porosity
    double por;
    int n_voxels = 0;
    int n_void = 0;
    for (int ii = 0; ii < par.voxel_geometry.size(); ii++) {
      n_voxels++;
      if (par.voxel_geometry[ii] == 0) n_void++;
    }
    por = (double)n_void / n_voxels;

    permeability.resize(9);
    for (int ii = 0; ii < 9; ii++) permeability[ii] = por * vel[ii];

  }
  else {
    double g_val[8], vel[4];
    int g_idx[8], g_jdx[8];

    // Assign Indices
    for (int i = 0; i < 2; i++) {
      // First block
      g_idx[i] = 0;
      g_jdx[i] = i;
      g_idx[2 + i] = 1;
      g_jdx[2 + i] = i;
      // Second block
      g_idx[4 + i] = 2;
      g_jdx[4 + i] = 2 + i;
      g_idx[6 + i] = 3;
      g_jdx[6 + i] = 2 + i;
    }

    // Compute averages
    // Compute averages
    compute_averages_x(par, pressure_ib_list, velocity_u, velocity_v, velocity_w, solution_xflow, vel[0], g_val[0]);
    compute_averages_y(par, pressure_ib_list, velocity_u, velocity_v, velocity_w, solution_xflow, vel[1], g_val[1]);
    compute_averages_x(par, pressure_ib_list, velocity_u, velocity_v, velocity_w, solution_yflow, vel[2], g_val[2]);
    compute_averages_y(par, pressure_ib_list, velocity_u, velocity_v, velocity_w, solution_yflow, vel[3], g_val[3]);

    for (int i = 0; i < 4; i++) {
      g_val[i + 4] = g_val[i];
    }


    // column major dense matrix for solve
    double *mat = (double *)calloc(4 * 4, sizeof(double));
    for (int ii = 0; ii < 8; ii++) {
      mat[(g_jdx[ii] * 4 + g_idx[ii])] = g_val[ii];
    }

    // solve linear system for K tensor
    solve_4x4( mat, &vel[0] );

    // determine porosity
    double por;
    int n_voxels = 0;
    int n_void = 0;
    for (int ii = 0; ii < par.voxel_geometry.size(); ii++) {
      n_voxels++;
      if (par.voxel_geometry[ii] == 0) n_void++;
    }
    por = (double)n_void / n_voxels;

    permeability.resize(4);
    for (int ii = 0; ii < 4; ii++) permeability[ii] = por * vel[ii];

  }
}
