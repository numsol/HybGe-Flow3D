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
  int NTHREADS = omp_get_max_threads();
  int block_size = ((int)cc_dof.size() % NTHREADS) ? (int)((cc_dof.size() / NTHREADS) + 1) : (int)(cc_dof.size() / NTHREADS);

  std::vector< std::vector< array_coo > > temp_arrays;
  temp_arrays.resize(NTHREADS);
  
}