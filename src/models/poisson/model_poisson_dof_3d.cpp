/* poisson dof 3d source */

// hgf includes
#include "model_poisson.hpp"

void
hgf::models::poisson::build_degrees_of_freedom_3d(const parameters& par, const hgf::mesh::voxel& msh)
{

  phi.resize(msh.els.size());

#pragma omp parallel for
  for (int cell = 0; cell < msh.els.size(); cell++) {
    degree_of_freedom dof_temp;
    // coordinates
    for (int dir = 0; dir < 3; dir++) {
      dof_temp.coords[dir] = 0;
      for (int ii = 0; ii < 8; ii++) dof_temp.coords[dir] += msh.els[cell].vtx[ii].coords[dir];
      dof_temp.coords[dir] = dof_temp.coords[dir] / 8;
    }
    // doftype
    dof_temp.doftype = 0;
    // cell numbers
    dof_temp.cell_numbers[0] = cell;
    dof_temp.cell_numbers[1] = -1;
    // neighbors
    for (int nbr = 0; nbr < 6; nbr++) dof_temp.neighbors[nbr] = msh.els[cell].fac[nbr].neighbor;
    // place phi dof
    phi[cell] = dof_temp;
  }
}


