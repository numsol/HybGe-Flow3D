/* poisson dof 2d source */

// hgf includes
#include "model_poisson.hpp"

void
hgf::models::poisson::build_degrees_of_freedom_2d(const parameters& par, const hgf::mesh::voxel& msh)
{
  phi.resize(msh.els.size());

#pragma omp parallel for
  for (int cell = 0; cell < msh.els.size(); cell++) {
    degree_of_freedom dof_temp;
    // coordinates
    dof_temp.coords[0] = 0.25 * \
      (msh.els[cell].vtx[0].coords[0] + msh.els[cell].vtx[1].coords[0] + \
        msh.els[cell].vtx[2].coords[0] + msh.els[cell].vtx[3].coords[0]);
    dof_temp.coords[1] = 0.25 * \
       (msh.els[cell].vtx[0].coords[1] + msh.els[cell].vtx[1].coords[1] + \
         msh.els[cell].vtx[2].coords[1] + msh.els[cell].vtx[3].coords[1]);
    // doftype
    dof_temp.doftype = 0;
    // cell numbers
    dof_temp.cell_numbers[0] = cell;
    dof_temp.cell_numbers[1] = -1;
    // neighbors
    for (int nbr = 0; nbr < 4; nbr++) {
      dof_temp.neighbors[nbr] = msh.els[cell].edg[nbr].neighbor;
    }
    // push_back
    phi[cell] = dof_temp;
  }
}

