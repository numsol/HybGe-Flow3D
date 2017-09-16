/* Example solves the homogeneous Poisson equation. The solution is saved for visualization. Build with included
   CMakeLists.txt, and use:
     poisson_homogeneous <path/to/problemfolder>
   Some example problem folders are included at examples/geometries.
*/

#include <vector>
#include <iostream>
#include <stdlib.h>
#include <omp.h>

#include "hgflow.hpp"

int
main( int argc, const char* argv[] )
{
  //----- POISSON example -----//
  std::cout << "\n//----Solving POISSON---//\n";
  //-- timers --//
  double begin, rebegin, para_time, mesh_time, build_time, solve_time, postp_time, total_time;

  begin = omp_get_wtime();
  //--- problem parameters ---//
  parameters par;
  hgf::init_parameters(par, argv[1]);

  para_time = omp_get_wtime() - begin;
  rebegin = omp_get_wtime();

  //--- mesh ---//
  // check mesh sanity
  hgf::mesh::geo_sanity(par);
  // build the mesh
  hgf::mesh::voxel msh;
  msh.build(par);

  mesh_time = omp_get_wtime() - rebegin;
  rebegin = omp_get_wtime();

  //--- poisson model ---//
  hgf::models::poisson poiss;
  // set alpha before build
  poiss.set_constant_scalar_alpha(par, 1.0);
  // build the degrees of freedom and the array -- if alpha is not already set, initializes to identity on each cell.
  poiss.build(par, msh);

  // set up boundary conditions
  poiss.setup_homogeneous_dirichlet_bc(par, msh);

  // force term set to a constant
  poiss.set_constant_force(par, 1);

  build_time = omp_get_wtime() - rebegin;
  rebegin = omp_get_wtime();

  // solve with Paralution
  hgf::solve::paralution::init_solver();
  hgf::solve::paralution::solve(par, poiss.coo_array, poiss.rhs, poiss.solution);
  hgf::solve::paralution::finalize_solver();

  solve_time = omp_get_wtime() - rebegin;
  rebegin = omp_get_wtime();

  // save the solution for visualization with paraview
  std::string file_name = "Solution";
  poiss.output_vtk(par, msh, file_name);

  postp_time = omp_get_wtime() - rebegin;
  total_time = omp_get_wtime() - begin;

  // print timings
  std::cout << "\n\n//------------------ Poisson Finished! Time reports ------------------//\n";
  std::cout << "Parameter setup time: " << para_time << "\n";
  std::cout << "Mesh time: " << mesh_time << "\n";
  std::cout << "Array construction time: " << build_time << "\n";
  std::cout << "Solver time: " << solve_time << "\n";
  std::cout << "Post processessing time: " << postp_time << "\n";
  std::cout << "Total time: " << total_time << "\n";

}
