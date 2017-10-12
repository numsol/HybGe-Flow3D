/* Example computes the upscaled permeability in the x-direction from the solution to Stokes and prints
   the result to console. The x-flow solution is saved for visualization. Build with included
   CMakeLists.txt, and use:
     permeability <path/to/problemfolder>
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
  //----- STOKES example -----//
  std::cout << "\n//----Solving STOKES----//\n";
  //-- timers --//
  double begin, rebegin, para_time, mesh_time, build_time, solve_time, postp_time, total_time;

  begin = omp_get_wtime();
  //--- problem parameters ---//
  parameters par;
  hgf::init_parameters(par, argv[1]);

  para_time = omp_get_wtime() - begin;
  rebegin = omp_get_wtime();

  //--- mesh ---//
  // check mesh sanity and remove dead pores
  hgf::mesh::geo_sanity(par);
  int pores_removed = 0;
  pores_removed = hgf::mesh::remove_dead_pores(par);
  if (pores_removed > 0) std::cout << "\n" << pores_removed << " dead pores removed.\n";
  // build the mesh
  hgf::mesh::voxel msh;
  msh.build(par);

  mesh_time = omp_get_wtime() - rebegin;
  rebegin = omp_get_wtime();

  //--- stokes model ---//
  hgf::models::stokes x_stks;
  // build the degrees of freedom and the array for interior cells, initializes viscosity to 1.0
  x_stks.build(par, msh);
  // if viscosity != 1, set after build
  // x_stks.viscosity = 0.5;

  // add penalty for immersed boundary cells
  double eta = 1e-5;
  x_stks.immersed_boundary(par, eta);

  // set up boundary conditions
  hgf_inflow inflow = HGF_INFLOW_CONSTANT;
  x_stks.setup_xflow_bc(par, msh, inflow);

  build_time = omp_get_wtime() - rebegin;
  rebegin = omp_get_wtime();

  // solve with Paralution
  hgf::solve::paralution::init_solver();
  // block diagonal preconditioner for 3d problem
  if (par.dimension == 3) { 
    hgf::solve::paralution::solve_ps_flow(par, x_stks.coo_array, x_stks.rhs, x_stks.solution_int, \
      (int)std::accumulate(x_stks.interior_u.begin(), x_stks.interior_u.end(), 0), \
      (int)std::accumulate(x_stks.interior_v.begin(), x_stks.interior_v.end(), 0), \
      (int)std::accumulate(x_stks.interior_w.begin(), x_stks.interior_w.end(), 0), \
      (int)x_stks.pressure.size());
  }
  // simple GMRES + ILU for 2d
  else { 
    hgf::solve::paralution::solve(par, x_stks.coo_array, x_stks.rhs, x_stks.solution_int);
  }
  hgf::solve::paralution::finalize_solver();

  solve_time = omp_get_wtime() - rebegin;
  rebegin = omp_get_wtime();

  // fill in solution vector with solution of linear system + boundary values
  x_stks.solution_build();  

  // compute permeability and print to console
  double permeability = hgf::multiscale::flow::compute_permeability_x( par, \
    x_stks.pressure_ib_list, x_stks.velocity_u, x_stks.velocity_v, x_stks.velocity_w, x_stks.solution );
  std::cout << "\nX permeability = " << permeability << "\n";

  // save the x-flow solution for visualization with paraview
  std::string file_name = "Solution_x_inflow";
  x_stks.output_vtk(par, msh, file_name);

  // save the full state of the flow simulation to .dat file
  std::string state_name = "Stokes_State";
  x_stks.write_state(par, msh, state_name);

  postp_time = omp_get_wtime() - rebegin;
  total_time = omp_get_wtime() - begin;

  // print timings
  std::cout << "\n\n//------------------ Stokes Finished! Time reports ------------------//\n";
  std::cout << "Parameter setup time: " << para_time << "\n";
  std::cout << "Mesh time: " << mesh_time << "\n";
  std::cout << "Array construction time: " << build_time << "\n";
  std::cout << "Solver time: " << solve_time << "\n";
  std::cout << "Post processessing time: " << postp_time << "\n";
  std::cout << "Total time: " << total_time << "\n";

}
