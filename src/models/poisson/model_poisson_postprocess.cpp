/* poisson postprocess source */

// hgf includes
#include "model_poisson.hpp"
#include <iomanip>

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

/** \brief hgf::models::poisson::output_vtk saves the solution to the Poisson equation to a file for VTK visualiztion.
 *
 * @param[in] par - parameters struct containing problem information, including problem directory.
 * @param[in] msh - mesh object containing a quadrilateral or hexagonal representation of geometry from problem folder addressed in parameters& par.
 * @param[in,out] file_name - string used to name the output file, which is placed in the problem directory contained in parameters& par.
 */
void
hgf::models::poisson::output_vtk(const parameters& par, const hgf::mesh::voxel& msh, std::string& file_name)
{

  if (par.dimension == 3) { // 3d output
    int nNodes = (int)msh.gtlNode.size() / 8;
    int nEls = (int)msh.els.size();
    // build an exclusive nodes vector
    std::vector<double> nodes(nNodes * 3);
#pragma omp parallel for
    for (int ii = 0; ii < nNodes; ii++) {
      if (msh.gtlNode[idx2(ii, 0, 8)]) {
        nodes[idx2(ii, 0, 3)] = msh.els[msh.gtlNode[idx2(ii, 0, 8)] - 1].vtx[5].coords[0];
        nodes[idx2(ii, 1, 3)] = msh.els[msh.gtlNode[idx2(ii, 0, 8)] - 1].vtx[5].coords[1];
        nodes[idx2(ii, 2, 3)] = msh.els[msh.gtlNode[idx2(ii, 0, 8)] - 1].vtx[5].coords[2];
      }
      else if (msh.gtlNode[idx2(ii, 1, 8)]) {
        nodes[idx2(ii, 0, 3)] = msh.els[msh.gtlNode[idx2(ii, 1, 8)] - 1].vtx[4].coords[0];
        nodes[idx2(ii, 1, 3)] = msh.els[msh.gtlNode[idx2(ii, 1, 8)] - 1].vtx[4].coords[1];
        nodes[idx2(ii, 2, 3)] = msh.els[msh.gtlNode[idx2(ii, 1, 8)] - 1].vtx[4].coords[2];
      }
      else if (msh.gtlNode[idx2(ii, 2, 8)]) {
        nodes[idx2(ii, 0, 3)] = msh.els[msh.gtlNode[idx2(ii, 2, 8)] - 1].vtx[7].coords[0];
        nodes[idx2(ii, 1, 3)] = msh.els[msh.gtlNode[idx2(ii, 2, 8)] - 1].vtx[7].coords[1];
        nodes[idx2(ii, 2, 3)] = msh.els[msh.gtlNode[idx2(ii, 2, 8)] - 1].vtx[7].coords[2];

      }
      else if (msh.gtlNode[idx2(ii, 3, 8)]) {
        nodes[idx2(ii, 0, 3)] = msh.els[msh.gtlNode[idx2(ii, 3, 8)] - 1].vtx[6].coords[0];
        nodes[idx2(ii, 1, 3)] = msh.els[msh.gtlNode[idx2(ii, 3, 8)] - 1].vtx[6].coords[1];
        nodes[idx2(ii, 2, 3)] = msh.els[msh.gtlNode[idx2(ii, 3, 8)] - 1].vtx[6].coords[2];

      }
      else if (msh.gtlNode[idx2(ii, 4, 8)]) {
        nodes[idx2(ii, 0, 3)] = msh.els[msh.gtlNode[idx2(ii, 4, 8)] - 1].vtx[1].coords[0];
        nodes[idx2(ii, 1, 3)] = msh.els[msh.gtlNode[idx2(ii, 4, 8)] - 1].vtx[1].coords[1];
        nodes[idx2(ii, 2, 3)] = msh.els[msh.gtlNode[idx2(ii, 4, 8)] - 1].vtx[1].coords[2];

      }
      else if (msh.gtlNode[idx2(ii, 5, 8)]) {
        nodes[idx2(ii, 0, 3)] = msh.els[msh.gtlNode[idx2(ii, 5, 8)] - 1].vtx[0].coords[0];
        nodes[idx2(ii, 1, 3)] = msh.els[msh.gtlNode[idx2(ii, 5, 8)] - 1].vtx[0].coords[1];
        nodes[idx2(ii, 2, 3)] = msh.els[msh.gtlNode[idx2(ii, 5, 8)] - 1].vtx[0].coords[2];

      }
      else if (msh.gtlNode[idx2(ii, 6, 8)]) {
        nodes[idx2(ii, 0, 3)] = msh.els[msh.gtlNode[idx2(ii, 6, 8)] - 1].vtx[3].coords[0];
        nodes[idx2(ii, 1, 3)] = msh.els[msh.gtlNode[idx2(ii, 6, 8)] - 1].vtx[3].coords[1];
        nodes[idx2(ii, 2, 3)] = msh.els[msh.gtlNode[idx2(ii, 6, 8)] - 1].vtx[3].coords[2];

      }
      else if (msh.gtlNode[idx2(ii, 7, 8)]) {
        nodes[idx2(ii, 0, 3)] = msh.els[msh.gtlNode[idx2(ii, 7, 8)] - 1].vtx[2].coords[0];
        nodes[idx2(ii, 1, 3)] = msh.els[msh.gtlNode[idx2(ii, 7, 8)] - 1].vtx[2].coords[1];
        nodes[idx2(ii, 2, 3)] = msh.els[msh.gtlNode[idx2(ii, 7, 8)] - 1].vtx[2].coords[2];

      }
    }
    // write to vtk file
    bfs::path output_path(par.problem_path / file_name.c_str());
    output_path += ".vtk";
    std::ofstream outstream;
    outstream.open(output_path.string());
    outstream << "# vtk DataFile Version 3.0\n";
    outstream << "vtk output\n";
    outstream << "ASCII\n\n";
    outstream << "DATASET UNSTRUCTURED_GRID\n";
    outstream << "POINTS " << nNodes << " double\n";
    for (int row = 0; row < nNodes; row++) {
      outstream << nodes[idx2(row, 0, 3)] << "\t";
      outstream << nodes[idx2(row, 1, 3)] << "\t";
      outstream << nodes[idx2(row, 2, 3)] << "\n";
    }
    outstream << "\n";
    outstream << "CELLS " << nEls << " " << 9 * nEls << "\n";
    for (int row = 0; row < nEls; row++) {
      outstream << 8 << "\t";
      outstream << msh.els[row].vtx[0].gnum << "\t";
      outstream << msh.els[row].vtx[1].gnum << "\t";
      outstream << msh.els[row].vtx[2].gnum << "\t";
      outstream << msh.els[row].vtx[3].gnum << "\t";
      outstream << msh.els[row].vtx[7].gnum << "\t";
      outstream << msh.els[row].vtx[6].gnum << "\t";
      outstream << msh.els[row].vtx[5].gnum << "\t";
      outstream << msh.els[row].vtx[4].gnum << "\t";
    }
    outstream << "\n";
    outstream << "CELL_TYPES " << nEls << "\n";
    for (int row = 0; row < nEls; row++) {
      outstream << 12 << "\n";
    }
    outstream << "\n";
    outstream << "CELL_DATA " << nEls << "\n";
    outstream << "SCALARS phi double\n";
    outstream << "LOOKUP_TABLE default\n";
    for (int row = 0; row < nEls; row++) {
      outstream << solution[row] << "\n";
    }
    outstream.close();
  }
  else { // 2d output
    int nNodes = (int)msh.gtlNode.size() / 4;
    int nEls = (int)msh.els.size();
    // build an exclusive nodes vector
    std::vector<double> nodes(nNodes * 2);
#pragma omp parallel for
    for (int ii = 0; ii < nNodes; ii++) {
      if (msh.gtlNode[idx2(ii, 0, 4)]) {
        nodes[idx2(ii, 0, 2)] = msh.els[msh.gtlNode[idx2(ii, 0, 4)] - 1].vtx[2].coords[0];
        nodes[idx2(ii, 1, 2)] = msh.els[msh.gtlNode[idx2(ii, 0, 4)] - 1].vtx[2].coords[1];
      }
      else if (msh.gtlNode[idx2(ii, 1, 4)]) {
        nodes[idx2(ii, 0, 2)] = msh.els[msh.gtlNode[idx2(ii, 1, 4)] - 1].vtx[3].coords[0];
        nodes[idx2(ii, 1, 2)] = msh.els[msh.gtlNode[idx2(ii, 1, 4)] - 1].vtx[3].coords[1];
      }
      else if (msh.gtlNode[idx2(ii, 2, 4)]) {
        nodes[idx2(ii, 0, 2)] = msh.els[msh.gtlNode[idx2(ii, 2, 4)] - 1].vtx[1].coords[0];
        nodes[idx2(ii, 1, 2)] = msh.els[msh.gtlNode[idx2(ii, 2, 4)] - 1].vtx[1].coords[1];
      }
      else if (msh.gtlNode[idx2(ii, 3, 4)]) {
        nodes[idx2(ii, 0, 2)] = msh.els[msh.gtlNode[idx2(ii, 3, 4)] - 1].vtx[0].coords[0];
        nodes[idx2(ii, 1, 2)] = msh.els[msh.gtlNode[idx2(ii, 3, 4)] - 1].vtx[0].coords[1];
      }
    }
    // write solution vtk file
    bfs::path output_path( par.problem_path / file_name.c_str() );
    output_path += ".vtk";
    std::ofstream outstream;
    outstream.open(output_path.string());
    outstream << "# vtk DataFile Version 3.0\n";
    outstream << "vtk output\n";
    outstream << "ASCII\n\n";
    outstream << "DATASET UNSTRUCTURED_GRID\n";
    outstream << "POINTS " << nNodes << " double\n";
    for (int row = 0; row < nNodes; row++) {
      outstream << nodes[idx2(row, 0, 2)] << "\t";
      outstream << nodes[idx2(row, 1, 2)] << "\t";
      outstream << 0.0 << "\n";
    }
    outstream << "\n";
    outstream << "CELLS " << nEls << " " << 5 * nEls << "\n";
    for (int row = 0; row < nEls; row++) {
      outstream << 4 << "\t";
      outstream << msh.els[row].vtx[0].gnum << "\t";
      outstream << msh.els[row].vtx[1].gnum << "\t";
      outstream << msh.els[row].vtx[2].gnum << "\t";
      outstream << msh.els[row].vtx[3].gnum << "\t";
    }
    outstream << "\n";
    outstream << "CELL_TYPES " << nEls << "\n";
    for (int row = 0; row < nEls; row++) {
      outstream << 9 << "\n";
    }
    outstream << "\n";
    outstream << "CELL_DATA " << nEls << "\n";
    outstream << "SCALARS phi double\n";
    outstream << "LOOKUP_TABLE default\n";
    for (int row = 0; row < nEls; row++) {
      outstream << solution[row] << "\n";
    }
    outstream.close();
  }
}


