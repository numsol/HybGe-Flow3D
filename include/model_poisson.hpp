#ifndef _MODELS_POISSON_H
#define _MODELS_POISSON_H

#include "hgflow.hpp"

namespace hgf
{
  namespace models 
  {
    /** \brief Contains functionality for setup and post-processessing the solution of the Poisson equation in 2d or 3d.
     * 
     */
    class poisson 
    {

      public:

        std::vector< degree_of_freedom > cc_dof;              /**< Degrees of freedom associated with cell centers in the Poisson model. */
        std::vector< int > immersed_boundary;                 /**< If immersed_boundary[i] == 1, then cell i is an immersed boundary cell. */
        std::vector< coo_array > coo_array;                   /**< Linear system associated with discretization of the Poisson equation. */
        std::vector< double > rhs;                            /**< Right-hand side vector (force). */
        std::vector< double > solution;                       /**< Solution vector. */
        void build(const parameters& par, const hgf::mesh::voxel& msh);
        void output_vtk(const parameters& par, const hgf::mesh::voxel& msh, std::string& file_name);
        void write_state(const parameters& par, const hgf::mesh::voxel& msh, std::string& file_name);
        void immersed_boundary(const parameters& par, double eta);
        void import_immersed_boundary(parameters& par, std::vector< int >& input_ib, double eta);       

      private:

        std::vector< boundary_nodes > boundary;                    
        void build_degrees_of_freedom_2d(const parameters& par, const hgf::mesh::voxel& msh);
        void build_array_2d(const parameters& par, const hgf::mesh::voxel& msh);
 
        void build_degrees_of_freedom_3d(const parameters& par, const hgf::mesh::voxel& msh);
        void build_array_3d(const parameters& par, const hgf::mesh::voxel& msh);

    };
  }
}