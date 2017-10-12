#ifndef _MODELS_POISSON_H
#define _MODELS_POISSON_H

#include "hgflow.hpp"

// system includes
#include <vector>
#include <array>
#include <iostream>
#include <omp.h>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <boost/filesystem.hpp>

#define mean_perm(aval, bval) (0.5 * (aval + bval))

namespace hgf
{
  /**
   * \brief Contains functions and classes that setup and post-process flow models.
   *
   */
  namespace models
  {
    /** \brief Contains functionality for setup and post-processessing the solution of the Poisson equation in 2d or 3d.
     * 
     */
    class poisson
    {

      public:
      
        std::vector< degree_of_freedom > phi;                         /**< Degrees of freedom. */
        std::vector< array_coo > coo_array;                           /**< Linear system associated to the Poisson problem stored in COO (coordinate) sparse format */
        std::vector< double > rhs;                                    /**< Right-hand side vector (force). */
        std::vector< double > solution;                               /**< Vector for storing full solution. */
        std::vector< std::vector< double > > alpha;                   /**< Coefficient tensor. */
        void build(const parameters& par, const hgf::mesh::voxel& msh);
        void output_vtk(const parameters& par, const hgf::mesh::voxel& msh, std::string& file_name);
        void set_constant_force(const parameters& par, const double& force_in);
        void set_constant_scalar_alpha(const parameters& par, const hgf::mesh::voxel& msh, const double& alpha_in);
        void set_constant_tensor_alpha(const parameters& par, const hgf::mesh::voxel& msh, const std::vector< double >& alpha_in);
        void setup_homogeneous_dirichlet_bc(const parameters& par, const hgf::mesh::voxel& msh);
        void set_nonhomogeneous_dirichlet_bc(const parameters& par, const hgf::mesh::voxel& msh, double (*f)( int dof_num, double coords[3] ) );
        void add_to_force( const std::vector< double >& rhs_to_add );
    
      private:

        int NTHREADS, block_size;
        
        void build_degrees_of_freedom_2d(const parameters& par, const hgf::mesh::voxel& msh);
        void dof_neighbors_2d(const parameters& par, const hgf::mesh::voxel& msh);
        void build_array_2d(const parameters& par, const hgf::mesh::voxel& msh);
        void homogeneous_dirichlet_2d(const parameters& par, const hgf::mesh::voxel& msh);

        void build_degrees_of_freedom_3d(const parameters& par, const hgf::mesh::voxel& msh);
        void dof_neighbors_3d(const parameters& par, const hgf::mesh::voxel& msh);
        void build_array_3d(const parameters& par, const hgf::mesh::voxel& msh);
        void homogeneous_dirichlet_3d(const parameters& par, const hgf::mesh::voxel& msh);
    };
  }
}

#endif
