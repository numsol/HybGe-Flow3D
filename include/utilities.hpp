/* utilities header */
#ifndef _UTILITIES_H
#define _UTILITIES_H

// system includes
#include <vector>
#include <boost/filesystem.hpp>

namespace bfs = boost::filesystem;

/** \brief HGF is a library designed to solve multiscale flow problems in complex domains.
 * 
 */
namespace hgf
{
  void
  init_parameters(parameters& par, const std::string& problem_path);

  /** \brief Contains utility functions.
   *
   */
  namespace utility
  { 
    bool
    find_file( const bfs::path& problem_path, \
               const std::string& file_name, \
               bfs::path& file_path);

    void
    load_parameters(parameters& par, const bfs::path& problem_path);

    void
    print_parameters(parameters& par);

    void
    import_voxel_geometry(parameters& par, const bfs::path& problem_path);

    bool
    check_symmetry(std::vector< array_coo >& array);

    void
    sort_array(std::vector< array_coo >& array);

    void
    unique_array(std::vector< array_coo >& array);
  }

  namespace mesh
  {
    void
    refine_voxel_uniform(parameters& par, int refine_len);

    int
    geo_sanity(parameters& par);

    int
    remove_dead_pores(parameters& par);
  }
}

#endif
