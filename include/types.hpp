#ifndef _TYPES_H
#define _TYPES_H

#include <boost/filesystem.hpp>

/** \brief Struct holding a variety of problem information.
 *
 */
struct parameters
{
  double viscosity;                              /**< For fluid flow problems, the viscosity of the fluid. */
  int dimension;                                 /**< Specifies if the problem is 2d or 3d. */
  double length;                                 /**< Specifies the length of the domain (x-direction). */
  double width;                                  /**< Specifies the width of the domain (y-direction). */
  double height;                                 /**< Specifies the height of the domain (z-direction). */
  int nx;					 /**< Specifies the x mesh dimension. */
  int ny;					 /**< Specifies the y mesh dimension. */
  int nz;					 /**< Specifies the z mesh dimension. */
  int solver_max_iterations;                     /**< Specifies the maximum iterations allowed in iterative solvers. */
  double solver_absolute_tolerance;              /**< Specifies the absolute error tolerance for iterative solvers. */
  double solver_relative_tolerance;              /**< Specifies the relative error tolerance for iterative solvers. */
  int solver_verbose;                            /**< Specifies the level of console output produced by iterative solvers. */
  std::vector< unsigned long > voxel_geometry;   /**< Vector storing a voxel geometry read from the Geometry.dat input file. */
  boost::filesystem::path problem_path;          /**< Path to folder containing Geometry.dat and Parameters.dat input files */
};

/** \brief Mesh vertex struct
 *
 */
struct vertex
{
  double coords[3];                       /**< Array of vertex coordinates. coords[0] gives the x coordinate, coords[1] gives the y coordinate, and coords[2] gives the z coordinate. */
  int gnum;                               /**< Integer specifying the global node number for this vertex. */
};

/** \brief Mesh edge struct
 *
 */
struct edge
{
  int vns[2];                             /**< Array giving the local vertex numbers connected by the edge. */
  int gnum;								                /**< Integer specifying the global edge number for this edge. */
  int neighbor;                           /**< For 2d problems, integer specifying the global cell number for the cell sharing this edge. -1 for boundary edge. */
  int bctype;                             /**< For 2d problems, integer specifying the boundary type for the edge. 0: interior, 1: dirichlet, 2: neumann. */
};

/** \brief Mesh face struct
 *
 */
struct face
{
  int vns[4];                             /**< Array giving the local vertex numbers connected by the face. */
  int gnum;                               /**< Integer specifying the global face number for this face. */
  int neighbor;                           /**< Integer specifying the global cell number for the cell sharing this face. -1 for a boundary face. */
  int bctype;                             /**< Integer specifying the boundary type for the face. 0: interior, 1: dirichlet, 2: neumann. */
};

/** \brief Hexahedral or quadrilateral cell struct 
 *
 */
struct qcell
{
  struct edge edg[12];                         /**< Array of edge structs detailing the edges in the cell. For a 2d problem, edg[0:3] desribes all edges. */
  struct face fac[6];                          /**< Array of face structs detailing the faces in the cell. Not used in 2d problems. */
  struct vertex vtx[8];                        /**< Array of vertex structs detailing the vertices in the cell. For a 2d problem, vtx[0:3] describes all vertices. */
  double dx;                                   /**< Length of the cell (x direction). */
  double dy;                                   /**< Width of the cell (y direction). */
  double dz;                                   /**< Height of the cell (z direction). */
};

/** \brief Struct describing a degree of freedom in a structured model.
 *
 */
struct degree_of_freedom
{
  int doftype;                   /**< Integer specifiying the type of degree of freedom: 0 for interior, 1 for edge, 2 for face. */
  double coords[3];              /**< Array of coordinates of the degree of freedom. coords[0] gives the x coordinate, coords[1] gives the y coordinate, and coords[2] gives the z coordinate. */
  int cell_numbers[2];           /**< Array of mesh cell numbers containing the degree of freedom. */
  int neighbors[6];              /**< Array listing the global number of neighboring degrees of freedom, i.e. DOFs which interact with this DOF in the model. */
};

/** \brief Struct describing a degree of freedom in an unstructured model. 
 *
 * To account for the lack of structure std::vectors are used in place of static arrays. Structure models should use degree_of_freedom instead.
 */
struct dof_unstructured
{
  int doftype;                     /**< Integer specifiying the type of degree of freedom: 0 for interior, 1 for edge, 2 for face. */
  double coords[3];                /**< Array of coordinates of the degree of freedom. coords[0] gives the x coordinate, coords[1] gives the y coordinate, and coords[2] gives the z coordinate. */
  std::vector<double> normals;     /**< Array containing normal direction vectors for faces (3d) or edges (2d). */
  std::vector<int> cell_numbers;   /**< Array of mesh cell numbers containing the degree of freedom. */
  std::vector<int> neighbors;      /**< Array listing the global number of neighboring degrees of freedom, i.e. DOFs which interact with this DOF in the model. */
};

/** \brief Struct for coordinate sparse data format.
 *
 */
struct array_coo
{
  int i_index;     /**< Integer determining the row of the array entry. */
  int j_index;     /**< Integer determining the column of the array entry. */
  double value;    /**< Double precision value of the array entry. */
};

/** \brief Struct for sorting an array of degrees of freedom by x-coordinate then by y-coordinate (y increases fastest). 
 *
 * z-coordinate is not considered. Intended use is for 2d models.
 */
struct byXbyY
{
  /** \brief Operator returns true if x coordinate of first argument is less than that of second, or if the x coordinate are equal
   *         and the y coordinate of the first argument is less than that of the second.
   */
  bool operator()(degree_of_freedom const &one, degree_of_freedom const &two)
  {
    return (one.coords[0] < two.coords[0] || \
      (one.coords[0] == two.coords[0] && one.coords[1] < two.coords[1]));
  }
};

/** \brief Struct for sorting an array of degrees of freedom by z-coordinate.
 *
 * y-coordinate increases fastest, then x-coordinate, then z-coordinate.
 */
struct byZbyXbyY
{
  /** \brief Operator returns true if z coordinate of first argument is less than that of second, or if the z coordinates are equal
   *         and the x coordinate of the first argument is less than that of the second, or if the z and x coordinates
   *         are equal and the y coordinate of the first argument ie less than that of the second.
   */
  bool operator()(degree_of_freedom const &one, degree_of_freedom const &two)
  {
    return (one.coords[2] < two.coords[2] || (one.coords[2] == two.coords[2] && one.coords[0] < two.coords[0]) \
      || (one.coords[2] == two.coords[2] && one.coords[0] == two.coords[0] && one.coords[1] < two.coords[1]));
  }
};

/** \brief Struct for sorting an array of degrees of freedom by y-coordinate, then by x-coordinate, then by z-coordinate.
 *
 * z-coordinate increases fastest, then x-coordinate, then y-coordinate.
 */
struct byYbyXbyZ
{
  /** \brief Operator returns true if y coordinate of first argument is less than that of second, or if the y coordinates are equal
   *         and the x coordinate of the first argument is less than that of the second, or if the y and x coordinates
   *         are equal and the z coordinate of the first argument ie less than that of the second.
   */
  bool operator()(degree_of_freedom const &one, degree_of_freedom const &two)
  {
    return (one.coords[1] < two.coords[1] || (one.coords[1] == two.coords[1] && one.coords[0] < two.coords[0]) \
      || (one.coords[1] == two.coords[1] && one.coords[0] == two.coords[0] && one.coords[2] < two.coords[2]));
  }
};

/** \brief Struct for sorting a coordinate sparse array by row index and then by column index. 
 *
 * Column index increases fastest, then row index.
 */
struct byIbyJ
{
  /** \brief Operator returns true if the row index of the first argument is less than that of the second, or if the row indices
   *         are equal and the column index of the first argument is less than that of the second.
   */
  bool operator()(array_coo const &one, array_coo const &two)
  {
    return (one.i_index < two.i_index || \
      (one.i_index == two.i_index && one.j_index < two.j_index));
  }
};

/** \brief Struct for sorting a coordinate sparse array by column index then by row index. 
 *
 * Row index increases fastest, then column index.
 */
struct byJbyI
{
  /** \brief Operator returns true if the column index of the first argument is less than that of the second, or if the column indices
   *         are equal and the row index of the first argument is less than that of the second.
   */
  bool operator()(array_coo const &one, array_coo const &two)
  {
    return (one.j_index < two.j_index || \
      (one.j_index == two.j_index && one.i_index < two.i_index));
  }
};

/** \brief Struct for tracking boundary node information.
 *
 */
struct boundary_nodes
{
  int type;              /**< Boundary type. */
  double value;          /**< Boundary value. */
};

#endif
