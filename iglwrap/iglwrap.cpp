#include <Eigen/Dense>
#include <stdlib.h>
#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <igl/copyleft/cgal/piecewise_constant_winding_number.h>
#include "iglwrap.h"
#include <boost/bind/bind.hpp>
using namespace boost::placeholders;
using namespace Eigen;

typedef Matrix<double,Dynamic,3> VertexMatrix;
typedef Matrix<int,Dynamic,3> FaceMatrix;

typedef Matrix<double, Dynamic, 3, RowMajor> vertices_t;
typedef Matrix<int, Dynamic, 3, RowMajor> faces_t;

int igl_mesh_boolean(
  int op,
  int nv1, int nf1, double *mv1, int *mf1,
  int nv2, int nf2, double *mv2, int *mf2,
  int *nv3, int *nf3, double **mv3, int **mf3, int **index) {
  Map<vertices_t> v1(mv1, nv1, 3), v2(mv2, nv2, 3);
  Map<faces_t> f1(mf1, nf1, 3), f2(mf2, nf2, 3);
  vertices_t v3;
  faces_t f3;
  VectorXi j;

//   std::cout << "matrix v1:\n" << v1 << "\nmatrix f1:" << f1 << "\n";
//   std::cout << "matrix v2:\n" << v2 << "\nmatrix f2:" << f2 << "\n";
//   std::cout << "calling mesh_boolean\n";
  igl::copyleft::cgal::mesh_boolean(v1,f1,v2,f2,igl::MeshBooleanType(op),
            v3, f3, j);
//   std::cout << "done!\n";
//   std::cout << "f3 = " << f3 << "\n";
//   std::cout << j << "\n";
  *mv3 = (double *) malloc(3*v3.rows()*sizeof(double));
  if(*mv3 == 0) {
    return -1;
  }
  *mf3 = (int *) malloc(3*f3.rows()*sizeof(int));
  if(*mf3 == 0) {
    free(mv3);
    return -1;
  }
  *index = (int *) malloc(f3.rows()*sizeof(int));
  if (*index == 0) {
    free(mv3);
    free(mf3);
    return -1;
  }
  memcpy(*mv3, v3.data(), 3*v3.rows()*sizeof(double));
  memcpy(*mf3, f3.data(), 3*f3.rows()*sizeof(int));
  memcpy(*index, j.data(), f3.rows()*sizeof(int));
  *nv3 = v3.rows();
  *nf3 = f3.rows();
//   std::cout << "copy operations done\n";
  return 0;
}

int igl_mesh_is_pwn(int nv, int nf, double *mv, int *mf) {
  Matrix<double, Dynamic, 3> v = Map<vertices_t> (mv, nv, 3);
  Matrix<int, Dynamic, 3> f = Map<faces_t> (mf, nf, 3);
  return (int)igl::copyleft::cgal::piecewise_constant_winding_number(v, f);
}
