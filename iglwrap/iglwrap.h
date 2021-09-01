#ifdef __cplusplus
extern "C" {
#endif

/* A mesh is represented as (nv, nf, mv, mf), where
 * int nv = number of vertices
 * int nf = number of faces
 * double *mv = matrix of vertex coordinates (row-wise)
 * int *mf = matrix of triangles (row-wise)
 *
 * Operation is:
 * 0: union
 * 1: intersection
 * 2: minus
 * 3: xor
 *
 * Places malloc()ed arrays (in row-major order) as mv3 and mf3.
 * Returns -1 on failure and 0 on success.
 */

int igl_mesh_boolean(
  int op,
  int nv1, int nf1, double *mv1, int *mf1,
  int nv2, int nf2, double *mv2, int *mf2,
  int *nv3, int *nf3, double **mv3, int **mf3, int **index);

int igl_mesh_is_pwn(int nv, int nf, double *mv, int *mf);

#ifdef __cplusplus
} // extern
#endif
