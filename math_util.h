/* -----------------------------------------------------------------------
   See COPYRIGHT.TXT and LICENSE.TXT for copyright and license information
   ----------------------------------------------------------------------- */
#ifndef _math_util_h_
#define _math_util_h_

#include <float.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif
#ifndef M_SQRT2
#define M_SQRT2         1.41421356237309504880
#endif
#ifndef M_SQRTPI
#define M_SQRTPI        1.77245385090551602792981
#endif
#ifndef M_TWOPI
#define M_TWOPI         (M_PI * 2.0)
#endif
#ifndef DBL_MAX
#define DBL_MAX         (1E+37)
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2       0.70710678118654752440  /* 1/sqrt(2) */
#endif

/* Returns integer data type */
#define ROUND_INT(x) (((x) >= 0) ? ((long)((x)+0.5)) : (long)(-(-(x)+0.5)))

/* Returns double data type -- note MSVC does not have C99 round(). */
#define ROUND(x) ((double) (ROUND_INT(x)))

/* Returns +1 or -1, depeding on sign.  Zero yeilds +1. */
#define SIGN(x) (((x) >= 0) ? (+1) : (-1))


/* Primatives */
static inline void vec2_add2 (double* v1, const double* v2) {
    v1[0] += v2[0]; v1[1] += v2[1];
}

template <class T>
static inline void vec3_add2 (T* v1, const T* v2) {
    v1[0] += v2[0]; v1[1] += v2[1]; v1[2] += v2[2];
}

template <class T>
static inline void vec3_add3 (T* v1, const T* v2, const T* v3) {
    v1[0] = v2[0] + v3[0]; v1[1] = v2[1] + v3[1]; v1[2] = v2[2] + v3[2];
}

template <class T>
static inline void vec3_copy (T* v1, const T* v2) {
    v1[0] = v2[0]; v1[1] = v2[1]; v1[2] = v2[2];
}

static inline void vec4_copy (double* v1, const double* v2) {
    v1[0] = v2[0]; v1[1] = v2[1]; v1[2] = v2[2]; v1[3] = v2[3];
}

template <class T>
static inline T vec3_dot (const T* v1, const T* v2) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

static inline double vec4_dot (const double* v1, const double* v2) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3];
}

template <class T>
static inline void vec3_scale2 (T* v1, T a) {
    v1[0] *= a; v1[1] *= a; v1[2] *= a;
}

template <class T>
static inline void vec3_scale3 (T* v1, const T* v2, T a) {
    v1[0] = a * v2[0]; v1[1] = a * v2[1]; v1[2] = a * v2[2];
}

static inline void vec3_sub2 (double* v1, const double* v2) {
    v1[0] -= v2[0]; v1[1] -= v2[1]; v1[2] -= v2[2];
}

template <class T>
static inline void vec3_sub3 (T* v1, const T* v2, const T* v3) {
    v1[0] = v2[0] - v3[0]; v1[1] = v2[1] - v3[1]; v1[2] = v2[2] - v3[2];
}

template <class T>
static inline void vec3_invert (T* v1) {
    vec3_scale2 (v1, T(-1.0));
}

template <class T>
static inline void vec_zero (T* v1, int n) {
    memset (v1, 0, n*sizeof(T));
}

/* Length & distance */
template <class T>
static inline T vec3_len (const T* v1) {
    return sqrt(vec3_dot(v1,v1));
}

template <class T>
static inline void vec3_normalize1 (T* v1) {
    vec3_scale2 (v1, 1 / vec3_len(v1));
}

static inline double vec3_dist (const double* v1, const double* v2) {
    double tmp[3];
    vec3_sub3 (tmp, v1, v2);
    return vec3_len(tmp);
}

/* Cross product */
template <class T>
static inline void vec3_cross (T* v1, const T* v2, const T* v3)
{
    v1[0] = v2[1] * v3[2] - v2[2] * v3[1];
    v1[1] = v2[2] * v3[0] - v2[0] * v3[2];
    v1[2] = v2[0] * v3[1] - v2[1] * v3[0];
}

/* Outer product */
static inline void vec_outer (double* v1, const double* v2, const double* v3, const int n)
{
    int i,j;
    for (j=0; j<n; j++) {
        for (i=0; i<n; i++) {
            v1[n*j + i] = v2[j] * v3[i];
        }
    }
}

/* Matrix ops */

/* Matrix element m[i,j] for matrix with c columns */
#define m_idx(m1,c,i,j) m1[i*c+j]

/* v1 = m2 * v3 */
static inline void mat43_mult_vec3 (double* v1, const double* m2, const double* v3) {
    v1[0] = vec4_dot(&m2[0], v3);
    v1[1] = vec4_dot(&m2[4], v3);
    v1[2] = vec4_dot(&m2[8], v3);
}

/* m1 = m2 * m3 */
template <class T>
static inline void mat_mult_mat (T* m1, 
				 const T* m2, int m2_rows, int m2_cols, 
				 const T* m3, int m3_rows, int m3_cols)
{
    int i,j,k;
    for (i = 0; i < m2_rows; i++) {
	for (j = 0; j < m3_cols; j++) {
	    T acc = 0.0;
	    for (k = 0; k < m2_cols; k++) {
		acc += m_idx(m2,m2_cols,i,k) * m_idx(m3,m3_cols,k,j);
	    }
	    m_idx(m1,m3_cols,i,j) = acc;
	}
    }
}

#endif
