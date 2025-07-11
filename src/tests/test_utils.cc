/*!
  \file   test_utils.cc
  \author Simon Pfreundschuh <simonpf@chalmers.se>
  \date   Sun May  3 21:56:19 2015
*/

#include "test_utils.h"

#include <matpack.h>

#include <cmath>

using std::abs;

//! Add noise to vector.
/*!

  \param[in,out] v The vector to add the noise to.
  \param[in] scale Range for the generated noise given by [0,scale].
*/
void add_noise(VectorView v, Numeric scale) {
  Rand<Numeric> rand(0, scale);

  for (Size i = 0; i < v.size(); i++) {
    v[i] += rand();
  }
}

//! Fill matrix with random values.
/*!

Fills the given matrix with random values of type Numeric in the range
[0, range], if positive is set to true, or [-range, range], if positive
is set to false.

  \param[in,out] A The matrix to be filled.
  \param[in] range The range of the values to draw the values from.
  \param[in] positive If true the matrix is filled with values from the interval
                  [0,range], otherwise the values are taken from the interval
                  [-range, range].
*/
void random_fill_matrix(StridedMatrixView A, Numeric range, bool positive) {
  Index m = A.nrows();
  Index n = A.ncols();

  Rand<Numeric> rand(positive ? 0 : -range, range);

  for (Index i = 0; i < m; i++) {
    for (Index j = 0; j < n; j++) {
      A[i, j] = (Numeric)rand();
    }
  }
}
void random_fill_matrix(StridedComplexMatrixView A, Numeric range, bool positive) {
  Index m = A.nrows();
  Index n = A.ncols();

  Rand<Numeric> rand(positive ? 0 : -range, range);

  for (Index i = 0; i < m; i++) {
    for (Index j = 0; j < n; j++) {
      A[i, j] = Complex((Numeric)rand(), (Numeric)rand());
    }
  }
}

//! Generate random sparse matrix
/*!

  Fills a sparse m-times-n matrix A with max{m,n} random values at random
  positions.

  \param[out] A The matrix to be filled.
  \param[in] range The range from which values are genereated. If positive = true
  the values are generated from the range [0,range], otherwise from the range
  [-range, range]
  \param[in] positive See above.
*/
void random_fill_matrix(Sparse& A, Numeric range, bool positive) {
  Index m = A.nrows();
  Index n = A.ncols();

  Rand<Numeric> random_number(positive ? 0 : -range, range);

  Index nelem = std::max(m, n);

  for (Index i = 0; i < nelem; i++) {
    Index m1, n1;
    m1 = rand() % m;
    n1 = rand() % n;

    A.rw(m1, n1) = random_number();
  }
}

//! Generate identical, random sparse and dense matrices.
/*!

  Fills a dense and a sparse m-times-n matrix A with max{m,n} random values at
  random positions.

  \param[out] A The dense matrix to be filled.
  \param[out] B The sparse matrix to be filled.
  \param[in] range The range from which values are genereated. If positive = true
  the values are generated from the range [0,range], otherwise from the range
  [-range, range]
  \param[in] positive See above.
*/
void random_fill_matrix(Matrix& A, Sparse& B, Numeric range, bool positive) {
  Index m = A.nrows();
  Index n = A.ncols();

  assert(B.nrows() == m);
  assert(B.ncols() == n);

  Rand<Numeric> random_number(positive ? 0 : -range, range);

  Index nelem = std::max(m, n);

  for (Index i = 0; i < nelem; i++) {
    Index m1, n1;
    m1 = rand() % m;
    n1 = rand() % n;

    A[m1, n1]    = random_number();
    B.rw(m1, n1) = A[m1, n1];
  }
}

//! Generate random, symmetric matrix.
/*!

  \param[out] A The matrix to be filled.
  \param[in] range The range from which values are genereated. If positive = true
  the values are generated from the range [0,range], otherwise from the range
  [-range, range]
  \param[in] positive See above.
*/
void random_fill_matrix_symmetric(StridedMatrixView A, Numeric range, bool positive) {
  random_fill_matrix(A, range, positive);
  Matrix M(A);
  A += transpose(M);
}
void random_fill_matrix_symmetric(StridedComplexMatrixView A,
                                  Numeric range,
                                  bool positive) {
  random_fill_matrix(A, range, positive);
  ComplexMatrix M(A);
  A += transpose(M);
}

//! Generate random, positive definite matrix.
/*!

  Generate a random, positive definite matrix by generating
  a positive semi-definite matrix and adding the identity matrix.

  \param[out] A The random, positive definite matrix.
  \param[in] range The range from which the random values are picked. If
  positive == true, the values are taken from the range [0,range], otherwise
  the are taken from the range [-range, range].
  \param[in] positive See above.
*/
void random_fill_matrix_pos_def(StridedMatrixView A, Numeric range, bool positive) {
  Index n = A.ncols();

  // Ensure that A is square.
  assert(A.ncols() == A.nrows());

  // Generate random, pos. semi-def. Matrix
  random_fill_matrix(A, range, positive);
  Matrix M(A);
  mult(A, M, transpose(M));

  // Add identity matrix.
  for (Index i = 0; i < n; i++) {
    A[i, i] += 1.0;
  }
}

//! Generate random, positive semi-definite matrix.
/*!

  Generate a random, positive semi-definite matrix by
  randomly generating a matrix and multiplying it by its transpose.

  \param[out] A The random, positive semi-definite matrix.
  \param[in] range The range from which the random values are picked. If
  positive == true, the values are taken from the range [0,range], otherwise
  the are taken from the range [-range, range].
  \param[in] positive See above.
*/
void random_fill_matrix_pos_semi_def(StridedMatrixView A,
                                     Numeric range,
                                     bool positive) {
  random_fill_matrix(A, range, positive);
  Matrix M(A);
  mult(A, M, transpose(M));
}

//! Fill vector with random values.
/*!

  Fills the given vector with random values of type Numeric drawn
  from the range [0, range], if positive is set to true, or from the
  range [-range, range], if positive == false.

  \param[in,out] v The vector to be filled.
  \param[in] range The range from which the values are taken.
  \param[in] positive If true, the values are taken from the interval [0, range],
                      otherwise from the range [-range, range].
*/
void random_fill_vector(VectorView v, Numeric range, bool positive) {
  Index n = v.size();

  Rand<Numeric> rand(positive ? 0 : -range, range);

  for (Index i = 0; i < n; i++) {
    v[i] = rand();
  }
}

//! Pick random random submatrix of size m times n.
/*!
  Randomly chooses a submatrix of the given matrix A and returns the
  corresponding StridedMatrixView.

  \param[in] A The matrix to choose the submatrix from.
  \param[in] m Number of rows of the submatrix.
  \param[in] n Number of columns of the submatrix.

  \return StridedConstMatrixView corresponding to a randomly chosen m-by-n submatrix.
*/
StridedMatrixView random_submatrix(StridedMatrixView A, int m, int n) {
  Index m0(A.nrows()), n0(A.ncols());
  assert((0 <= m) && (m <= m0));
  assert((0 <= n) && (m <= n0));

  Rand<Index> rand_m(0, (m0 - m - 1)), rand_n(0, (n0 - n - 1));
  Index m1, n1;
  m1 = rand_m();
  n1 = rand_n();

  Range r1(m1, m), r2(n1, n);
  return A[r1, r2];
}

//! Generate random sub-range of the range [0, n-1].
/*!
  Generate random Range object such that 0 <= extent <= n
  and 0 <= start < n - extent.

  \param n The range [0, n-1] to pick the sub-range from.
  \return Random sub-range.
*/
Range random_range(Index n) {
  Rand<Index> extent(1, n);
  Index e = extent();

  Index s = 0;
  if (0 <= (n - e - 1)) {
    Rand<Index> start(0, n - e - 1);
    s = start();
  }
  return Range(s, e);
}

//! Maximum element-wise error of two vectors.
/*!
  If relative is true, the maximum element-wise error is computed.
  Otherwise the absolute error is computed.

  \param[in] v1 The first vector.
  \param[in] v2 The reference vector used to normalize the relative error.
  \param[in] relative If true the relative error is computed, otherwise the absolute
                  error is computed.

  \return The maximum relative or absolute element-wise error.
*/
Numeric get_maximum_error(ConstVectorView v1,
                          ConstVectorView v2,
                          bool relative) {
  Index n = std::min(v1.size(), v2.size());

  Numeric max = 0.0, err = 0.0;

  for (Index i = 0; i < n; i++) {
    err = 0.0;

    if (relative) {
      if (v2[i] != 0.0) {
        err = abs((v2[i] - v1[i]) / v2[i]);
      }

    } else {
      err = abs(v2[i] - v2[i]);
    }

    if (err > max) {
      max = err;
    }
  }
  return err;
}

//! Maximum element-wise error of two matrices.
/*!
  If relative is true, the maximum element-wise error is computed.
  Otherwise the absolute error is computed.

  \param[in] A1 The first matrix.
  \param[in] A2 The reference matrix used to normalize the relative error.
  \param[in] relative If true the relative error is computed, otherwise the absolute
                  error is computed.

  \return The maximum relative or absolute element-wise error.
*/
Numeric get_maximum_error(StridedConstMatrixView A1,
                          StridedConstMatrixView A2,
                          bool relative) {
  Index m = std::min(A1.nrows(), A2.nrows());
  Index n = std::min(A1.ncols(), A2.ncols());

  Numeric max = 0.0, err = 0.0;

  for (Index i = 0; i < m; i++) {
    for (Index j = 0; j < n; j++) {
      err = 0.0;

      if (relative) {
        if (A2[i, j] != 0.0) {
          err = abs((A2[i, j] - A1[i, j]) / A2[i, j]);
        }

      } else {
        err = A2[i, j] - A1[i, j];
      }

      if (err > max) {
        max = err;
      }
    }
  }

  return max;
}
Numeric get_maximum_error(StridedConstComplexMatrixView A1,
                          StridedConstComplexMatrixView A2,
                          bool relative) {
  Index m = std::min(A1.nrows(), A2.nrows());
  Index n = std::min(A1.ncols(), A2.ncols());

  Numeric max = 0.0, err = 0.0;

  for (Index i = 0; i < m; i++) {
    for (Index j = 0; j < n; j++) {
      err = 0.0;

      if (relative) {
        if (A2[i, j].real() != 0.0 && A2[i, j].imag() != 0.0) {
          err = abs((A2[i, j] - A1[i, j]) / A2[i, j]);
        }

      } else {
        err = abs(A2[i, j] - A1[i, j]);
      }

      if (err > max) {
        max = err;
      }
    }
  }

  return max;
}
