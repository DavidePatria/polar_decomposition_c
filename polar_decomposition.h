

/**
 * DOCUMENTATION
 *
 * Reference : 
 *
 *  An algorithm to compute the polar decomposition of a 3 × 3 matrix,
 *  Nicholas J. Higham, Vanni Noferini, 2016.
 *  Paper
 *: https://www.researchgate.net/publication/296638898_An_algorithm_to_compute_the_polar_decomposition_of_a_3_3_matrix/link/56e29d8c08ae03f02790a388/download
 *  Source code : https://github.com/higham/polar-decomp-3by3
 *
 * Let A be a non-singular 3×3 matrice, like the ones used in channel mixer or
 *in cameras input profiles. Such matrices define transforms between RGB and XYZ
 *spaces depending on the vector base transform. The vector base is the 3 RGB
 *primaries, defined from/to XYZ reference (absolute) space. Converting between
 * color spaces is then only a change of coordinates for the pixels color
 *vector, depending on how the primaries rotate and rescale in XYZ.
 *
 * RGB spaces conversions are therefore linear maps from old RGB to XYZ to new
 *RGB. Geometrically, linear maps can be interpreted as a combination of
 *scalings (homothety), rotations and shear mapping (transvection).
 *
 * But they also have an interesting property : 
 *
 *   For any 3×3 invertible matrice A describing a linear map, the general
 *linear map can be decomposed as a single 3D rotation around a particular 3D
 *vector.
 *
 *   That is, there is a factorization of A = Q * H, where Q is
 *   the matrice of rotation around a axis of vector H.
 *
 *
 * This is interesting for us, on the GUI side. 3×3 matrices (9 params) are not
 *intuitive to users, and the visual result of a single coefficient change is
 *hard to predict. This method allows us to reduce 9 input parameters to :
 *  * 6 : 3 angles of rotation, and the 3D coordinates of the (non-unit)
 *rotation axis vector,
 *  * 7 : 3 angles of rotation, the 3D coordinates of the unit rotation axis
 *vector, and a scaling factor for this vector.
 *
 * Usually, this is achieved by using HSL spaces, which suck because they work
 *only for bounded signals in [ 0 ; 1 ]. Also, they are not colorspaces, not
 *connected to either physics or psychology, so they are bad. Anyone saying
 * otherwise is fake news.
 *
 * The present method generalizes the HSL approach to XYZ, LMS and weird spaces,
 *with none of the drawbacks of the the cheapo lazy-ass maths-disabled HSL
 *bullshit. It's great. You should try it some time. Simply the best.
 *
 **/

// define the type for the variables that were previously TYPE as more precision was not though to be necessary
#define TYPE double

// define which ABS to use based on the chosen type
#if TYPE == double
#define ABS(n) fabs(n)
#define SQRT(n) sqrt(n)
#elif TYPE == TYPE
#define ABS(n) ABS(n)
#define SQRT(n) sqrtf(n)
#endif

#include "svd.h"
#include <complex.h>
#include <stdbool.h>
#include <stdio.h>
// is it necessary?
// #include "QR_decomp.h"

// Note : if you review this code using the original Matlab implementation,
// remember Matlab indexes arrays from 1, while C starts at 0, so every index
// needs to be shifted by -1.

// absolute value element wise for a 2x2 matrix
static inline void abs_matrix_22(TYPE matrix_vals[][2], TYPE matrix_abs[][2]) {
	int i, j;

	for (i = 0; i < 2; i++) {
		for (j = 0; j < 2; j++) {
			// printf("prima: %f\n dopo:  %f\n", matrix[i][j],
			// ABS(matrix[i][j]) );
			matrix_abs[i][j] = ABS(matrix_vals[i][j]);
			// printf("%f\n", matricella[i][j] );
		}
	}
	// printf("ok 3\n");
	printf("\n");
}

//==============================================================================

// transform from double index to single one so generic matrix multiply can be used with the result
// res is in row major order
void flatten_matrix_4x4(TYPE matrix[4][4], TYPE* flattened) {
	for (size_t i = 0; i < 4; i++) {
		for (size_t j = 0; j < 4; j++) {
			flattened[i * 4 + j] = matrix[i][j];
		}
	}
}

//==============================================================================

// multiply two matrices of any order, assuming that the multiplication is
// possible, so no checks done
void matrix_multiply(const TYPE* A, const TYPE* B, TYPE* RES, int rows_A, int cols_A, int cols_B) {
	int ii, jj, kk;
	TYPE product;

	product = 0;

	for (kk = 0; kk < rows_A; kk++) {
		for (jj = 0; jj < cols_A; jj++) {
			for (ii = 0; ii < cols_A; ii++) {
				product += A[cols_A * kk + ii] * B[cols_B * ii + jj];
			}
			RES[kk * cols_A + jj] = product;
			product               = 0;
		}
	}
}

// get the value of the greater element of a 2x2 matrix
TYPE max_val_matrix(TYPE* matrice, size_t size) {
	TYPE max;

	for (size_t i = 0; i < size; i++) {
		if (matrice[i] > max) {
			max = matrice[i];
		}
	}
	return max;
}

static inline void normalize_array(TYPE* vector, int size) {
	int ii;
	TYPE norm = 0;

	// size = sizeof(vector)/sizeof(vector[0]);

	// printf("lunghezza nella funzione: %d\n", size );
	// printf("lunghezza: %d\n", sizeof(vector));

	// computing the square of the norm of v
	for (ii = 0; ii < 4; ii++) {
		norm += vector[ii] * vector[ii];
	}
	// squared root of the norm of v
	norm = SQRT(norm);

	// printf("norm is: %f\n", norm);

	// dividing each element of v for the norm
	for (ii = 0; ii < size; ii++) {
		vector[ii] /= norm;
	}
}

void polar_decomposition(TYPE A[3][3], TYPE Q[3][3], TYPE H[3][3]) {
	// Frobenius / L2 norm of the matrice - aka we sum the squares of each
	// matrice element and take the sqrt
	const TYPE norm =
	    SQRT(A[0][0] * A[0][0] + A[0][1] * A[0][1] + A[0][2] * A[0][2] + A[1][0] * A[1][0] + A[1][1] * A[1][1] +
	         A[1][2] * A[1][2] + A[2][0] * A[2][0] + A[2][1] * A[2][1] + A[2][2] * A[2][2]);

	// Normalize the matrice A in-place, so A norm is 1
	for (size_t i = 0; i < 3; i++)
		for (size_t j = 0; j < 3; j++)
			A[i][j] /= norm;

	// Compute the conditionning of the matrice
	TYPE m, b;

	m = A[1][1] * A[2][2] - A[1][2] * A[2][1];
	b = m * m;
	m = A[1][0] * A[2][2] - A[1][2] * A[2][0];
	b += m * m;
	m = A[1][0] * A[2][1] - A[1][1] * A[2][0];
	b += m * m;

	m = A[0][0] * A[2][1] - A[0][1] * A[2][0];
	b += m * m;
	m = A[0][0] * A[2][2] - A[0][2] * A[2][0];
	b += m * m;
	m = A[0][1] * A[2][2] - A[0][2] * A[2][1];
	b += m * m;

	m = A[0][1] * A[1][2] - A[0][2] * A[1][1];
	b += m * m;
	m = A[0][0] * A[1][2] - A[0][2] * A[1][0];
	b += m * m;
	m = A[0][0] * A[1][1] - A[0][1] * A[1][0];
	b += m * m;

	b = -4.f * b + 1.f;

	TYPE d;
	bool subspa = 0;

	// copy of A
	TYPE AA[3][3] = {{A[0][0], A[0][1], A[0][2]}, {A[1][0], A[1][1], A[1][2]}, {A[2][0], A[2][1], A[2][2]}};

	size_t r = 0, c = 0;
	TYPE dd = 1.f;

	/*
  if(b - 1.f + 1e-4f) > 0.f)
  {
	we could use the quick path if perf is critical.
	It's not implemented here yet because we don't use this function for each
  pixel and for compactness of the code. See the original matlab code in case
  you need high perf.
  }
  else
  */

	// General / slow path

	// Search index (r, c) of the max element in matrice
	if (ABS(A[1][0]) > ABS(A[0][0]))
		r = 1; // c = 0
	if (ABS(A[2][0]) > ABS(A[r][c]))
		r = 2; // c = 0
	if (ABS(A[0][1]) > ABS(A[r][c]))
		r = 0, c = 1;
	if (ABS(A[1][1]) > ABS(A[r][c]))
		r = 1, c = 1;
	if (ABS(A[2][1]) > ABS(A[r][c]))
		r = 2, c = 1;
	if (ABS(A[0][2]) > ABS(A[r][c]))
		r = 0, c = 2;
	if (ABS(A[1][2]) > ABS(A[r][c]))
		r = 1, c = 2;
	if (ABS(A[2][2]) > ABS(A[r][c]))
		r = 2, c = 2;

	int k;
	// BUONO FINO A QUI

	if (r > 0) {
		// invert lines 0 and r
		TYPE temp_r[3] = {AA[r][0], AA[r][1], AA[r][2]};
		TYPE temp_0[3] = {AA[0][0], AA[0][1], AA[0][2]};

		for (size_t k = 0; k < 3; k++) {
			AA[0][k] = temp_r[k];
		};
		for (size_t k = 0; k < 3; k++)
			AA[r][k] = temp_0[k];
		dd = -dd;
	}

	if (c > 0) {
		// invert columns 0 and c
		TYPE temp_c[3] = {AA[0][c], AA[1][c], AA[2][c]};
		TYPE temp_0[3] = {AA[0][0], AA[1][0], AA[2][0]};
		for (size_t k = 0; k < 3; k++)
			AA[k][0] = temp_c[k];
		for (size_t k = 0; k < 3; k++)
			AA[k][c] = temp_0[k];
		dd = -dd;
	}

	TYPE U[3] = {AA[0][0], 0.f, 0.f};

	TYPE m0        = AA[0][1] / AA[0][0];
	TYPE m1        = AA[0][2] / AA[0][0];
	TYPE AAA[2][2] = {{AA[1][1] - AA[1][0] * m0, AA[1][2] - AA[1][0] * m1},
	                  {AA[2][1] - AA[2][0] * m0, AA[2][2] - AA[2][0] * m1}};

	r = 0, c = 0;
	if (ABS(AA[1][0]) > ABS(AA[0][0]))
		r = 1; // c = 0
	if (ABS(AA[0][1]) > ABS(AA[r][c]))
		r = 0, c = 1;
	if (ABS(AA[1][1]) > ABS(AA[r][c]))
		r = 1, c = 1;

	if (r == 1)
		dd = -dd;
	if (c > 0)
		dd = -dd;

	U[1] = AA[r][c];

	// fixed from U(2). needs check
	if (U[1] == 0)
		U[2] = 0;
	else
		U[2] = AA[2 - r][2 - c] - AA[r][2 - c] * AA[2 - r][c] / U[1];

	d  = dd;
	dd = dd * U[0] * U[1] * U[2];

	if (U[0] < 0)
		d = -d;
	if (U[1] < 0)
		d = -d;
	if (U[2] < 0)
		d = -d;

	TYPE AU = ABS(U[2]);

	TYPE nit;

	if (AU > 6.607e-8f) {
		nit = 16.8f + 2.f * log10f(AU);
		nit = ceilf(15.f / nit);
	} else {
		subspa = true;
	}

	if (d == 0)
		d = 1.f;

	dd = 8.f * d * dd;

	TYPE t = A[0][0] + A[1][1] + A[2][2];

	TYPE B[4][4] = {{t, A[1][2] - A[2][1], A[2][0] - A[0][2], A[0][1] - A[1][0]},
	                {0.f, 2.f * A[0][0] - t, A[0][1] + A[1][0], A[0][2] + A[2][0]},
	                {0.f, 0.f, 2.f * A[1][1] - t, A[1][2] + A[2][1]},
	                {0.f, 0.f, 0.f, 2.f * A[2][2] - t}};

	for (size_t i = 0; i < 4; ++i)
		for (size_t j = 0; j < 4; ++j)
			B[i][j] /= d;

	B[1][0] = B[0][1];
	B[2][0] = B[0][2];
	B[3][0] = B[0][3];

	B[2][1] = B[1][2];
	B[3][1] = B[1][3];
	B[3][2] = B[2][3];

	// Find largest eigenvalue
	TYPE x;

	// ISSUE: possibile cause of the problem
	if (b >= -0.3332f) {
		// Use analytic formula if matrice is well conditioned
		double complex Delta0 = 1.f + 3. * b;
		double complex Delta1 = -1. + (27. / 16.) * dd * dd + 9. * b;
		double complex phi    = (Delta1 / Delta0) / csqrt(Delta0);
		double complex SS     = (4. / 3.) * (1. + ccosf(cacosf(phi) / 3.) * csqrt(Delta0));
		double complex S      = csqrt(SS) / 2.;

		x = (TYPE)(creal(S) + 0.5 * sqrt(fmaxf(0., creal(-SS + 4. + dd / S))));
	} else {
		// Use Newton if matrice is ill conditioned
		// We use double precision temporarily because the solution can
		// degenerate faster in single precision
		double x_temp = sqrt(3.);
		double xold   = 3;
		while ((xold - x_temp) > 1e-12) {
			xold       = x_temp;
			double px  = x_temp * (x_temp * (x_temp * x_temp - 2.) - dd) + b;
			double dpx = x_temp * (4. * x_temp * x_temp - 4.) - dd;
			x_temp     = x_temp - px / dpx;
		}
		x = (TYPE)x_temp;
	}
	// X OUT OF HERE SEEMS TO BE DIFFERENT FROM THE MATLAB ONE

	// Again, don't do the quick path
	TYPE BB[4][4];

	for (size_t i = 0; i < 4; ++i)
		for (size_t j = 0; j < 4; ++j) {
			BB[i][j] -= B[i][j];
			if (i == j)
				BB[i][j] += x; // add x on the diagonal
		}

	// different from matlab because this is an index and matlab starts from 1
	size_t p[4] = {0, 1, 2, 3};

	TYPE L[4][4] = {{1.f, 0.f, 0.f, 0.f}, {0.f, 1.f, 0.f, 0.f}, {0.f, 0.f, 1.f, 0.f}, {0.f, 0.f, 0.f, 1.f}};

	// not sure about this one initializing
	// double D[4][4] = {{0.f}};
	TYPE D[4][4] = {{0.f, 0.f, 0.f, 0.f}, {0.f, 0.f, 0.f, 0.f}, {0.f, 0.f, 0.f, 0.f}, {0.f, 0.f, 0.f, 0.f}};

	// First step
	r = 3;
	if (BB[3][3] < BB[2][2])
		r = 2;
	if (BB[r][r] < BB[1][1])
		r = 1;
	if (BB[r][r] > BB[0][0]) {
		// ASWME: not entirely sure about this one, in octave the result is
		// different p([1 r(1)]) = [r(1) 1]
		p[0] = r;
		p[r] = 0;

		// BB = BB(p,p);
		for (size_t i = 0; i < 4; ++i)
			for (size_t j = 0; j < 4; ++j)
				BB[i][j] = BB[p[i]][p[j]];
	}

	// according to what the matlab code does D(1) = B(1,1) assigns the first
	// element of D, despite D being 4x4
	D[0][0] = BB[0][0];

	L[1][0] = BB[1][0] / D[0][0];
	L[2][0] = BB[2][0] / D[0][0];
	L[3][0] = BB[3][0] / D[0][0];

	BB[1][1] -= L[1][0] * BB[0][1];
	BB[2][1] -= L[1][0] * BB[0][2];
	BB[1][2] = BB[2][1];

	BB[3][1] -= L[1][0] * BB[0][3];
	BB[1][3] = BB[3][1];
	BB[2][2] -= L[2][0] * BB[0][2];

	BB[3][2] -= L[2][1] * BB[1][3];
	BB[2][3] = D[3][2];
	BB[3][3] -= L[3][0] * BB[0][3];

	// Second step
	r = 2;

	if (BB[2][2] < BB[1][1])
		r = 1;

	if (BB[r][r] > BB[1][1]) {
		p[1] = r;
		p[r] = 2;

		// BB([2 r(1)],:) = BB([r(1) 2],:);
		//  another scope so the variables can be reused outside of it
		// get row number 1 (which is 2 in matlab)
		// TYPE temp_q[4] = {BB[1][0], BB[1][1],BB[1][2],BB[1][3]};
		TYPE temp[4];
		// now row 1 can be overwirtten
		for (size_t i = 0; i < 4; i++) {
			temp[i]  = BB[1][i];
			BB[1][i] = BB[r][i];
			BB[r][i] = temp[i];
		}
		// get column number 1 (which is 2 in matlab)
		for (size_t i = 0; i < 4; i++) {
			temp[i]  = BB[i][1];
			BB[i][1] = BB[i][r];
			BB[i][r] = temp[i];
		}
		// now swap the rows of L in the same way
		for (size_t i = 0; i < 3; i++) {
			temp[i] = L[1][i];
			L[1][i] = L[r][i];
			L[r][i] = temp[i];
		}
		// now swap the columns of L in the same way
		// get column number 1 (which is 2 in matlab)
		for (size_t i = 0; i < 3; i++) {
			temp[i] = L[i][1];
			L[i][1] = L[i][r];
			L[i][r] = temp[i];
		}
	}

	D[1][1] = BB[1][1];

	L[2][1] = BB[2][1] / D[1][1];
	L[3][1] = BB[3][1] / D[1][1];

	D[2][2] = BB[2][2] - L[2][1] * BB[1][2];
	D[3][2] = BB[3][2] - L[2][1] * BB[1][3];
	D[2][3] = D[3][2];
	D[3][3] = BB[3][3] - L[3][1] * BB[1][3];

	TYPE DD = D[2][2] * D[3][3] - D[2][3] * D[2][3];

	// skip DD == 0

	TYPE ID[2][2] = {{D[3][3], D[2][3]}, {D[2][3], D[2][2]}};

	// SKIP SUBSPA == 1
	// going directly for else, subspa = false

	TYPE v[4] = {L[1][0] * L[3][1] + L[2][0] * L[3][2] - L[1][0] * L[3][2] * L[2][1] - L[3][0],
	             L[3][2] * L[2][1] - L[3][1], -L[3][2], 1};
	// this would already be defined if the other case for subspa were done
	TYPE IL[4][4] = {
	    {1, 0, 0, 0}, {-L[1][0], 1, 0, 0}, {L[1][0] * L[2][1] - L[2][0], -L[2][1], 1, 0}, {v[0], v[1], v[2], v[3]}};

	// normalize array
	normalize_array(v, 4);

	for (size_t i = 0; i < nit; i++) {
		// flat array for the result of product
		TYPE RES[4];
		// matrix_multiply(IL_flat, v, RES, 4, 4, 1);
		matrix_multiply(*IL, v, RES, 4, 4, 1);
		// now copy the result of multiplication into the array
		for (size_t i = 0; i < 4; i++) {
			v[i] = RES[i];
		}

		v[0] /= D[0][0];
		v[1] /= D[1][1];
		// first row of ID by v[2:-1]
		v[2] = (ID[0][0] * v[2] + ID[0][1] * v[3]) / DD;
		v[3] = (ID[1][0] * v[2] + ID[1][1] * v[3]) / DD;

		// TYPE RES[4] = {};
		// this should work given that v can be "transposed" by simply changing sizes in the function
		matrix_multiply(v, *IL, RES, 1, 4, 4);
		// assign the result of multiplication to v
		for (size_t i = 0; i < 4; i++) {
			v[i] = RES[i];
		}
		// this does not matter, right?
		// v = v';
		normalize_array(v, 4);
	}

	// not 100% sure about this one, though is seems the counterintuitive behaviour of the matlab code
	// TODO: check
	v[0] = v[p[0]];
	v[1] = v[p[1]];
	v[0] = v[p[2]];
	v[0] = v[p[3]];

	TYPE v11, v22, v33, v12, v03, v13, v02, v01, v23;

	v01 = 2 * v[0] * v[1];
	v02 = 2 * v[0] * v[2];
	v03 = 2 * v[0] * v[3];

	v11 = 2 * v[1] * v[1];
	v22 = 2 * v[2] * v[2];
	v33 = 2 * v[3] * v[3];

	v12 = 2 * v[1] * v[2];
	v13 = 2 * v[1] * v[3];
	v23 = 2 * v[2] * v[3];

	Q[0][0] = 1 - v22 - v33;
	Q[0][1] = v12 + v03;
	Q[0][2] = v13 - v02;

	Q[1][0] = v12 - v03;
	Q[1][1] = 1 - v11 - v33;
	Q[1][2] = v01 + v23;

	Q[2][0] = v02 + v13;
	Q[2][1] = v23 - v01;
	Q[2][2] = 1 - v11 - v22;

	if (d == -1) {
		// invert th sign of every element of Q
		for (size_t i = 0; i < 3; i++) {
			for (size_t j = 0; j < 3; j++) {
				Q[i][j] *= -1;
			}
		}
	};

	// since Q' is only needed for this product it can be stored into a temporary to be used
	TYPE Q_t[3][3];

	for (size_t i = 0; i < 3; i++) {
		for (size_t j = 0; j < 3; j++) {
			Q_t[j][i] = Q[i][j];
		}
	}

	// tostore the flattened matrix
	// TYPE Q_t_flat[9];
	// flatten_matrix_3x3(temp, Q_t_flat);
	// to be unflattened to be assigned to the result
	// TYPE H_flat[9];

	// TYPE A_flat[9];
	// flatten_matrix_3x3(A, A_flat);

	matrix_multiply(*Q_t, *A, *H, 3, 3, 3);
	// multiply by the norm again
	for (size_t jj = 0; jj < 9; jj++) {
		(*H)[jj] *= norm;
		// also A was normalized at the beginning
		(*A)[jj] *= norm;
	}

	// unflatten_matrix_3x3(H, H_flat);
	// unflatten_matrix_3x3(H, H_flat);
	// unflatten_matrix_3x3(A, A_flat);
};
