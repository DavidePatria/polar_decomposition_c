

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
// an if so the type can be choosen with a gcc flag but would fallback to this if not defined
#ifndef TYPE
#define TYPE double
#endif

// define which ABS to use based on the chosen type
#if TYPE == double
#define ABS(n) fabs(n)
#define SQRT(n) sqrt(n)
#elif TYPE == TYPE
#define ABS(n) fabs(n)
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

void swap_rows(TYPE* m, size_t rows, size_t cols, size_t row0, size_t row1);
void swap_cols(TYPE* m, size_t rows, size_t cols, size_t col0, size_t col1);
void abs_matrix(TYPE* matrix_abs, size_t num_el);
// multiply two matrices of any order, assuming that the multiplication is
// possible, so no checks done
void matrix_multiply(const TYPE* A, const TYPE* B, TYPE* RES, int rows_A, int cols_A, int cols_B);
// get the value of the greater element of a 2x2 matrix
TYPE max_val_matrix(TYPE* matrice, size_t size);
static inline void normalize_array(TYPE* vector, int size);
void stampa_matrice(TYPE* m, size_t rows, size_t cols);
void compute_null_space(TYPE* nullspace, const TYPE a, const TYPE b, const TYPE c);
void copy(TYPE* dest, TYPE* source, size_t num_el);

// shamelessly copied function
void orthonormalize_v_with_qr(TYPE* v0, TYPE* v1, const TYPE v00, const TYPE v10, const TYPE v01, const TYPE v11);
void multiply_il_v(TYPE* RES, const TYPE IL01, const TYPE IL02, const TYPE IL03, const TYPE IL12, const TYPE IL13,
                   const TYPE* v);
void multiply_id_v(TYPE* RES, const TYPE ID00, const TYPE ID11, const TYPE* ID, const TYPE* v);
void multiply_v_il(TYPE* RES, const TYPE* v, const TYPE IL01, const TYPE IL02, const TYPE IL03, const TYPE IL12,
                   const TYPE IL13);
void orthonormalize_v_with_qr_single(TYPE* v0, TYPE* v1);
void multiply_minus_v_d(TYPE* RES, TYPE* v, TYPE* D);
TYPE vector_scalar(const TYPE* a, const TYPE* b, const size_t num_el);

//==============================================================================
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

	b = -4.0 * b + 1.0;

	TYPE d;
	bool subspa = 0;

	// copy of A
	TYPE AA[3][3] = {{A[0][0], A[0][1], A[0][2]}, {A[1][0], A[1][1], A[1][2]}, {A[2][0], A[2][1], A[2][2]}};

	size_t r = 0, c = 0;
	TYPE dd = 1.0;

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

	if (r > 0) {
		// invert lines 0 and r
		swap_rows(*AA, 3, 3, 0, r);
		dd = -dd;
	}

	if (c > 0) {
		swap_cols(*AA, 3, 3, 0, c);
		dd = -dd;
	}

	TYPE U[3] = {AA[0][0], 0.0, 0.0};

	TYPE m0        = AA[0][1] / AA[0][0];
	TYPE m1        = AA[0][2] / AA[0][0];
	TYPE AAA[2][2] = {{AA[1][1] - AA[1][0] * m0, AA[1][2] - AA[1][0] * m1},
	                  {AA[2][1] - AA[2][0] * m0, AA[2][2] - AA[2][0] * m1}};

	r = 0, c = 0;
	if (ABS(AAA[1][0]) > ABS(AAA[0][0]))
		r = 1; // c = 0
	if (ABS(AAA[0][1]) > ABS(AAA[r][c]))
		r = 0, c = 1;
	if (ABS(AAA[1][1]) > ABS(AAA[r][c]))
		r = 1, c = 1;

	if (r == 1)
		dd = -dd;
	if (c > 0)
		dd = -dd;

	U[1] = AAA[r][c];

	// fixed from U(2). needs check
	if (U[1] == 0)
		U[2] = 0;
	else
		// l'espressione qui sotto fa zero e non dovrebbe, perché dopo quando AU viene confrontata è falso ma dovrebbe
		// essere vero
		U[2] = AAA[1 - r][1 - c] - AAA[r][1 - c] * AAA[1 - r][c] / U[1];

	d  = dd;
	dd = dd * U[0] * U[1] * U[2];

	if (U[0] <= 0)
		d = -d;
	if (U[1] <= 0)
		d = -d;
	if (U[2] <= 0)
		d = -d;

	TYPE AU = ABS(U[1]);

	TYPE nit;

	if (AU > 6.607e-8) {
		nit = 16.80 + 2.0 * log10(AU);
		nit = ceil(15.0 / nit);
	} else {
		subspa = true;
	}

	if (d == 0)
		d = 1.0;

	dd = 8.0 * d * dd;

	TYPE t = A[0][0] + A[1][1] + A[2][2];

	TYPE B[4][4] = {{t, A[1][2] - A[2][1], A[2][0] - A[0][2], A[0][1] - A[1][0]},
	                {0.0, 2.0 * A[0][0] - t, A[0][1] + A[1][0], A[0][2] + A[2][0]},
	                {0.0, 0.0, 2.0 * A[1][1] - t, A[1][2] + A[2][1]},
	                {0.0, 0.0, 0.0, 2.0 * A[2][2] - t}};

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
		TYPE complex Delta0 = 1.f + 3. * b;
		TYPE complex Delta1 = -1. + (27. / 16.) * dd * dd + 9. * b;
		TYPE complex phi    = (Delta1 / Delta0) / csqrt(Delta0);
		TYPE complex SS     = (4. / 3.) * (1. + ccos(cacos(phi) / 3.) * csqrt(Delta0));
		TYPE complex S      = csqrt(SS) / 2.;

		x = (TYPE)(creal(S) + 0.5 * sqrt(fmax(0., creal(-SS + 4. + dd / S))));
	} else {
		// Use Newton if matrice is ill conditioned
		// We use double precision temporarily because the solution can
		// degenerate faster in single precision
		TYPE x_temp = sqrt(3.);
		TYPE xold   = 3;
		while ((xold - x_temp) > 1e-12) {
			xold     = x_temp;
			TYPE px  = x_temp * (x_temp * (x_temp * x_temp - 2.) - dd) + b;
			TYPE dpx = x_temp * (4. * x_temp * x_temp - 4.) - dd;
			x_temp   = x_temp - px / dpx;
		}
		x = (TYPE)x_temp;
	}

	// Again, don't do the quick path
	TYPE BB[4][4];

	for (size_t i = 0; i < 4; ++i) {
		for (size_t j = 0; j < 4; ++j) {
			BB[i][j] -= B[i][j];
			if (i == j)
				BB[i][j] += x; // add x on the diagonal
		}
	}

	// different from matlab because this is an index and matlab starts from 1
	size_t p[4]  = {0, 1, 2, 3};
	TYPE L[4][4] = {{1.0, 0.0, 0.0, 0.0}, {0.0, 1.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 1.0}};
	TYPE D[4][4] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};

	// First step
	r = 3;
	if (BB[3][3] < BB[2][2])
		r = 2;
	if (BB[r][r] < BB[1][1])
		r = 1;
	if (BB[r][r] > BB[0][0]) {
		p[0] = r;
		p[r] = 0;

		TYPE temp[4][4] = {0.0};

		// for (int i = 0; i < 4; i++)
		// 	for (int j = 0; j < 4; j++) {
		// 		temp[i][j] = BB[i][j];
		// 	}

		copy(*temp, *BB, 16);

		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; j++) {
				BB[i][j] = temp[p[i]][p[j]];
			}
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

	BB[3][2] -= L[2][0] * BB[0][3]; // wrong result
	BB[2][3] = BB[3][2];
	BB[3][3] -= L[3][0] * BB[0][3];

	// Second step
	r = 2;

	if (BB[2][2] < BB[1][1])
		r = 1;

	// does not work for simple case one
	if (BB[r][r] > BB[0][0]) {
		p[1] = r;
		p[r] = 2;

		// BB([2 r(1)],:) = BB([r(1) 2],:);
		swap_rows(*BB, 4, 4, 1, r);
		// get column number 1 (which is 2 in matlab)
		swap_cols(*BB, 4, 4, 1, r);
		// now swap the rows of L in the same way
		swap_rows(*L, 4, 4, 1, r);
		// now swap the columns of L in the same way
		// get column number 1 (which is 2 in matlab)
		swap_cols(*L, 4, 4, 1, r);
	}

	D[1][1] = BB[1][1];

	L[2][1] = BB[2][1] / D[1][1];
	L[3][1] = BB[3][1] / D[1][1];

	D[2][2] = BB[2][2] - L[2][1] * BB[1][2];
	D[3][2] = BB[3][2] - L[2][1] * BB[1][3];
	D[2][3] = D[3][2];
	D[3][3] = BB[3][3] - L[3][1] * BB[1][3];

	TYPE DD = D[2][2] * D[3][3] - D[2][3] * D[2][3];
	TYPE v[4];

	// TODO: test if case
	// it is still unclear how to manage the definition of v according to all the if and else
	if (DD == 0) {
		printf("DD is zero\n");
		// DD is symmetric, instead of doing max(abs(D)) this is enough
		const bool all_zero = (D[2][2] == 0 && D[3][3] == 0 && D[3][2] == 0);
		if (all_zero == 0) {
			// TYPE v[4] = {L[1][0] * L[3][1] - L[3][0], -L[3][1], 0, 1};
			v[0] = L[1][0] * L[3][1] - L[3][0];
			v[1] = L[3][1];
			v[2] = 0;
			v[2] = 1;
		} else {
			TYPE nullspace[3];
			compute_null_space(nullspace, D[2][2], D[3][2], D[3][3]);
			v[0] = nullspace[0] * L[1][0] * L[2][1] + nullspace[1] * L[1][0] * L[3][1] - L[3][0];
			v[1] = -nullspace[0] * L[2][1] - nullspace[1] * L[3][1];
			v[2] = nullspace[0];
			v[3] = nullspace[1];
		}
	} else {

		TYPE ID[2][2] = {{D[3][3], -D[2][3]}, {-D[2][3], D[2][2]}};

		// this inizialitaion values are the ones for inside suabspa==false, but are set here to avoid the assle of
		// reassigning them after one by one

		// shamelessly copied from the cpp implementation, doing it on my own was too airy
		if (subspa) {

			const TYPE IL01 = -L[1][0];                    //(0,1);
			const TYPE IL02 = L[1][0] * L[2][1] - L[2][1]; // L(0,1) * L(1,2) - L(0,2);
			const TYPE IL12 = L[2][1];                     //-L(1,2);
			const TYPE IL03 = L[1][0] * L[3][1] - L[3][0]; // L(0,1) * L(1,3) - L(0,3);
			const TYPE IL13 = -L[3][1];                    //-L(1,3);

			const TYPE v00 = IL02;
			const TYPE v10 = IL03;
			const TYPE v01 = IL12;
			const TYPE v11 = IL13;

			const TYPE ID00 = 1 / D[0][0];
			const TYPE ID11 = 1 / D[1][1];

			TYPE v0[4], v1[4];

			orthonormalize_v_with_qr(v0, v1, v00, v10, v01, v11);

			for (int i = 0; i < 2; ++i) {
				multiply_il_v(v0, IL01, IL02, IL03, IL12, IL13, v0);
				multiply_il_v(v0, IL01, IL02, IL03, IL12, IL13, v1);

				multiply_id_v(v0, ID00, ID11, *ID, v0);
				multiply_id_v(v1, ID00, ID11, *ID, v1);

				multiply_v_il(v0, v0, IL01, IL02, IL03, IL12, IL13);
				multiply_v_il(v1, v1, IL01, IL02, IL03, IL12, IL13);
			}

			orthonormalize_v_with_qr_single(v0, v1);

			// TODO: conclude the focking work

			TYPE v0_temp[4];
			multiply_v_il(v0_temp, v0, L[1][0], L[2][0], L[3][0], L[2][1], L[3][1]);
			TYPE v1_temp[4];
			multiply_v_il(v1_temp, v1, L[1][0], L[2][0], L[3][0], L[2][1], L[3][1]);

			TYPE H0[4];
			multiply_minus_v_d(H0, v0_temp, *D);
			TYPE H1[4];
			multiply_minus_v_d(H1, v1_temp, *D);

			const TYPE H00 = vector_scalar(H0, v0_temp, 4);
			const TYPE H10 = vector_scalar(H0, v1_temp, 4);
			const TYPE H01 = vector_scalar(H1, v0_temp, 4);
			const TYPE H11 = vector_scalar(H1, v1_temp, 4);

			if (ABS(H10) < (TYPE)(1.0e-15)) {
				if (H00 > H10) {
					copy(v, v0, 4);
				} else {
					copy(v, v1, 4);
				}
			} else {
				const TYPE r = (H00 - H11) / (2 * H10);
				const int s  = (H10 < 0 ? -1 : 1);
				const TYPE f = r + s * SQRT(1 + r * r);
				v[0]         = v0[0] * f + v1[0];
				v[1]         = v0[1] * f + v1[1];
				v[2]         = v0[2] * f + v1[2];
				v[3]         = v0[3] * f + v1[3];
			}
		} else {
			// this inizialitaion values are the ones for inside suabspa==false, but are set here to avoid the assle of
			// reassigning them after one by one TYPE IL[4][4] = { {1, 0, 0, 0}, {-L[1][0], 1, 0, 0}, {L[1][0] * L[2][1]
			// - L[2][0], -L[2][1], 1, 0}, {v[0], v[1], v[2], v[3]}};

			// TYPE v[4] = {L[1][0] * L[3][1] + L[2][0] * L[3][2] - L[1][0] * L[3][2] * L[2][1] - L[3][0],
			//              L[3][2] * L[2][1] - L[3][1], -L[3][2], 1};
			v[0] = L[1][0] * L[3][1] + L[2][0] * L[3][2] - L[1][0] * L[3][2] * L[2][1] - L[3][0];
			v[1] = L[3][2] * L[2][1] - L[3][1];
			v[2] = -L[3][2];
			v[3] = 1;

			TYPE IL[4][4] = {{1, 0, 0, 0},
			                 {-L[1][0], 1, 0, 0},
			                 {L[1][0] * L[2][1] - L[2][0], -L[2][1], 1, 0},
			                 {v[0], v[1], v[2], v[3]}};
			// this would already be defined if the other case for subspa were done

			// normalize array
			normalize_array(v, 4);

			for (size_t i = 0; i < nit; i++) {
				// flat array for the result of product
				TYPE RES[4];
				// matrix_multiply(IL_flat, v, RES, 4, 4, 1);
				matrix_multiply(*IL, v, RES, 4, 4, 1);
				// now copy the result of multiplication into the array
				copy(v, RES, 4);

				v[0] /= D[0][0];
				v[1] /= D[1][1];
				v[2] = (ID[0][0] * RES[2] + ID[0][1] * RES[3]) / DD;
				v[3] = (ID[1][0] * RES[2] + ID[1][1] * RES[3]) / DD;

				// TYPE RES[4] = {};
				// this should work given that v can be "transposed" by simply changing sizes in the function
				matrix_multiply(v, *IL, RES, 1, 4, 4);
				// assign the result of multiplication to v
				copy(v, RES, 4);
				normalize_array(v, 4);
			}
		} // end else (that is susbspa==false)
	}     // end else DD == 0

	//============================================================================
	// copy the vector so it does not change during permutation
	TYPE* temp = (TYPE*)calloc(1, 4 * sizeof(TYPE));
	// for (size_t i = 0; i < 4; i++) {
	// 	temp[i] = v[i];
	// }
	copy(temp, v, 4);

	v[0] = temp[p[0]];
	v[1] = temp[p[1]];
	v[2] = temp[p[2]];
	v[3] = temp[p[3]];

	free(temp);
	//============================================================================

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
		// invert the sign of every element of Q
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

	matrix_multiply(*Q_t, *A, *H, 3, 3, 3);
	// multiply by the norm again
	for (size_t jj = 0; jj < 9; jj++) {
		(*H)[jj] *= norm;
		// also A was normalized at the beginning
		(*A)[jj] *= norm;
	}
};

//==============================================================================

void abs_matrix(TYPE* matrix_abs, size_t num_el) {
	int i, j;

	for (i = 0; i < 2; i++) {
		matrix_abs[i] = ABS(matrix_abs[i]);
		// printf("%f\n", matricella[i][j] );
	}
	// printf("ok 3\n");
	printf("\n");
};

void matrix_multiply(const TYPE* A, const TYPE* B, TYPE* RES, int rows_A, int cols_A, int cols_B) {
	TYPE product;

	// k is the row of A that is being multiplied
	for (size_t k = 0; k < rows_A; ++k) {
		// is has to be multiplied to each of the columns of B
		for (size_t j = 0; j < cols_B; ++j) {
			// use a separate variable to be sure to start from 0 with the product
			product = 0;
			for (size_t i = 0; i < cols_A; ++i) {
				product += A[cols_A * k + i] * B[j + cols_B * i];
			}
			RES[k * cols_B + j] = product;
		}
	}
}

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

void swap_cols(TYPE* m, size_t rows, size_t cols, size_t col0, size_t col1) {
	// vector with the same size of one row input matrix (that is the number of columns)
	TYPE* temp = (TYPE*)calloc(cols, sizeof(TYPE));

	for (int i = 0; i < rows; i++) {
		temp[i]            = m[col0 + cols * i];
		m[col0 + cols * i] = m[col1 + cols * i];
		m[col1 + cols * i] = temp[i];
	}
	free(temp);
}

void swap_rows(TYPE* m, size_t rows, size_t cols, size_t row0, size_t row1) {
	// vector with the same size of one row input matrix (that is the number of columns)
	TYPE* temp = (TYPE*)calloc(cols, sizeof(TYPE));

	for (int i = 0; i < cols; i++) {
		temp[i]            = m[cols * row0 + i];
		m[cols * row0 + i] = m[cols * row1 + i];
		m[cols * row1 + i] = temp[i];
	}
	free(temp);
}

void stampa_matrice(TYPE* m, size_t rows, size_t cols) {
	int i;

	for (i = 0; i < rows * cols; i++) {
		printf("%f ", m[i]);
		if (i % cols == cols - 1) {
			printf("\n");
		}
	}
	// un bel a capo prima di chiudere
	printf("\n");
}

void compute_null_space(TYPE* nullspace, const TYPE a, const TYPE b, const TYPE c) {
	// chceck that determinant is zero and do something about it
	if (a * c - b * b == 0) {
	}

	if (a != 0) {
		nullspace[0] = b;
		nullspace[1] = -a;
	} else {
		// check on these too
		// assert(a == 0);
		// assert(b == 0);
		// assert(c != 0);
		nullspace[0] = c;
		nullspace[1] = -b;
	}
}

void copy(TYPE* dest, TYPE* source, size_t num_el) {
	for (size_t psi = 0; psi < num_el; psi++) {
		dest[psi] = source[psi];
	}
}

void orthonormalize_v_with_qr(TYPE* v0, TYPE* v1, const TYPE v00, const TYPE v10, const TYPE v01, const TYPE v11) {
	v0[0] = v00;
	v0[1] = v01;
	v0[2] = 1;
	v0[3] = 0;
	normalize_array(v0, 4);

	v1[0] = v10 + v01 * v01 * v10 - v00 * v01 * v11;
	v1[1] = v11 - v00 * v01 * v10 + v00 * v00 * v11;
	v1[2] = -v00 * v10 - v01 * v11;
	v1[3] = v00 * v00 + v01 * v01 + 1;
	normalize_array(v1, 4);
}

void multiply_il_v(TYPE* RES, const TYPE IL01, const TYPE IL02, const TYPE IL03, const TYPE IL12, const TYPE IL13,
                   const TYPE* v) {
	RES[0] = v[0];
	RES[1] = v[0] * IL01 + v[1];
	RES[2] = v[0] * IL02 + v[1] * IL12 + v[2];
	RES[3] = v[0] * IL03 + v[1] * IL13 + v[3];
}

// ID is 2x2
void multiply_id_v(TYPE* RES, const TYPE ID00, const TYPE ID11, const TYPE* ID, const TYPE* v) {
	RES[0] = v[0] * ID00;
	RES[1] = v[1] * ID11;
	RES[2] = v[2] * ID[0] + v[3] * ID[1];
	RES[3] = v[2] * ID[2] + v[3] * ID[3];
}

void multiply_v_il(TYPE* RES, const TYPE* v, const TYPE IL01, const TYPE IL02, const TYPE IL03, const TYPE IL12,
                   const TYPE IL13) {
	RES[0] = v[0] + v[1] * IL01 + v[2] * IL02 + v[3] * IL03;
	RES[1] = v[1] + v[2] * IL12 + v[3] * IL13;
	RES[2] = v[2];
	RES[3] = v[3];
}

void orthonormalize_v_with_qr_single(TYPE* v0, TYPE* v1) {
	// The factorization was obtained symbolically by WolframAlpha
	// by running the following query:
	//
	// QRDecomposition[{{a,e},{b,f},{c,g},{d,h}}]

	normalize_array(v0, 4);

	// To avoid numerical stability issues when multiplying too big values in the solution,
	// we scale down the second vector.
	TYPE factor = v1[0];
	if (factor < ABS(v1[1]))
		factor = ABS(v1[1]);
	if (factor < ABS(v1[2]))
		factor = ABS(v1[2]);
	if (factor < ABS(v1[3]))
		factor = ABS(v1[3]);
	// TODO: find and equivalent in C
	// factor = 1 / (factor + std::numeric_limits<TReal>::min());
	factor = 1 / (factor);

	const TYPE a = v0[0];
	const TYPE b = v0[1];
	const TYPE c = v0[2];
	const TYPE d = v0[3];
	const TYPE e = v1[0] * factor;
	const TYPE f = v1[1] * factor;
	const TYPE g = v1[2] * factor;
	const TYPE h = v1[3] * factor;

	// OPTME: We could reuse some of the multiplications.
	// ANSME: Is there a more numerically stable way to do this?
	v1[0] = b * b * e + c * c * e + d * d * e - a * b * f - a * c * g - a * d * h;
	v1[1] = -a * b * e + a * a * f + c * c * f + d * d * f - b * c * g - b * d * h;
	v1[2] = -a * c * e - b * c * f + a * a * g + b * b * g + d * d * g - c * d * h;
	v1[3] = -a * d * e - b * d * f - c * d * g + a * a * h + b * b * h + c * c * h;
	normalize_array(v1, 4);
}

void multiply_minus_v_d(TYPE* RES, TYPE* v, TYPE* D) {
	RES[0] = -v[0] * D[0];
	RES[1] = -v[1] * D[5];
	RES[2] = -v[2] * D[5] - v[3] * D[14];
	RES[3] = -v[2] * D[11] - v[3] * D[15];
}

// scalar product of two vectors
TYPE vector_scalar(const TYPE* a, const TYPE* b, const size_t num_el) {
	TYPE somma = {0.0};
	for (size_t jod = 0; jod < num_el; jod++) {
		somma += a[jod] * b[jod];
	}
	return somma;
}
