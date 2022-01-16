

/**
 * DOCUMENTATION
 *
 * Reference : 
 *
 *  An algorithm to compute the polar decomposition of a 3 × 3 matrix,
 *  Nicholas J. Higham, Vanni Noferini, 2016.
 *  Paper : https://www.researchgate.net/publication/296638898_An_algorithm_to_compute_the_polar_decomposition_of_a_3_3_matrix/link/56e29d8c08ae03f02790a388/download
 *  Source code : https://github.com/higham/polar-decomp-3by3
 *
 * Let A be a non-singular 3×3 matrice, like the ones used in channel mixer or in cameras input profiles.
 * Such matrices define transforms between RGB and XYZ spaces depending on the vector base transform.
 * The vector base is the 3 RGB primaries, defined from/to XYZ reference (absolute) space. Converting between
 * color spaces is then only a change of coordinates for the pixels color vector, depending on how the primaries
 * rotate and rescale in XYZ.
 *
 * RGB spaces conversions are therefore linear maps from old RGB to XYZ to new RGB. Geometrically, linear maps
 * can be interpreted as a combination of scalings (homothety), rotations and shear mapping (transvection).
 *
 * But they also have an interesting property : 
 *
 *   For any 3×3 invertible matrice A describing a linear map, the general linear map can be decomposed as
 *   a single 3D rotation around a particular 3D vector.
 *
 *   That is, there is a factorization of A = Q * H, where Q is
 *   the matrice of rotation around a axis of vector H.
 *
 *
 * This is interesting for us, on the GUI side. 3×3 matrices (9 params) are not intuitive to users, and the visual result of a
 * single coefficient change is hard to predict. This method allows us to reduce 9 input parameters to :
 *  * 6 : 3 angles of rotation, and the 3D coordinates of the (non-unit) rotation axis vector,
 *  * 7 : 3 angles of rotation, the 3D coordinates of the unit rotation axis vector, and a scaling factor for this vector.
 *
 * Usually, this is achieved by using HSL spaces, which suck because they work only for bounded signals in [ 0 ; 1 ].
 * Also, they are not colorspaces, not connected to either physics or psychology, so they are bad. Anyone saying
 * otherwise is fake news.
 *
 * The present method generalizes the HSL approach to XYZ, LMS and weird spaces, with none of the drawbacks of the
 * the cheapo lazy-ass maths-disabled HSL bullshit. It's great. You should try it some time. Simply the best.
 *
 **/

#include <stdio.h>
#include <complex.h>
#include "svd.h"
// is it necessary?
#include <stdbool.h>

// Note : if you review this code using the original Matlab implementation,
// remember Matlab indexes arrays from 1, while C starts at 0, so every index needs to be shifted by -1.

static inline void normalize_array(float *vector) {
	int ii;
	float norm = 0;
	// computing the square of the norm of v
	for(ii=0; ii<4; ii++){
		norm += vector[ii]*vector[ii]; 
	}
	// squared root of the norm of vector
	norm = sqrt(norm);

	//dividing each element of v for the norm
	for(ii=0; ii<4; ii++){
		vector[ii]/norm;
	}
}

//absolute value element wise for a 2x2 matrix
static inline void abs_matrix_22(float matrix_vals[][2], float matrix_abs[][2]) {
  int i, j;
    
  for(i = 0; i < 2; i++){
	  for(j = 0; j < 2; j++){
		  //printf("prima: %f\n dopo:  %f\n", matrix[i][j], fabsf(matrix[i][j]) );
		  matrix_abs[i][j] = fabsf(matrix_vals[i][j]); 
		  //printf("%f\n", matricella[i][j] );
        }
    }
  //printf("ok 3\n");
  printf("\n");
}

// multiply two matrices of any order, assuming that the multiplication is possible
void matrix_multiply( double *A, double *B, double *RES, int rows_A, int cols_A, int cols_B) {
	int ii, jj, kk;
	int rows_B;
	double product;

	if(rows_B != cols_A) {
		printf("operation not allowed for wrong dimensions\n");
	}
	
	// for simplicity of definition. if this below isn't the case the product isn't possible
	rows_B = cols_A;

	product = 0;

	for(kk=0; kk< rows_A; kk++) {
		for(jj=0; jj< cols_A; jj++) {
			for(ii=0; ii< cols_A; ii++) {
				product += A[cols_A*kk+ii]*B[cols_B*ii+jj];
			}
			RES[kk*cols_A+jj] = product;
			product = 0;
		}
	}
}


// get the value of the greater element of a 2x2 matrix
// if anyone knows how to make the function work with a 'float **matrice' I'm interested
static inline double max_val_matrix_22(float matrice[][2]) {
	int i,j;
	float max;
	max = matrice[0][0];
	for(i = 0; i < 2; i++) {
		for(j = 0; j < 2; j++){
			if (matrice[i][j] > max) {
				max = matrice[i][j];
			}
		}
	}
	return max;
}


// redefined below? che cazzo succede?
//static inline void polar_decomposition(float A[3][3], float Q[3][3], float H[3][3])
// {
//  // Frobenius / L2 norm of the matrice - aka we sum the squares of each matrice element and take the sqrt
//	float max;
//	int i,j;
//    const float norm = sqrtf(A[0][0] * A[0][0] + A[0][1] * A[0][1]  + A[0][2] * A[0][2] );
//	for(i = 0; i < 2; i++) {
//		for(j = 0; j < 2; j++){
//			if (matrice[i][j] > max) {
//				max = matrice[i][j];
//			}
//		}
//	}
//	return max;
//}

static inline void normalize_array(float *vector, int size) {
	int ii;
	float norm = 0;

	//size = sizeof(vector)/sizeof(vector[0]);

	//printf("lunghezza nella funzione: %d\n", size ); 
	//printf("lunghezza: %d\n", sizeof(vector)); 

	// computing the square of the norm of v
	for(ii=0; ii<4; ii++){
		norm += vector[ii]*vector[ii]; 
	}
	// squared root of the norm of v
	norm = sqrt(norm);

	//printf("norm is: %f\n", norm);

	//dividing each element of v for the norm
	for(ii=0; ii<size; ii++){
		vector[ii] /= norm;
	}
}


static inline void polar_decomposition(float A[3][3], float Q[3][3], float H[3][3])
{
	// Frobenius / L2 norm of the matrice - aka we sum the squares of each matrice element and take the sqrt
	const float norm = sqrtf(A[0][0] * A[0][0] + A[0][1] * A[0][1]  + A[0][2] * A[0][2] +
	                         A[1][0] * A[1][0] + A[1][1] * A[1][1]  + A[1][2] * A[1][2] +
	                         A[2][0] * A[2][0] + A[2][1] * A[2][1]  + A[2][2] * A[2][2]);
	
	
	// Normalize the matrice A in-place, so A norm is 1
	for(size_t i = 0; i < 3; i++)
		for(size_t j = 0; j < 3; j++)
			A[i][j] /= norm;
	
	// Compute the conditionning of the matrice
	float m, b;
	
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
	
	float d;
	bool subspa = false;
	
	// copy of A
	float AA[3][3] = { { A[0][0], A[0][1], A[0][2] },
	                   { A[1][0], A[1][1], A[1][2] },
	                   { A[2][0], A[2][1], A[2][2] } };
	
	
	size_t r = 0, c = 0;
	float dd = 1.f;
	
	/*
	if(b - 1.f + 1e-4f) > 0.f)
	{
	  we could use the quick path if perf is critical.
	  It's not implemented here yet because we don't use this function for each pixel and
	  for compactness of the code.
	  See the original matlab code in case you need high perf.
	}
	else
	*/
	
	// General / slow path
	
	// Search index (r, c) of the max element in matrice
	if(fabsf(A[1][0]) > fabsf(A[0][0])) r = 1; // c = 0
	if(fabsf(A[2][0]) > fabsf(A[r][c])) r = 2; // c = 0
	if(fabsf(A[0][1]) > fabsf(A[r][c])) r = 0, c = 1;
	if(fabsf(A[1][1]) > fabsf(A[r][c])) r = 1, c = 1;
	if(fabsf(A[2][1]) > fabsf(A[r][c])) r = 2, c = 1;
	if(fabsf(A[0][2]) > fabsf(A[r][c])) r = 0, c = 2;
	if(fabsf(A[1][2]) > fabsf(A[r][c])) r = 1, c = 2;
	if(fabsf(A[2][2]) > fabsf(A[r][c])) r = 2, c = 2;
	
	
	
	
	
	int k;
	
	if(r > 0)
	{
		// invert lines 0 and r
		float temp_r[3] = { AA[r][0], AA[r][1], AA[r][2] };
		float temp_0[3] = { AA[0][0], AA[0][1], AA[0][2] };
		for(size_t k = 0 ; k < 3; ++k) {
			AA[0][k] = temp_r[k];
		};
		for(size_t k = 0; k < 3; ++k) AA[r][k] = temp_0[k];
		dd = -dd;
	}
	
	if(c > 0)
	{
		// invert columns 0 and c
		float temp_c[3] = { AA[0][c], AA[1][c], AA[2][c] };
		float temp_0[3] = { AA[0][0], AA[1][0], AA[2][0] };
		for(size_t k = 0; k < 3; ++c) AA[k][0] = temp_c[k];
		for(size_t k = 0; k < 3; ++c) AA[k][c] = temp_0[k];
		dd = -dd;
	}
	
	float U[3] = { AA[0][0], 0.f, 0.f };
	
	float m0 = AA[0][1] / AA[0][0];
	float m1 = AA[0][2] / AA[0][0];
	float AAA[2][2] = { { AA[1][1] - AA[1][0] * m0, AA[1][2] - AA[1][0] * m1 },
	                    { AA[2][1] - AA[2][0] * m0, AA[2][2] - AA[2][0] * m1 } };
	
	r = 0, c = 0;
	if(fabsf(AA[1][0]) > fabsf(AA[0][0])) r = 1; // c = 0
	if(fabsf(AA[0][1]) > fabsf(AA[r][c])) r = 0, c = 1;
	if(fabsf(AA[1][1]) > fabsf(AA[r][c])) r = 1, c = 1;
	
	if(r == 1) dd = -dd;
	if(c > 0) dd = -dd;
	
	U[1] = AA[r][c];
	
	// fixed from U(2). needs check
	if(U[1] == 0) U[2] = 0;
	else U[2] = AA[2-r][2-c] - AA[r][2-c] * AA[2-r][c] / U[1];
	
	d = dd;
	dd = dd * U[0] * U[1] * U[2];
	
	if(U[0] < 0) d = -d;
	if(U[1] < 0) d = -d;
	if(U[2] < 0) d = -d;
	
	float AU = fabsf(U[2]);
	
	float nit;
	
	if(AU > 6.607e-8f)
	{
	    nit = 16.8f + 2.f * log10f(AU);
	    nit = ceilf(15.f / nit);
	}
	else
	{
		// corrected as above to get values from stdbool.h
		subspa = true;
	}
	
	if(d == 0) d = 1.f;
	
	dd = 8.f * d * dd;
	
	float t = A[0][0] + A[1][1] + A[2][2];
	
	float B[4][4] = { {   t, A[1][2] - A[2][1], A[2][0] - A[0][2], A[0][1] - A[1][0] },
	                  { 0.f, 2.f * A[0][0] - t, A[0][1] + A[1][0], A[0][2] + A[2][0] },
	                  { 0.f,               0.f, 2.f * A[1][1] - t, A[1][2] + A[2][1] },
	                  { 0.f,               0.f,               0.f, 2.f * A[2][2] - t } };
	
	for(size_t i = 0; i < 4; ++i)
		for(size_t j = 0; j < 4; ++j)
			B[i][j] /= d;
	
	
	B[1][0] = B[0][1];
	B[2][0] = B[0][2];
	B[3][0] = B[0][3];
	B[2][1] = B[1][2];
	B[3][1] = B[1][3];
	B[3][2] = B[2][3];
	
	// Find largest eigenvalue
	float x;
	
	if(b >= -0.3332f)
	{
		// Use analytic formula if matrice is well conditioned
		double complex Delta0 = 1.f + 3. * b;
		double complex Delta1 = -1. + (27. / 16.) * dd * dd + 9. * b;
		double complex phi = (Delta1 / Delta0) / csqrt(Delta0);
		double complex SS  =  (4. / 3.) * (1. + ccosf(cacosf(phi) / 3.) * csqrt(Delta0));
		double complex S = csqrt(SS) / 2.;
		
		x = (float)(creal(S) + 0.5 * sqrt(max(0., creal(-SS + 4. + dd / S))));
	}
	else
	{
		// Use Newton if matrice is ill conditioned
		// We use double precision temporarily because the solution can degenerate faster in single precision
		double x_temp = sqrt(3.);
		double xold = 3;
		while((xold - x_temp) > 1e-12)
		{
			xold = x_temp;
			double px = x_temp * (x_temp * (x_temp * x_temp - 2.) - dd) + b;
			double dpx = x_temp * (4. * x_temp * x_temp - 4.) - dd;
			x_temp = x_temp - px / dpx;
		}
		x = (float)x_temp;
	}
	
	// Again, don't do the quick path
	float BB[4][4];
	
	for(size_t i = 0; i < 4; ++i)
		for(size_t j = 0; j < 4; ++j)
		{
			BB[i][j] -= B[i][j];
			if(i == j) BB[i][j] + x; // add x on the diagonal
		}
	
	size_t p[4] = { 0, 1, 2, 3 };
	
	float L[4][4] = { { 1.f, 0.f, 0.f, 0.f },
	                  { 0.f, 1.f, 0.f, 0.f },
	                  { 0.f, 0.f, 1.f, 0.f },
	                  { 0.f, 0.f, 0.f, 1.f } };
	
	float D[4][4] = { { 0.f } };
	
	// First step
	r = 3;
	if(BB[3][3] < BB[2][2]) r = 2;
	if(BB[r][r] < BB[1][1]) r = 1;
	if(BB[r][r] > BB[0][0])
	{
		// p([1 r(1)]) = [r(1) 1]
		p[0] = r;
		p[r] = 0;
		
		// BB = BB(p,p);
		for(size_t i = 0; i < 4; ++i)
			for(size_t j = 0; j < 4; ++j)
				BB[i][j] = BB[p[i]][p[j]];
	}
	
	float D = BB[0][0];
	
	L[1][0] = BB[1][0];
	L[2][0] = BB[2][0];
	L[3][0] = BB[3][0];
	
	BB[1][1] = BB[1][1] - L[1][0] * BB[0][1];
	BB[2][1] = BB[2][1] - L[1][0] * BB[0][2];
	BB[1][2] = BB[2][1];
	
	BB[3][1] = BB[3][1] - L[1][0] * BB[0][3];
	BB[1][3] = BB[3][1];
	D[3][2] = BB[3][2]-L[2][1]*BB[1][3];
	D[2][3] = D[3][2];
	D[3][3] = BB[3][3]-L[3][1]*BB[1][3];
	
	
	
	
	
	
	 //=============================================================================
	
	// is this the right size for v?
	double v[4];	
	
	DD = D[2][2]*D[3][3]-D[2][3]*D[2][3];
	
	// slamming the cock for now about the case where DD==0
	//if (DD == 0) {
	//    // variables needed for the computations below
	//    double D_m[2][2];
	//    double abs_matrix_22[2][2];
	//    double max;
	//	// south-est minor of D
	//	D_m[0][0] = D[1][1];
	//	D_m[0][1] = D[1][2];
	//	D_m[1][1] = D[2][2];
	//	D_m[1][0] = D[2][1];
	//	// absolute value of each element written on matrix_abs
	//	abs_matrix_22( D_m, matrix_abs);	
	//	// max among all the elements
	//	max = max_val_matrix_22(abs_matrix);
	
	//	if(max == 0) {
	//		v[0] = L[1][0]*L[3][1]-L[3][0];
	//		v[1] = -L[3][1];
	//		v[2] = 0;
	//		v[3] = 1;
	
	//		int ii;
	//		double norm = 0;
	
	//		// computing the square of the norm of v
	//		for(ii=0; ii<4; ii++){
	//			norm += v[ii]*v[ii]; 
	//		}
	//		// squared root of the norm of v
	//		norm = sqrt(norm);
	
	//		//dividing each element of v for the norm
	//		for(ii=0; ii<4; ii++){
	//			v[ii]/norm;
	//		}
	//	} else {
	//		// NOT COMPLETED AS 'null' IS A PAIN IN THE ASS
	//		v = L'\[0;0;null(D(2:3,2:3))];
	//		v = v/norm(v);
	//	  }
	//} else {
	//	//ID = [D(4,4) -D(3,4); -D(3,4) D(3,3)];
	//	// can the matrix be directly assignes like this?
	//	double ID[2][2] = {{D[3][3], D[2][3]}, {D[2][3], D[2][2]}};
	
	double ID[2][2] = {{D[3][3], D[2][3]}, {D[2][3], D[2][2]}};
	
		


// =============================================================================
//      if max(abs(D(2:3,2:3))) == 0, v = [L(1,0)*L(3,1)-L(3,0);-L(3,1);0;1];
//  		v = v/norm(v);
//      else
//          v = L'\[0;0;null(D(2:3,2:3))];v = v/norm(v);
//      end        
//  else
//  ID = [D(3,3) -D(2,3); -D(2,3) D(2,2)];
// =============================================================================



// =============================================================================
  
	if (subspa == 1) {
	    //v = [L(2,1)*L(3,2)-L(3,1) L(2,1)*L(4,2)-L(4,1);-L(3,2) -L(4,2);1 0;0 1];
	  
		v[0][0] = L[1][0]*L[2][1]-L[2][0];
		v[0][1] = L[1][0]*L[3][1]-L[3][0];
		v[1][0] = -L[2][1];
		v[1][1] = -L[3][1];
		v[2][0] = 1;
		v[2][1] = 0;
		v[3][0] = 0
		v[3][1] = 1;
		
		// first use of L here, it needs to be defined, after the dimensions are clear
		  //IL = [1 0 0 0;-L(2,1) 1 0 0;v'];
		
		IL[0][0] = 1;
		IL[0][1] = 0;
		IL[0][2] = 0;
		IL[0][3] = 0;
		
		IL[1][0] = -L[1,0];
		IL[1][1] = 0;
		IL[1][2] = 0;
		IL[1][3] = 0;
		
		// assigning v as transposed already
		IL[2][0] = v[0][0];
		IL[2][1] = v[0][1];
		IL[2][2] = v[0][2];
		IL[2][3] = v[0][3];
		
		IL[3][0] = v[1][0];
		IL[3][1] = v[1][1];
		IL[3][2] = v[1][2];
		IL[3][3] = v[1][3];
		
		
		//==========================================================================
		
		// The tilde is used to get part of the output of qr and discrd the remaining
		  //[v ~] = qr(v,0); //->cost in flops if implemented by hand: 37 M+24 A+4 O
		// Since v and Q (resulting from qr decomp) have the same order, the v = Q.
		// The behaviour of the implemented function for QR works for the present case
		// since if A(mxn) has m>n the output is the same.
		int row = 4;
		int col = 2;
		
		double Q[row][col];
		
		QR_decomp(v,Q,R,row,col);
		
		// assigning Q values to v. For the present case size(Q) = size(v)
		
		int ii;
		for (ii=0; ii<row ;ii++) {
			for (jj=0; jj<col ;jj++) {
				v[ii][jj] = Q[ii][jj];
			}
		}
		
		// IL is 4x4
		
		double v_m[4][2];
		
		matrix_multiplication(IL, v, v_m, 4,4,2);
		
		  //v = IL*v;  // it looks faster to multiply than to solve lin syst (even though should be same flops)
		
		  //v(1,:) = v(1,:)/D(1,1);
		for(ii=0; ii<2; ii++) {
			v_m[0][ii] /= D[0][0];
		}
		
		  //v(2,:) = v(2,:)/D(2,2);
		for(ii=0; ii<2; ii++) {
			v_m[1][ii] /= D[1][1];
		}
		
		// a part of v (which is v_m in the current case) needs to be extracted and later concatenated
		
		  //v(3:4,:) = ID*v(3:4,:)/DD(1);
		
		//matrix_multiplication(ID, v_m, v_mm, 0,0,0);
		// result matrix for the next multiplication, which is a submatrix of v[4,2], so [2,2]
		double v_res[2][2];
		// submatrix extracted from v
		double v_mm[2][2];
		
		v_mm[0][0] = v[2][0];
		v_mm[0][1] = v[2][1];
		v_mm[1][0] = v[3][0];
		v_mm[1][1] = v[3][1];
		
		// ID has been previously defined as [2][2]
		matrix_multiplication(ID, v_m, v_res, 2,2,2);
		
		// changing the values of the original matrix and diving v_res at the same time
		v[2][0] = v_res[0][0]/DD[0];
		v[2][1] = v_res[0][1]/DD[0];
		v[3][0] = v_res[1][0]/DD[0];
		v[3][1] = v_res[1][1]/DD[0];
		
		
		// v_mm and v_res could now be freed, if they where dinamically allocated, but since I am stupido ai chènnòt du só, plis fri de mètrisis fór mi.
		
		
		// transposed matrix
		double v_t[2,4];
		
		v_t[0][0] = v[0][0];
		v_t[0][1] = v[1][0];
		v_t[0][2] = v[2][0];
		v_t[0][3] = v[3][0];
		
		v_t[1][0] = v[0][1];
		v_t[1][1] = v[1][1];
		v_t[1][2] = v[2][1];
		v_t[1][3] = v[3][1];
		
		
		// v = v'*IL;
		// put the result in v
		// need for a temporary matrix since v changes size
		// this also may be freed
		double v_res1[2][4];
		matrix_multiply( v_t, IL, v_res1, 2,4,4);
		// after this v is 2x4
		
		// traspose may be freed now and another one is required
		
		// aritransposed matrix
		// v = v';
		double v_tt[2,4];
		
		v_tt[0][0] = v_res1[0][0];
		v_tt[1][0] = v_res1[0][1];
		v_tt[2][0] = v_res1[0][2];
		v_tt[3][0] = v_res1[0][3];
		
		v_tt[0][1] = v_res1[1][0];
		v_tt[1][1] = v_res1[1][1];
		v_tt[2][1] = v_res1[1][2];
		v_tt[3][1] = v_res1[1][3];
		
		// next multiplication, the v on the right hand is transposed
		// from the previous operation
		  // v = IL*v;
		// equivalent to this
		// v = IL*v'
		double v_res2[4][2];
		matrix_multiply( IL, v_tt, v_res2, 4,4,2);
		//v the result is now 4x2
		
		// dumb to not use anpther for cycle but I'm already loosing my sight
		  // v(1,:) = v(1,:)/D(1,1);
		  // v(2,:) = v(2,:)/D(2,2);
		for(j=0; j<2; j++) {
			v_res2[0,j] /= D[0,0];
			v_res2[1,j] /= D[1,1];
		};
		
		// ID is 2x2
		// double sub_v[2][2];
		
		// v(3:4,:) = ID*v(3:4,:)/DD(1);
		for(j=0; j<2; j++) {
			v[2,j] = ID[2,j]*D[0,0];
			v[3,j] = ID[3,j]*D[1,1];
		};
	
	
	
	
	  /*
	  
	
	  v = v'*IL;v = v';
	          [v ~] = qr(v,0);
	  H = v'*L;H = -H*D*H'; %Cheaper
	
	  if (fabsf(H(1,2))<1e-15) {
	      if ( [H(1][1] > H[1][2] ) {
	          v = v(:,1);
	  	else {
	          v = v(:,2);
	  	  ];
	  } else {
	        r = (H[1][1]-H[2][2]) / (2*H([1][2]);
	        v = v*[r+sign(H(1,2))*sqrt(1+r(1)*r(1));1];
	        v = v/norm(v);
	  };
	
	//else
	//    v = [L(2,1)*L(4,2)+L(3,1)*L(4,3)-L(2,1)*L(4,3)*L(3,2)-L(4,1);
	//         L(4,3)*L(3,2)-L(4,2); -L(4,3) ;1];
	//    IL = [1 0 0 0; -L(2,1) 1 0 0; L(2,1)*L(3,2)-L(3,1) -L(3,2) 1 0; v'];
	//    nv = realsqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3)+v(4)*v(4));
	//    v = v/nv(1);
	
	//    for (it = 1:nit) {
	//        v = IL*v;
	//        v(1) = v(1)/D(1,1);v(2) = v(2)/D(2,2);v(3:4) = ID*v(3:4)/DD(1);
	//        v = v'*IL;v = v';
	//        nv = realsqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3)+v(4)*v(4));
	//        v = v/nv(1);
	//	};
	//};
	end
	v(p) = v;
	
	*/
 }
