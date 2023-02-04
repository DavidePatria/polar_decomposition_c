#include "polar_decomposition.h"
#include <stdio.h>

void print_float(TYPE num);
void check_orthonormal(TYPE Q[3][3]);
void stampa_3x3(TYPE Q[3][3]);
void copy(TYPE* dest, TYPE* source, size_t num_el);

#ifndef TEST_N
#define TEST_N 1
#endif

int main() {

	TYPE A[3][3];
	// simple case one
	switch (TEST_N) {
	case 1: {
		TYPE temp[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
		copy(*A, *temp, 9);
	}
	case 2: {
		TYPE temp[3][3] = {{0.1, 0.1, 0.3}, {0.2, 0.1, 0.2}, {0.3, 0.0, 0.1}};
		copy(*A, *temp, 9);
	}
	// TODO: check and use here the maledict test mentioned in the cpp version
	case 3: {
		TYPE temp[3][3] = {{-1, 2, 0}, {1.3, 0.4, 3.4}, {0.3, 2.4, 3.4}};
		copy(*A, *temp, 9);
	}
	}

	TYPE Q[3][3];
	TYPE H[3][3];

	// stampa_3x3(A);
	polar_decomposition(A, Q, H);

	// stampa_3x3(A);
	stampa_3x3(Q);
	stampa_3x3(H);

	matrix_multiply(*(Q), *(H), *(A), 3, 3, 3);

	// stampa_3x3(A);

	check_orthonormal(Q);

	return 0;
}

void stampa_3x3(TYPE Q[3][3]) {
	for (size_t riga = 0; riga < 3; riga++) {
		for (size_t col = 0; col < 3; col++) {
			printf("%f ", Q[riga][col]);
			if (col == 2)
				printf("\n");
		}
	}
	printf("\n");
}

void check_orthonormal(TYPE Q[3][3]) {
	TYPE scalar01 = 0;
	TYPE scalar12 = 0;
	TYPE scalar02 = 0;

	for (size_t i = 0; i < 3; i++) {
		scalar01 += Q[0][i] * Q[1][i];
		scalar12 += Q[1][i] * Q[2][i];
		scalar02 += Q[0][i] * Q[2][i];
	}
	print_float(scalar01);
	print_float(scalar12);
	print_float(scalar02);
}

void print_float(TYPE num) { printf("scalar12: %f\n", num); }
