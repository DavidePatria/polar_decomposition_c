#include "polar_decomposition.h"
#include <stdio.h>

void stampa_3x3(float Q[3][3]) {
	for (size_t riga = 0; riga < 3; riga++) {
		for (size_t col = 0; col < 3; col++) {
			printf("%f ", Q[riga][col]);
			if (col == 2)
				printf("\n");
		}
	}
	printf("\n");
}

int main() {

	float A[3][3] = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
	float Q[3][3];
	float H[3][3];

	// stampa_3x3(A);
	polar_decomposition(A, Q, H);

	stampa_3x3(A);
	stampa_3x3(Q);
	stampa_3x3(H);

	float Q_flat[9];
	flatten_matrix_3x3(Q, Q_flat);

	float H_flat[9];
	flatten_matrix_3x3(H, H_flat);

	float A_flat[9];
	flatten_matrix_3x3(A, A_flat);

	matrix_multiply(Q_flat, H_flat, A_flat, 3,3, 3);

	stampa_3x3(A);


	return 0;
}
