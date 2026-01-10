#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <A.fits> <B.fits>\n", argv[0]);
        return 1;
    }

    double *A, *B;
    long Na, Pa, Nb, Pb;
    int xa, ya;

    read_fits(argv[1], &A, &Na, &Pa, &xa, &ya, NULL);
    read_fits(argv[2], &B, &Nb, &Pb, &xa, &ya, NULL);

    if (Na != Nb || Pa != Pb) {
        fprintf(stderr, "Dimensions mismatch: %ldx%ld vs %ldx%ld\n", Na, Pa, Nb, Pb);
        return 1;
    }

    printf("Correlation Matrix Diagonals (Ai vs Bi):\n");
    for (long k = 0; k < Pa; k++) {
        double sumA = 0, sumB = 0, sumAB = 0, sumA2 = 0, sumB2 = 0;
        for (long i = 0; i < Na; i++) {
            double valA = A[i * Pa + k];
            double valB = B[i * Pb + k];
            sumA += valA;
            sumB += valB;
        }
        double meanA = sumA / Na;
        double meanB = sumB / Na;

        for (long i = 0; i < Na; i++) {
            double valA = A[i * Pa + k] - meanA;
            double valB = B[i * Pb + k] - meanB;
            sumAB += valA * valB;
            sumA2 += valA * valA;
            sumB2 += valB * valB;
        }
        double corr = sumAB / sqrt(sumA2 * sumB2);
        printf("  Mode %ld: r = %.6f, stdA = %.6f, stdB = %.6f\n", k, corr, sqrt(sumA2/(Na-1)), sqrt(sumB2/(Na-1)));
    }

    free(A); free(B);
    return 0;
}
