#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Calculate Pearson correlation or covariance between two vectors of length N
double calc_stat(double *vecA, double *vecB, long N, int do_cov) {
    double sumA = 0, sumB = 0, sumAB = 0, sumA2 = 0, sumB2 = 0;

    for (long i = 0; i < N; i++) {
        sumA += vecA[i];
        sumB += vecB[i];
    }
    double meanA = sumA / N;
    double meanB = sumB / N;

    for (long i = 0; i < N; i++) {
        double valA = vecA[i] - meanA;
        double valB = vecB[i] - meanB;
        sumAB += valA * valB;
        sumA2 += valA * valA;
        sumB2 += valB * valB;
    }

    if (do_cov) {
        return sumAB / N;
    } else {
        if (sumA2 == 0 || sumB2 == 0) return 0.0;
        return sumAB / sqrt(sumA2 * sumB2);
    }
}

int main(int argc, char *argv[]) {
    int do_cov = 0;
    int arg_offset = 0;

    if (argc > 1 && strcmp(argv[1], "-cov") == 0) {
        do_cov = 1;
        arg_offset = 1;
    }

    if (argc - arg_offset != 3) {
        fprintf(stderr, "Usage: %s [-cov] <A.fits> <B.fits>\n", argv[0]);
        return 1;
    }

    double *A, *B;
    long Na, Pa, Nb, Pb;
    int xa, ya;

    read_fits(argv[1 + arg_offset], &A, &Na, &Pa, &xa, &ya, NULL);
    read_fits(argv[2 + arg_offset], &B, &Nb, &Pb, &xa, &ya, NULL);

    if (Na != Nb) {
        fprintf(stderr, "Sample count mismatch: %ld vs %ld\n", Na, Nb);
        return 1;
    }

    // Allocate temp vectors for extraction
    double *vecA = malloc(Na * sizeof(double));
    double *vecB = malloc(Na * sizeof(double));

    printf("| A \\ B |");
    for (long j = 0; j < Pb; j++) {
        printf(" B%ld |", j);
    }
    printf("\n|");
    for (long j = 0; j <= Pb; j++) {
        printf("---|");
    }
    printf("\n");

    for (long i = 0; i < Pa; i++) {
        printf("| **A%ld** |", i);
        // Extract Col A_i
        for (long k = 0; k < Na; k++) vecA[k] = A[k * Pa + i];

        for (long j = 0; j < Pb; j++) {
            // Extract Col B_j
            for (long k = 0; k < Na; k++) vecB[k] = B[k * Pb + j];

            double val = calc_stat(vecA, vecB, Na, do_cov);

            // Highlight off-diagonals if non-zero (e.g. > 0.01)
            // But main diagonal is i==j (if Pa==Pb)

            if (fabs(val) < 1e-4) {
                printf(" %.4f |", 0.0);
            } else {
                printf(" **%.4f** |", val);
            }
        }
        printf("\n");
    }

    free(vecA); free(vecB);
    free(A); free(B);
    return 0;
}
