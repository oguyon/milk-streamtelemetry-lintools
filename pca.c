#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cblas.h>
#include <lapacke.h>
#include "common.h"

int main(int argc, char *argv[]) {
    int ncpu = 0;
    int use_float = 0;
    int arg_offset = 0;

    // Argument parsing
    int i = 1;
    while (i < argc) {
        if (strcmp(argv[i], "-ncpu") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: -ncpu requires an argument.\n");
                return 1;
            }
            ncpu = atoi(argv[i+1]);
            if (ncpu < 1) {
                fprintf(stderr, "Error: ncpu must be >= 1. Got %d.\n", ncpu);
                return 1;
            }
            i += 2;
            arg_offset += 2;
        } else if (strcmp(argv[i], "-float") == 0) {
            use_float = 1;
            i += 1;
            arg_offset += 1;
        } else {
            break; // Positional argument found
        }
    }

    if (argc - arg_offset != 5) {
        fprintf(stderr, "Usage: %s [-ncpu <n>] [-float] <npca> <input.fits> <modes.fits> <coeffs.fits>\n", argv[0]);
        return 1;
    }

    if (ncpu > 0) {
        openblas_set_num_threads(ncpu);
        printf("Using %d CPU cores.\n", ncpu);
    }
    if (use_float) {
        printf("Using Single Precision (float).\n");
    } else {
        printf("Using Double Precision (double).\n");
    }

    int npca = atoi(argv[1 + arg_offset]);
    const char *infile = argv[2 + arg_offset];
    const char *modes_file = argv[3 + arg_offset];
    const char *coeffs_file = argv[4 + arg_offset];

    if (npca < 1) {
        fprintf(stderr, "Error: npca must be >= 1. Got %d.\n", npca);
        return 1;
    }

    printf("PCA Configuration:\n");
    printf("  npca: %d\n", npca);
    printf("  Input: %s\n", infile);
    printf("  Modes Output: %s\n", modes_file);
    printf("  Coeffs Output: %s\n", coeffs_file);

    if (use_float) {
        float *X = NULL;
        long N, P;
        int xa, ya;
        int naxis;

        printf("Reading %s...\n", infile);
        read_fits_float(infile, &X, &N, &P, &xa, &ya, &naxis);
        printf("  Dimensions: %ld samples x %ld pixels (%d x %d)\n", N, P, xa, ya);

        if (npca > N) {
            fprintf(stderr, "Warning: npca (%d) > N (%ld). Reducing npca to N.\n", npca, N);
            npca = (int)N;
        }
        if (npca > P) {
            fprintf(stderr, "Warning: npca (%d) > P (%ld). Reducing npca to P.\n", npca, P);
            npca = (int)P;
        }

        // SVD
        long K = (N < P) ? N : P;
        float *S = (float *)malloc(K * sizeof(float));
        float *U = (float *)malloc(N * K * sizeof(float));
        float *Vt = (float *)malloc(K * P * sizeof(float));

        printf("Computing SVD (using sgesdd)...\n");
        lapack_int info = LAPACKE_sgesdd(LAPACK_ROW_MAJOR, 'S', N, P, X, P,
                                         S, U, K, Vt, P);
        if (info != 0) {
            fprintf(stderr, "LAPACKE_sgesdd failed with error code %d\n", info);
            return 1;
        }

        // Coeffs = U * S (first npca columns)
        float *Coeffs = (float *)malloc(N * npca * sizeof(float));
        for(long i=0; i<N; i++) {
            for(int j=0; j<npca; j++) {
                Coeffs[i * npca + j] = U[i * K + j] * S[j];
            }
        }

        // Check sign of the first mode
        double coeff_mean = 0.0;
        for(long i=0; i<N; i++) {
            coeff_mean += Coeffs[i * npca + 0];
        }
        coeff_mean /= N;

        if (coeff_mean < 0) {
            printf("Flipping sign of Mode 0 to ensure positive average coefficient.\n");
            // Flip Coeffs column 0
            for(long i=0; i<N; i++) {
                Coeffs[i * npca + 0] *= -1.0f;
            }
            // Flip Vt row 0 (which becomes Mode 0)
            for(long j=0; j<P; j++) {
                Vt[0 * P + j] *= -1.0f;
            }
        }

        // Modes = Vt (first npca rows)
        float *Modes = (float *)malloc(npca * P * sizeof(float));
        for(int i=0; i<npca; i++) {
            for(long j=0; j<P; j++) {
                Modes[i * P + j] = Vt[i * P + j];
            }
        }

        printf("Writing Modes to %s...\n", modes_file);
        write_fits_3d_float(modes_file, Modes, xa, ya, npca);

        printf("Writing Coeffs to %s...\n", coeffs_file);
        write_fits_2d_float(coeffs_file, Coeffs, npca, N);

        printf("Done.\n");

        free(X); free(S); free(U); free(Vt);
        free(Coeffs); free(Modes);

    } else {
        // Double precision (original code)
        double *X = NULL;
        long N, P;
        int xa, ya;
        int naxis;

        printf("Reading %s...\n", infile);
        read_fits(infile, &X, &N, &P, &xa, &ya, &naxis);
        printf("  Dimensions: %ld samples x %ld pixels (%d x %d)\n", N, P, xa, ya);

        if (npca > N) {
            fprintf(stderr, "Warning: npca (%d) > N (%ld). Reducing npca to N.\n", npca, N);
            npca = (int)N;
        }
        if (npca > P) {
            fprintf(stderr, "Warning: npca (%d) > P (%ld). Reducing npca to P.\n", npca, P);
            npca = (int)P;
        }

        // SVD
        long K = (N < P) ? N : P;
        double *S = (double *)malloc(K * sizeof(double));
        double *U = (double *)malloc(N * K * sizeof(double));
        double *Vt = (double *)malloc(K * P * sizeof(double));

        printf("Computing SVD (using dgesdd)...\n");
        // Use dgesdd for speed (Divide and Conquer)
        lapack_int info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'S', N, P, X, P,
                                         S, U, K, Vt, P);
        if (info != 0) {
            fprintf(stderr, "LAPACKE_dgesdd failed with error code %d\n", info);
            return 1;
        }

        // Coeffs = U * S (first npca columns)
        double *Coeffs = (double *)malloc(N * npca * sizeof(double));
        for(long i=0; i<N; i++) {
            for(int j=0; j<npca; j++) {
                Coeffs[i * npca + j] = U[i * K + j] * S[j];
            }
        }

        // Check sign of the first mode
        double coeff_mean = 0.0;
        for(long i=0; i<N; i++) {
            coeff_mean += Coeffs[i * npca + 0];
        }
        coeff_mean /= N;

        if (coeff_mean < 0) {
            printf("Flipping sign of Mode 0 to ensure positive average coefficient.\n");
            // Flip Coeffs column 0
            for(long i=0; i<N; i++) {
                Coeffs[i * npca + 0] *= -1.0;
            }
            // Flip Vt row 0 (which becomes Mode 0)
            for(long j=0; j<P; j++) {
                Vt[0 * P + j] *= -1.0;
            }
        }

        // Modes = Vt (first npca rows)
        double *Modes = (double *)malloc(npca * P * sizeof(double));
        for(int i=0; i<npca; i++) {
            for(long j=0; j<P; j++) {
                Modes[i * P + j] = Vt[i * P + j];
            }
        }

        printf("Writing Modes to %s...\n", modes_file);
        write_fits_3d(modes_file, Modes, xa, ya, npca);

        printf("Writing Coeffs to %s...\n", coeffs_file);
        write_fits_2d(coeffs_file, Coeffs, npca, N);

        printf("Done.\n");

        free(X); free(S); free(U); free(Vt);
        free(Coeffs); free(Modes);
    }

    return 0;
}
