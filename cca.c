#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cblas.h>
#include <lapacke.h>
#include "common.h"

// Forward declarations
void perform_cca(double *X, long N, long P, double *Y, long Q, int nvec, double **A_vec, double **B_vec);
void perform_pca(double *X, long N, long P, int npca, double **Coeffs, double **Modes, const char* out_filename, int xa, int ya);

int main(int argc, char *argv[]) {
    int npca = 0;
    int arg_offset = 0;

    if (argc >= 6 && strcmp(argv[1], "-npca") == 0) {
        npca = atoi(argv[2]);
        if (npca < 1) {
            fprintf(stderr, "Error: npca must be >= 1. Got %d.\n", npca);
            return 1;
        }
        arg_offset = 2;
    } else if (argc != 4) {
        fprintf(stderr, "Usage: %s [-npca <n>] <nvec> <A.fits> <B.fits>\n", argv[0]);
        return 1;
    }

    if (argc - arg_offset != 4) {
         fprintf(stderr, "Usage: %s [-npca <n>] <nvec> <A.fits> <B.fits>\n", argv[0]);
         return 1;
    }

    int nvec = atoi(argv[1 + arg_offset]);
    const char *fileA = argv[2 + arg_offset];
    const char *fileB = argv[3 + arg_offset];

    printf("CCA Configuration:\n");
    printf("  nvec: %d\n", nvec);
    if (npca > 0) printf("  npca: %d\n", npca);
    printf("  File A: %s\n", fileA);
    printf("  File B: %s\n", fileB);

    double *X = NULL, *Y = NULL;
    long Na, Pa, Nb, Qb;
    int xa, ya, xb, yb;

    // Read Data
    printf("Reading %s...\n", fileA);
    read_fits(fileA, &X, &Na, &Pa, &xa, &ya);
    printf("  Dimensions: %ld samples x %ld pixels (%d x %d)\n", Na, Pa, xa, ya);

    printf("Reading %s...\n", fileB);
    read_fits(fileB, &Y, &Nb, &Qb, &xb, &yb);
    printf("  Dimensions: %ld samples x %ld pixels (%d x %d)\n", Nb, Qb, xb, yb);

    if (Na != Nb) {
        fprintf(stderr, "Error: Number of samples (N) must match between A and B.\n");
        free(X); free(Y);
        return 1;
    }
    long N = Na;

    if (nvec > N) {
        fprintf(stderr, "Warning: nvec (%d) > N (%ld). Reducing nvec to N.\n", nvec, N);
        nvec = (int)N;
    }

    // Center Data
    printf("Centering data...\n");
    center_columns(X, N, Pa);
    center_columns(Y, N, Qb);

    double *A_vec_cca = NULL, *B_vec_cca = NULL;

    if (npca > 0) {
        // PCA Mode
        if (npca > N) {
            fprintf(stderr, "Warning: npca (%d) > N (%ld). Reducing npca to N.\n", npca, N);
            npca = (int)N;
        }
        if (npca > Pa) npca = (int)Pa;
        if (npca > Qb) npca = (int)Qb;

        printf("Performing PCA on A (keeping %d modes)...\n", npca);
        double *CoeffsA = NULL, *ModesA = NULL;
        perform_pca(X, N, Pa, npca, &CoeffsA, &ModesA, "pcaA.fits", xa, ya);

        printf("Performing PCA on B (keeping %d modes)...\n", npca);
        double *CoeffsB = NULL, *ModesB = NULL;
        perform_pca(Y, N, Qb, npca, &CoeffsB, &ModesB, "pcaB.fits", xb, yb);

        double *Wa = NULL, *Wb = NULL;
        if (nvec > npca) {
             fprintf(stderr, "Warning: nvec (%d) > npca (%d). Reducing nvec to npca.\n", nvec, npca);
             nvec = npca;
        }

        printf("Performing CCA on PCA coefficients...\n");
        perform_cca(CoeffsA, N, npca, CoeffsB, npca, nvec, &Wa, &Wb);

        // Reconstruct
        A_vec_cca = (double *)malloc(nvec * Pa * sizeof(double));
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    nvec, Pa, npca,
                    1.0, Wa, npca,
                    ModesA, Pa,
                    0.0, A_vec_cca, Pa);

        B_vec_cca = (double *)malloc(nvec * Qb * sizeof(double));
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    nvec, Qb, npca,
                    1.0, Wb, npca,
                    ModesB, Qb,
                    0.0, B_vec_cca, Qb);

        free(CoeffsA); free(ModesA);
        free(CoeffsB); free(ModesB);
        free(Wa); free(Wb);

    } else {
        // Raw CCA
        if (nvec > Pa) {
            fprintf(stderr, "Warning: nvec (%d) > Pa (%ld). Reducing nvec to Pa.\n", nvec, Pa);
            nvec = (int)Pa;
        }
        if (nvec > Qb) {
            fprintf(stderr, "Warning: nvec (%d) > Qb (%ld). Reducing nvec to Qb.\n", nvec, Qb);
            nvec = (int)Qb;
        }

        perform_cca(X, N, Pa, Y, Qb, nvec, &A_vec_cca, &B_vec_cca);
    }

    printf("Writing ccaA.fits...\n");
    write_fits_3d("ccaA.fits", A_vec_cca, xa, ya, nvec);

    printf("Writing ccaB.fits...\n");
    write_fits_3d("ccaB.fits", B_vec_cca, xb, yb, nvec);

    printf("Done.\n");

    free(X); free(Y);
    free(A_vec_cca); free(B_vec_cca);
    return 0;
}

void perform_pca(double *X, long N, long P, int npca, double **Coeffs, double **Modes, const char* out_filename, int xa, int ya) {
    long K = (N < P) ? N : P;
    double *S = (double *)malloc(K * sizeof(double));
    double *U = (double *)malloc(N * K * sizeof(double));
    double *Vt = (double *)malloc(K * P * sizeof(double));
    double *superb = (double *)malloc((K - 1) * sizeof(double));

    LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'S', 'S', N, P, X, P,
                   S, U, K, Vt, P, superb);

    *Coeffs = (double *)malloc(N * npca * sizeof(double));
    for(long i=0; i<N; i++) {
        for(int j=0; j<npca; j++) {
            (*Coeffs)[i * npca + j] = U[i * K + j] * S[j];
        }
    }

    *Modes = (double *)malloc(npca * P * sizeof(double));
    for(int i=0; i<npca; i++) {
        for(long j=0; j<P; j++) {
            (*Modes)[i * P + j] = Vt[i * P + j];
        }
    }

    write_fits_3d(out_filename, *Modes, xa, ya, npca);

    free(S); free(U); free(Vt); free(superb);
}

void perform_cca(double *X, long N, long P, double *Y, long Q, int nvec, double **A_vec, double **B_vec) {
    long Kx = (N < P) ? N : P;
    long Ky = (N < Q) ? N : Q;

    double *Rx = (double *)malloc(Kx * P * sizeof(double));
    double *Ry = (double *)malloc(Ky * Q * sizeof(double));
    double *tauX = (double *)malloc(Kx * sizeof(double));
    double *tauY = (double *)malloc(Ky * sizeof(double));

    LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, N, P, X, P, tauX);
    memset(Rx, 0, Kx * P * sizeof(double));
    for (long i = 0; i < Kx; i++) {
        for (long j = i; j < P; j++) {
            Rx[i * P + j] = X[i * P + j];
        }
    }
    LAPACKE_dorgqr(LAPACK_ROW_MAJOR, N, Kx, Kx, X, P, tauX);

    LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, N, Q, Y, Q, tauY);
    memset(Ry, 0, Ky * Q * sizeof(double));
    for (long i = 0; i < Ky; i++) {
        for (long j = i; j < Q; j++) {
            Ry[i * Q + j] = Y[i * Q + j];
        }
    }
    LAPACKE_dorgqr(LAPACK_ROW_MAJOR, N, Ky, Ky, Y, Q, tauY);

    double *M = (double *)malloc(Kx * Ky * sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                Kx, Ky, N,
                1.0, X, P,
                Y, Q,
                0.0, M, Ky);

    double *S = (double *)malloc(((Kx < Ky) ? Kx : Ky) * sizeof(double));
    double *U = (double *)malloc(Kx * Kx * sizeof(double));
    double *Vt = (double *)malloc(Ky * Ky * sizeof(double));
    double *superb = (double *)malloc((((Kx < Ky) ? Kx : Ky) - 1) * sizeof(double));

    LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', Kx, Ky, M, Ky,
                   S, U, Kx, Vt, Ky, superb);

    printf("Canonical Correlations:\n");
    for(int i=0; i<nvec; i++) {
        if (i < ((Kx < Ky) ? Kx : Ky)) {
            printf("  Mode %d: %g\n", i, S[i]);
        }
    }

    *A_vec = (double *)malloc(nvec * P * sizeof(double));
    *B_vec = (double *)malloc(nvec * Q * sizeof(double));

    double *RHS_A = (double *)calloc(P * nvec, sizeof(double));
    for (int j = 0; j < nvec; j++) {
        for (int i = 0; i < Kx; i++) {
            RHS_A[i * nvec + j] = U[i * Kx + j];
        }
    }

    double *Rx_tmp = (double *)malloc(Kx * P * sizeof(double));
    memcpy(Rx_tmp, Rx, Kx * P * sizeof(double));
    double *S_Rx = (double *)malloc(((Kx < P) ? Kx : P) * sizeof(double));
    lapack_int rank_Rx;
    LAPACKE_dgelss(LAPACK_ROW_MAJOR, Kx, P, nvec, Rx_tmp, P, RHS_A, nvec, S_Rx, -1, &rank_Rx);

    for (int j = 0; j < nvec; j++) {
        for (long i = 0; i < P; i++) {
            (*A_vec)[j * P + i] = RHS_A[i * nvec + j];
        }
    }
    free(S_Rx); free(Rx_tmp); free(RHS_A);

    double *RHS_B = (double *)calloc(Q * nvec, sizeof(double));
    for (int j = 0; j < nvec; j++) {
        for (int i = 0; i < Ky; i++) {
            RHS_B[i * nvec + j] = Vt[j * Ky + i];
        }
    }

    double *Ry_tmp = (double *)malloc(Ky * Q * sizeof(double));
    memcpy(Ry_tmp, Ry, Ky * Q * sizeof(double));
    double *S_Ry = (double *)malloc(((Ky < Q) ? Ky : Q) * sizeof(double));
    lapack_int rank_Ry;
    LAPACKE_dgelss(LAPACK_ROW_MAJOR, Ky, Q, nvec, Ry_tmp, Q, RHS_B, nvec, S_Ry, -1, &rank_Ry);

    for (int j = 0; j < nvec; j++) {
        for (long i = 0; i < Q; i++) {
            (*B_vec)[j * Q + i] = RHS_B[i * nvec + j];
        }
    }
    free(S_Ry); free(Ry_tmp); free(RHS_B);

    free(Rx); free(Ry);
    free(tauX); free(tauY);
    free(M); free(S); free(U); free(Vt); free(superb);
}
