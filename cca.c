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
    int ncpu = 0;
    int arg_offset = 0;

    // We need to loop over arguments because we might have multiple flags now (-npca, -ncpu)
    // Basic argument parsing loop
    int i = 1;
    while (i < argc) {
        if (strcmp(argv[i], "-npca") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: -npca requires an argument.\n");
                return 1;
            }
            npca = atoi(argv[i+1]);
            if (npca < 1) {
                fprintf(stderr, "Error: npca must be >= 1. Got %d.\n", npca);
                return 1;
            }
            i += 2;
            arg_offset += 2;
        } else if (strcmp(argv[i], "-ncpu") == 0) {
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
        } else {
            break; // Positional argument found
        }
    }

    if (argc - arg_offset != 4) {
         fprintf(stderr, "Usage: %s [-npca <n>] [-ncpu <n>] <nvec> <A.fits> <B.fits>\n", argv[0]);
         return 1;
    }

    if (ncpu > 0) {
        openblas_set_num_threads(ncpu);
        printf("Using %d CPU cores.\n", ncpu);
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
    int naxisA, naxisB;

    // Read Data
    printf("Reading %s...\n", fileA);
    read_fits(fileA, &X, &Na, &Pa, &xa, &ya, &naxisA);
    printf("  Dimensions: %ld samples x %ld pixels (%d x %d), NAXIS=%d\n", Na, Pa, xa, ya, naxisA);

    printf("Reading %s...\n", fileB);
    read_fits(fileB, &Y, &Nb, &Qb, &xb, &yb, &naxisB);
    printf("  Dimensions: %ld samples x %ld pixels (%d x %d), NAXIS=%d\n", Nb, Qb, xb, yb, naxisB);

    if (Na != Nb) {
        fprintf(stderr, "Error: Number of samples (N) must match between A and B.\n");
        free(X); free(Y);
        return 1;
    }
    long N = Na;

    // Detect if input is Coefficients (2D) or Image Cube (3D)
    int is_coeffs = 0;
    if (naxisA == 2 && naxisB == 2) {
        is_coeffs = 1;
        printf("Detected input as coefficients (2D).\n");
    } else if (naxisA == 3 && naxisB == 3) {
        is_coeffs = 0;
        printf("Detected input as image cubes (3D).\n");
    } else {
        fprintf(stderr, "Error: Mismatched input formats (NAXIS A=%d, B=%d).\n", naxisA, naxisB);
        free(X); free(Y);
        return 1;
    }

    if (nvec > N) {
        fprintf(stderr, "Warning: nvec (%d) > N (%ld). Reducing nvec to N.\n", nvec, N);
        nvec = (int)N;
    }

    // Center Data?
    // User requested "Do not de-average the input" for PCA.
    // We apply this globally. The first mode (of PCA or CCA?) will likely capture the mean.
    // printf("Centering data...\n");
    // center_columns(X, N, Pa);
    // center_columns(Y, N, Qb);

    double *A_vec_cca = NULL, *B_vec_cca = NULL;

    if (is_coeffs) {
        // Input are coefficients. npca argument should ideally be ignored or matched?
        // If user passed -npca, it's irrelevant here because we don't do PCA on coefficients (usually).
        // Or do we? "Modify the cca program so it takes as input the modal coefficients... The output vectors are then also modal coefficients."
        // This implies CCA on Coeffs.
        if (npca > 0) {
            fprintf(stderr, "Warning: -npca ignored because input is already coefficients.\n");
        }

        if (nvec > Pa) {
            fprintf(stderr, "Warning: nvec (%d) > Pa (%ld). Reducing nvec to Pa.\n", nvec, Pa);
            nvec = (int)Pa;
        }
        if (nvec > Qb) {
            fprintf(stderr, "Warning: nvec (%d) > Qb (%ld). Reducing nvec to Qb.\n", nvec, Qb);
            nvec = (int)Qb;
        }

        perform_cca(X, N, Pa, Y, Qb, nvec, &A_vec_cca, &B_vec_cca);

        // Output is 2D coefficients (nvec x P)
        printf("Writing ccaA.fits (Coeffs)...\n");
        // write_fits_2d: width=P, height=nvec.
        write_fits_2d("ccaA.fits", A_vec_cca, Pa, nvec);

        printf("Writing ccaB.fits (Coeffs)...\n");
        write_fits_2d("ccaB.fits", B_vec_cca, Qb, nvec);

    } else if (npca > 0) {
        // PCA Mode on 3D images
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

        printf("Writing ccaA.fits...\n");
        write_fits_3d("ccaA.fits", A_vec_cca, xa, ya, nvec);

        printf("Writing ccaB.fits...\n");
        write_fits_3d("ccaB.fits", B_vec_cca, xb, yb, nvec);

    } else {
        // Raw CCA on 3D images
        if (nvec > Pa) {
            fprintf(stderr, "Warning: nvec (%d) > Pa (%ld). Reducing nvec to Pa.\n", nvec, Pa);
            nvec = (int)Pa;
        }
        if (nvec > Qb) {
            fprintf(stderr, "Warning: nvec (%d) > Qb (%ld). Reducing nvec to Qb.\n", nvec, Qb);
            nvec = (int)Qb;
        }

        perform_cca(X, N, Pa, Y, Qb, nvec, &A_vec_cca, &B_vec_cca);

        printf("Writing ccaA.fits...\n");
        write_fits_3d("ccaA.fits", A_vec_cca, xa, ya, nvec);

        printf("Writing ccaB.fits...\n");
        write_fits_3d("ccaB.fits", B_vec_cca, xb, yb, nvec);
    }

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

    lapack_int info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'S', N, P, X, P,
                                     S, U, K, Vt, P);
    if (info != 0) {
        fprintf(stderr, "LAPACKE_dgesdd failed with error code %d\n", info);
        exit(1);
    }

    *Coeffs = (double *)malloc(N * npca * sizeof(double));
    for(long i=0; i<N; i++) {
        for(int j=0; j<npca; j++) {
            (*Coeffs)[i * npca + j] = U[i * K + j] * S[j];
        }
    }

    // Check sign of the first mode
    double coeff_mean = 0.0;
    for(long i=0; i<N; i++) {
        coeff_mean += (*Coeffs)[i * npca + 0];
    }
    coeff_mean /= N;

    if (coeff_mean < 0) {
        printf("Flipping sign of Mode 0 to ensure positive average coefficient.\n");
        // Flip Coeffs column 0
        for(long i=0; i<N; i++) {
            (*Coeffs)[i * npca + 0] *= -1.0;
        }
        // Flip Vt row 0 (which becomes Mode 0)
        for(long j=0; j<P; j++) {
            Vt[0 * P + j] *= -1.0;
        }
    }

    *Modes = (double *)malloc(npca * P * sizeof(double));
    for(int i=0; i<npca; i++) {
        for(long j=0; j<P; j++) {
            (*Modes)[i * P + j] = Vt[i * P + j];
        }
    }

    write_fits_3d(out_filename, *Modes, xa, ya, npca);

    free(S); free(U); free(Vt);
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

    lapack_int info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'A', Kx, Ky, M, Ky,
                   S, U, Kx, Vt, Ky);
    if (info != 0) {
        fprintf(stderr, "LAPACKE_dgesdd failed with error code %d\n", info);
        exit(1);
    }

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
    free(M); free(S); free(U); free(Vt);
}
