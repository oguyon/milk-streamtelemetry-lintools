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
void perform_cca_float(float *X, long N, long P, float *Y, long Q, int nvec, float **A_vec, float **B_vec);
void perform_pca_float(float *X, long N, long P, int npca, float **Coeffs, float **Modes, const char* out_filename, int xa, int ya);

void print_help(const char *progname) {
    fprintf(stderr, "Usage: %s [options] <nvec> <A.fits> <B.fits>\n", progname);
    fprintf(stderr, "\n");
    fprintf(stderr, "Performs Canonical Correlation Analysis (CCA) between two datasets.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  <nvec>          Number of canonical vectors to compute.\n");
    fprintf(stderr, "  <A.fits>        Input dataset A (3D cube or 2D coefficients).\n");
    fprintf(stderr, "  <B.fits>        Input dataset B (3D cube or 2D coefficients).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -npca <n>       Pre-process inputs with PCA (keep n modes) before CCA.\n");
    fprintf(stderr, "  -ncpu <n>       Set number of CPU threads (OpenBLAS).\n");
    fprintf(stderr, "  -float          Use single precision (float) instead of double.\n");
    fprintf(stderr, "  -shift <n>      Time shift dataset B by n steps (positive or negative).\n");
    fprintf(stderr, "\n");
}

int main(int argc, char *argv[]) {
    int npca = 0;
    int ncpu = 0;
    int use_float = 0;
    int arg_offset = 0;
    int shift = 0;

    // We need to loop over arguments because we might have multiple flags now (-npca, -ncpu, -float)
    // Basic argument parsing loop
    int i = 1;
    while (i < argc) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_help(argv[0]);
            return 0;
        } else if (strcmp(argv[i], "-npca") == 0) {
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
        } else if (strcmp(argv[i], "-float") == 0) {
            use_float = 1;
            i += 1;
            arg_offset += 1;
        } else if (strcmp(argv[i], "-shift") == 0) {
             if (i + 1 >= argc) {
                fprintf(stderr, "Error: -shift requires an argument.\n");
                return 1;
            }
            shift = atoi(argv[i+1]);
            i += 2;
            arg_offset += 2;
        } else {
            break; // Positional argument found
        }
    }

    if (argc - arg_offset != 4) {
         if (argc - arg_offset < 4) {
             fprintf(stderr, "Error: Missing required arguments. Expected 3 positional arguments, found %d.\n", argc - arg_offset - 1);
         } else {
             fprintf(stderr, "Error: Too many arguments.\n");
         }
         print_help(argv[0]);
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
    if (shift != 0) {
        printf("Time shifting second series by %d steps.\n", shift);
    }

    int nvec = atoi(argv[1 + arg_offset]);
    const char *fileA = argv[2 + arg_offset];
    const char *fileB = argv[3 + arg_offset];

    if (nvec < 1) {
        fprintf(stderr, "Error: nvec must be >= 1. Got %d.\n", nvec);
        return 1;
    }

    printf("CCA Configuration:\n");
    printf("  nvec: %d\n", nvec);
    if (npca > 0) printf("  npca: %d\n", npca);
    printf("  File A: %s\n", fileA);
    printf("  File B: %s\n", fileB);

    double *X = NULL, *Y = NULL;
    long Na, Pa, Nb, Qb;
    int xa, ya, xb, yb;
    int naxisA, naxisB;

    if (use_float) {
        float *X_f = NULL, *Y_f = NULL;
        float *X_orig = NULL, *Y_orig = NULL;

        printf("Reading %s...\n", fileA);
        read_fits_float(fileA, &X_f, &Na, &Pa, &xa, &ya, &naxisA);
        X_orig = X_f;
        printf("  Dimensions: %ld samples x %ld pixels (%d x %d), NAXIS=%d\n", Na, Pa, xa, ya, naxisA);

        printf("Reading %s...\n", fileB);
        read_fits_float(fileB, &Y_f, &Nb, &Qb, &xb, &yb, &naxisB);
        Y_orig = Y_f;
        printf("  Dimensions: %ld samples x %ld pixels (%d x %d), NAXIS=%d\n", Nb, Qb, xb, yb, naxisB);

        // Shift and Truncate Logic
        long t_start = (0 > -shift) ? 0 : -shift; // max(0, -shift)
        long t_end = (Na < (Nb - shift)) ? Na : (Nb - shift); // min(Na, Nb - shift)
        long N_new = t_end - t_start;

        if (N_new < 1) {
            fprintf(stderr, "Error: Overlap after shift is %ld (invalid).\n", N_new);
            free(X_orig); free(Y_orig);
            return 1;
        }

        if (N_new < Na || N_new < Nb) {
            printf("  Truncating/Shifting: N_overlap = %ld (Original: %ld, %ld)\n", N_new, Na, Nb);
        }

        // Adjust pointers
        X_f = X_orig + (t_start * Pa);
        Y_f = Y_orig + ((t_start + shift) * Qb);
        long N = N_new;

        int is_coeffs = 0;
        if (naxisA == 2 && naxisB == 2) {
            is_coeffs = 1;
            printf("Detected input as coefficients (2D).\n");
        } else if (naxisA == 3 && naxisB == 3) {
            is_coeffs = 0;
            printf("Detected input as image cubes (3D).\n");
        } else {
            fprintf(stderr, "Error: Mismatched input formats (NAXIS A=%d, B=%d).\n", naxisA, naxisB);
            free(X_orig); free(Y_orig);
            return 1;
        }

        if (nvec > N) {
            fprintf(stderr, "Warning: nvec (%d) > N (%ld). Reducing nvec to N.\n", nvec, N);
            nvec = (int)N;
        }

        float *A_vec_cca = NULL, *B_vec_cca = NULL;

        if (is_coeffs) {
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

            perform_cca_float(X_f, N, Pa, Y_f, Qb, nvec, &A_vec_cca, &B_vec_cca);

            printf("Writing ccaA.fits (Coeffs)...\n");
            write_fits_2d_float("ccaA.fits", A_vec_cca, Pa, nvec);

            printf("Writing ccaB.fits (Coeffs)...\n");
            write_fits_2d_float("ccaB.fits", B_vec_cca, Qb, nvec);

        } else if (npca > 0) {
            if (npca > N) {
                fprintf(stderr, "Warning: npca (%d) > N (%ld). Reducing npca to N.\n", npca, N);
                npca = (int)N;
            }
            if (npca > Pa) npca = (int)Pa;
            if (npca > Qb) npca = (int)Qb;

            printf("Performing PCA on A (keeping %d modes)...\n", npca);
            float *CoeffsA = NULL, *ModesA = NULL;
            perform_pca_float(X_f, N, Pa, npca, &CoeffsA, &ModesA, "pcaA.fits", xa, ya);

            printf("Performing PCA on B (keeping %d modes)...\n", npca);
            float *CoeffsB = NULL, *ModesB = NULL;
            perform_pca_float(Y_f, N, Qb, npca, &CoeffsB, &ModesB, "pcaB.fits", xb, yb);

            float *Wa = NULL, *Wb = NULL;
            if (nvec > npca) {
                 fprintf(stderr, "Warning: nvec (%d) > npca (%d). Reducing nvec to npca.\n", nvec, npca);
                 nvec = npca;
            }

            printf("Performing CCA on PCA coefficients...\n");
            perform_cca_float(CoeffsA, N, npca, CoeffsB, npca, nvec, &Wa, &Wb);

            // Reconstruct
            A_vec_cca = (float *)malloc(nvec * Pa * sizeof(float));
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                        nvec, Pa, npca,
                        1.0f, Wa, npca,
                        ModesA, Pa,
                        0.0f, A_vec_cca, Pa);

            B_vec_cca = (float *)malloc(nvec * Qb * sizeof(float));
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                        nvec, Qb, npca,
                        1.0f, Wb, npca,
                        ModesB, Qb,
                        0.0f, B_vec_cca, Qb);

            free(CoeffsA); free(ModesA);
            free(CoeffsB); free(ModesB);
            free(Wa); free(Wb);

            printf("Writing ccaA.fits...\n");
            write_fits_3d_float("ccaA.fits", A_vec_cca, xa, ya, nvec);

            printf("Writing ccaB.fits...\n");
            write_fits_3d_float("ccaB.fits", B_vec_cca, xb, yb, nvec);

        } else {
            if (nvec > Pa) {
                fprintf(stderr, "Warning: nvec (%d) > Pa (%ld). Reducing nvec to Pa.\n", nvec, Pa);
                nvec = (int)Pa;
            }
            if (nvec > Qb) {
                fprintf(stderr, "Warning: nvec (%d) > Qb (%ld). Reducing nvec to Qb.\n", nvec, Qb);
                nvec = (int)Qb;
            }

            perform_cca_float(X_f, N, Pa, Y_f, Qb, nvec, &A_vec_cca, &B_vec_cca);

            printf("Writing ccaA.fits...\n");
            write_fits_3d_float("ccaA.fits", A_vec_cca, xa, ya, nvec);

            printf("Writing ccaB.fits...\n");
            write_fits_3d_float("ccaB.fits", B_vec_cca, xb, yb, nvec);
        }

        printf("Done.\n");

        free(X_orig); free(Y_orig);
        free(A_vec_cca); free(B_vec_cca);

    } else {
        // Double precision (original code)
        double *X_orig = NULL, *Y_orig = NULL;

        printf("Reading %s...\n", fileA);
        read_fits(fileA, &X, &Na, &Pa, &xa, &ya, &naxisA);
        X_orig = X;
        printf("  Dimensions: %ld samples x %ld pixels (%d x %d), NAXIS=%d\n", Na, Pa, xa, ya, naxisA);

        printf("Reading %s...\n", fileB);
        read_fits(fileB, &Y, &Nb, &Qb, &xb, &yb, &naxisB);
        Y_orig = Y;
        printf("  Dimensions: %ld samples x %ld pixels (%d x %d), NAXIS=%d\n", Nb, Qb, xb, yb, naxisB);

        // Shift and Truncate Logic
        long t_start = (0 > -shift) ? 0 : -shift; // max(0, -shift)
        long t_end = (Na < (Nb - shift)) ? Na : (Nb - shift); // min(Na, Nb - shift)
        long N_new = t_end - t_start;

        if (N_new < 1) {
            fprintf(stderr, "Error: Overlap after shift is %ld (invalid).\n", N_new);
            free(X_orig); free(Y_orig);
            return 1;
        }

        if (N_new < Na || N_new < Nb) {
            printf("  Truncating/Shifting: N_overlap = %ld (Original: %ld, %ld)\n", N_new, Na, Nb);
        }

        // Adjust pointers
        X = X_orig + (t_start * Pa);
        Y = Y_orig + ((t_start + shift) * Qb);
        long N = N_new;

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
            free(X_orig); free(Y_orig);
            return 1;
        }

        if (nvec > N) {
            fprintf(stderr, "Warning: nvec (%d) > N (%ld). Reducing nvec to N.\n", nvec, N);
            nvec = (int)N;
        }

        double *A_vec_cca = NULL, *B_vec_cca = NULL;

        if (is_coeffs) {
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

            printf("Writing ccaA.fits (Coeffs)...\n");
            write_fits_2d("ccaA.fits", A_vec_cca, Pa, nvec);

            printf("Writing ccaB.fits (Coeffs)...\n");
            write_fits_2d("ccaB.fits", B_vec_cca, Qb, nvec);

        } else if (npca > 0) {
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

        free(X_orig); free(Y_orig);
        free(A_vec_cca); free(B_vec_cca);
    }
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

void perform_pca_float(float *X, long N, long P, int npca, float **Coeffs, float **Modes, const char* out_filename, int xa, int ya) {
    long K = (N < P) ? N : P;
    float *S = (float *)malloc(K * sizeof(float));
    float *U = (float *)malloc(N * K * sizeof(float));
    float *Vt = (float *)malloc(K * P * sizeof(float));

    lapack_int info = LAPACKE_sgesdd(LAPACK_ROW_MAJOR, 'S', N, P, X, P,
                                     S, U, K, Vt, P);
    if (info != 0) {
        fprintf(stderr, "LAPACKE_sgesdd failed with error code %d\n", info);
        exit(1);
    }

    *Coeffs = (float *)malloc(N * npca * sizeof(float));
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
            (*Coeffs)[i * npca + 0] *= -1.0f;
        }
        // Flip Vt row 0 (which becomes Mode 0)
        for(long j=0; j<P; j++) {
            Vt[0 * P + j] *= -1.0f;
        }
    }

    *Modes = (float *)malloc(npca * P * sizeof(float));
    for(int i=0; i<npca; i++) {
        for(long j=0; j<P; j++) {
            (*Modes)[i * P + j] = Vt[i * P + j];
        }
    }

    write_fits_3d_float(out_filename, *Modes, xa, ya, npca);

    free(S); free(U); free(Vt);
}

void perform_cca_float(float *X, long N, long P, float *Y, long Q, int nvec, float **A_vec, float **B_vec) {
    long Kx = (N < P) ? N : P;
    long Ky = (N < Q) ? N : Q;

    float *Rx = (float *)malloc(Kx * P * sizeof(float));
    float *Ry = (float *)malloc(Ky * Q * sizeof(float));
    float *tauX = (float *)malloc(Kx * sizeof(float));
    float *tauY = (float *)malloc(Ky * sizeof(float));

    LAPACKE_sgeqrf(LAPACK_ROW_MAJOR, N, P, X, P, tauX);
    memset(Rx, 0, Kx * P * sizeof(float));
    for (long i = 0; i < Kx; i++) {
        for (long j = i; j < P; j++) {
            Rx[i * P + j] = X[i * P + j];
        }
    }
    LAPACKE_sorgqr(LAPACK_ROW_MAJOR, N, Kx, Kx, X, P, tauX);

    LAPACKE_sgeqrf(LAPACK_ROW_MAJOR, N, Q, Y, Q, tauY);
    memset(Ry, 0, Ky * Q * sizeof(float));
    for (long i = 0; i < Ky; i++) {
        for (long j = i; j < Q; j++) {
            Ry[i * Q + j] = Y[i * Q + j];
        }
    }
    LAPACKE_sorgqr(LAPACK_ROW_MAJOR, N, Ky, Ky, Y, Q, tauY);

    float *M = (float *)malloc(Kx * Ky * sizeof(float));
    cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                Kx, Ky, N,
                1.0f, X, P,
                Y, Q,
                0.0f, M, Ky);

    float *S = (float *)malloc(((Kx < Ky) ? Kx : Ky) * sizeof(float));
    float *U = (float *)malloc(Kx * Kx * sizeof(float));
    float *Vt = (float *)malloc(Ky * Ky * sizeof(float));

    lapack_int info = LAPACKE_sgesdd(LAPACK_ROW_MAJOR, 'A', Kx, Ky, M, Ky,
                   S, U, Kx, Vt, Ky);
    if (info != 0) {
        fprintf(stderr, "LAPACKE_sgesdd failed with error code %d\n", info);
        exit(1);
    }

    printf("Canonical Correlations:\n");
    for(int i=0; i<nvec; i++) {
        if (i < ((Kx < Ky) ? Kx : Ky)) {
            printf("  Mode %d: %g\n", i, S[i]);
        }
    }

    *A_vec = (float *)malloc(nvec * P * sizeof(float));
    *B_vec = (float *)malloc(nvec * Q * sizeof(float));

    float *RHS_A = (float *)calloc(P * nvec, sizeof(float));
    for (int j = 0; j < nvec; j++) {
        for (int i = 0; i < Kx; i++) {
            RHS_A[i * nvec + j] = U[i * Kx + j];
        }
    }

    float *Rx_tmp = (float *)malloc(Kx * P * sizeof(float));
    memcpy(Rx_tmp, Rx, Kx * P * sizeof(float));
    float *S_Rx = (float *)malloc(((Kx < P) ? Kx : P) * sizeof(float));
    lapack_int rank_Rx;
    LAPACKE_sgelss(LAPACK_ROW_MAJOR, Kx, P, nvec, Rx_tmp, P, RHS_A, nvec, S_Rx, -1, &rank_Rx);

    for (int j = 0; j < nvec; j++) {
        for (long i = 0; i < P; i++) {
            (*A_vec)[j * P + i] = RHS_A[i * nvec + j];
        }
    }
    free(S_Rx); free(Rx_tmp); free(RHS_A);

    float *RHS_B = (float *)calloc(Q * nvec, sizeof(float));
    for (int j = 0; j < nvec; j++) {
        for (int i = 0; i < Ky; i++) {
            RHS_B[i * nvec + j] = Vt[j * Ky + i];
        }
    }

    float *Ry_tmp = (float *)malloc(Ky * Q * sizeof(float));
    memcpy(Ry_tmp, Ry, Ky * Q * sizeof(float));
    float *S_Ry = (float *)malloc(((Ky < Q) ? Ky : Q) * sizeof(float));
    lapack_int rank_Ry;
    LAPACKE_sgelss(LAPACK_ROW_MAJOR, Ky, Q, nvec, Ry_tmp, Q, RHS_B, nvec, S_Ry, -1, &rank_Ry);

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
