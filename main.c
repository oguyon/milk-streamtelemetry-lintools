#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fitsio.h>
#include <cblas.h>
#include <lapacke.h>

// Macro to check CFITSIO status
#define CHECK_STATUS(status) if (status) { fits_report_error(stderr, status); exit(status); }

// Function declarations
void read_fits(const char *filename, double **data, long *N, long *P, int *xa, int *ya);
void write_fits(const char *filename, double *data, int xa, int ya, int nvec);
void center_columns(double *data, long N, long P);
void perform_cca(double *X, long N, long P, double *Y, long Q, int nvec, double **A_vec, double **B_vec);
void perform_pca(double *X, long N, long P, int npca, double **Coeffs, double **Modes, const char* out_filename, int xa, int ya);

int main(int argc, char *argv[]) {
    // Parse arguments
    // Usage: milk-streamtelemetry-cca [-npca <n>] <nvec> <A.fits> <B.fits>

    int npca = 0;
    int arg_offset = 0;

    if (argc >= 6 && strcmp(argv[1], "-npca") == 0) {
        npca = atoi(argv[2]);
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
        // Also check P
        if (npca > Pa) {
             npca = (int)Pa;
        }
        if (npca > Qb) {
             npca = (int)Qb;
        }

        printf("Performing PCA on A (keeping %d modes)...\n", npca);
        double *CoeffsA = NULL, *ModesA = NULL;
        perform_pca(X, N, Pa, npca, &CoeffsA, &ModesA, "pcaA.fits", xa, ya);

        printf("Performing PCA on B (keeping %d modes)...\n", npca);
        double *CoeffsB = NULL, *ModesB = NULL;
        perform_pca(Y, N, Qb, npca, &CoeffsB, &ModesB, "pcaB.fits", xb, yb);

        // Perform CCA on Coefficients
        // CoeffsA: N x npca
        // CoeffsB: N x npca
        // CCA returns vectors in PCA space: Wa (npca x nvec), Wb (npca x nvec)

        double *Wa = NULL, *Wb = NULL;

        if (nvec > npca) {
             fprintf(stderr, "Warning: nvec (%d) > npca (%d). Reducing nvec to npca.\n", nvec, npca);
             nvec = npca;
        }

        printf("Performing CCA on PCA coefficients...\n");
        perform_cca(CoeffsA, N, npca, CoeffsB, npca, nvec, &Wa, &Wb);

        // Reconstruct Spatial Canonical Vectors
        // A_vec = ModesA^T * Wa
        // ModesA is npca x Pa (stored row major: Mode0, Mode1...).
        // We want A_vec (nvec x Pa) for writing.
        // Wait. ModesA was written as images. Size (xa, ya, npca).
        // In perform_pca, we store ModesA as npca x Pa row-major (Mode 0 is first row).
        // Wa is npca x nvec?
        // perform_cca output: "A_vec" (P x nvec) stored as nvec x P in memory for writing.
        // Here, "P" is `npca`. So Wa is nvec x npca.
        // Wa (nvec x npca) * ModesA (npca x Pa) = A_vec (nvec x Pa).
        // This gives us the nvec spatial vectors. Correct.

        // Dimensions:
        // Wa: nvec x npca.
        // ModesA: npca x Pa.
        // Output: nvec x Pa.

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
        // Original Raw CCA

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

    // Write Outputs
    printf("Writing ccaA.fits...\n");
    write_fits("ccaA.fits", A_vec_cca, xa, ya, nvec);

    printf("Writing ccaB.fits...\n");
    write_fits("ccaB.fits", B_vec_cca, xb, yb, nvec);

    printf("Done.\n");

    free(X);
    free(Y);
    free(A_vec_cca);
    free(B_vec_cca);

    return 0;
}

// ... (read_fits, write_fits, center_columns remain same)
// ... (perform_cca remains same)

// New function: perform_pca
// Computes SVD of X (N x P).
// Returns Coeffs (N x npca) and Modes (npca x P).
// Writes Modes to FITS.
void perform_pca(double *X, long N, long P, int npca, double **Coeffs, double **Modes, const char* out_filename, int xa, int ya) {
    // SVD: X = U S Vt
    // We want first npca components.
    // Economy SVD returns min(N, P) singular values/vectors.
    long K = (N < P) ? N : P;

    // Allocate S, U, Vt
    double *S = (double *)malloc(K * sizeof(double));
    double *U = (double *)malloc(N * K * sizeof(double)); // Row major: N x K
    double *Vt = (double *)malloc(K * P * sizeof(double)); // Row major: K x P
    double *superb = (double *)malloc((K - 1) * sizeof(double));

    // X is modified by dgesvd
    // Make a copy if needed? X is passed from main, can be destroyed?
    // perform_pca is called on X. X is not used afterwards in main if PCA path is taken (except to be freed).
    // But dgesvd destroys input 'A'.

    LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'S', 'S', N, P, X, P,
                   S, U, K, Vt, P, superb);

    // Coeffs = U * S (first npca columns)
    // Coeffs size: N x npca.
    *Coeffs = (double *)malloc(N * npca * sizeof(double));

    // Fill Coeffs
    for(long i=0; i<N; i++) {
        for(int j=0; j<npca; j++) {
            (*Coeffs)[i * npca + j] = U[i * K + j] * S[j];
        }
    }

    // Modes = Vt (first npca rows)
    // Modes size: npca x P.
    *Modes = (double *)malloc(npca * P * sizeof(double));

    // Vt is K x P. Copy first npca rows.
    for(int i=0; i<npca; i++) {
        for(long j=0; j<P; j++) {
            (*Modes)[i * P + j] = Vt[i * P + j];
        }
    }

    // Write Modes to FITS
    // write_fits expects (xa, ya, nvec). Data layout: nvec x (xa*ya).
    // Our Modes matrix matches exactly: npca x P.
    write_fits(out_filename, *Modes, xa, ya, npca);

    free(S); free(U); free(Vt); free(superb);
}


void read_fits(const char *filename, double **data, long *N, long *P, int *xa, int *ya) {
    fitsfile *fptr;
    int status = 0;
    int naxis;
    long naxes[3] = {0, 0, 0};

    fits_open_file(&fptr, filename, READONLY, &status);
    CHECK_STATUS(status);

    fits_get_img_dim(fptr, &naxis, &status);
    CHECK_STATUS(status);

    if (naxis != 3) {
        fprintf(stderr, "Error: %s is not a 3D FITS cube (NAXIS=%d).\n", filename, naxis);
        fits_close_file(fptr, &status);
        exit(1);
    }

    fits_get_img_size(fptr, 3, naxes, &status);
    CHECK_STATUS(status);

    *xa = (int)naxes[0];
    *ya = (int)naxes[1];
    *N = naxes[2];
    *P = (*xa) * (*ya);

    long total_elements = (*N) * (*P);
    *data = (double *)malloc(total_elements * sizeof(double));
    if (*data == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for %ld doubles.\n", total_elements);
        exit(1);
    }

    // Read all elements
    // FITS stores data in FORTRAN order (first index varies fastest).
    // C stores row-major.
    // However, fits_read_img reads sequentially from the file.
    // File structure:
    // Image 0: Row 0 (xa pixels), Row 1 (xa pixels)... Row ya-1.
    // Image 1: ...
    // So the buffer will contain Image 0, then Image 1, etc.
    //
    // We want X to be N x P.
    // Row 0 of X: Image 0 flattened.
    // Row 1 of X: Image 1 flattened.
    // This matches the file layout perfectly.
    // X[i * P + j] is the j-th pixel of the i-th image.

    long fpixel[3] = {1, 1, 1};
    fits_read_pix(fptr, TDOUBLE, fpixel, total_elements, NULL, *data, NULL, &status);
    CHECK_STATUS(status);

    fits_close_file(fptr, &status);
    CHECK_STATUS(status);
}

void write_fits(const char *filename, double *data, int xa, int ya, int nvec) {
    fitsfile *fptr;
    int status = 0;
    long naxes[3];
    naxes[0] = xa;
    naxes[1] = ya;
    naxes[2] = nvec;

    // Delete file if exists
    remove(filename);

    fits_create_file(&fptr, filename, &status);
    CHECK_STATUS(status);

    fits_create_img(fptr, DOUBLE_IMG, 3, naxes, &status);
    CHECK_STATUS(status);

    long total_elements = (long)xa * ya * nvec;
    long fpixel[3] = {1, 1, 1};
    fits_write_pix(fptr, TDOUBLE, fpixel, total_elements, data, &status);
    CHECK_STATUS(status);

    fits_close_file(fptr, &status);
    CHECK_STATUS(status);
}

void center_columns(double *data, long N, long P) {
    // Data is N x P (N rows, P columns).
    // We need to subtract mean from each column.
    // X[i * P + j]

    double *means = (double *)calloc(P, sizeof(double));
    if (!means) exit(1);

    // Compute means
    for (long i = 0; i < N; i++) {
        for (long j = 0; j < P; j++) {
            means[j] += data[i * P + j];
        }
    }
    for (long j = 0; j < P; j++) {
        means[j] /= N;
    }

    // Subtract means
    for (long i = 0; i < N; i++) {
        for (long j = 0; j < P; j++) {
            data[i * P + j] -= means[j];
        }
    }
    free(means);
}

void perform_cca(double *X, long N, long P, double *Y, long Q, int nvec, double **A_vec, double **B_vec) {
    // 1. QR Decomposition of X -> Qx Rx
    long Kx = (N < P) ? N : P;
    long Ky = (N < Q) ? N : Q;

    // Allocate Rx and Ry
    // Rx: Kx x P.
    double *Rx = (double *)malloc(Kx * P * sizeof(double));
    double *Ry = (double *)malloc(Ky * Q * sizeof(double));

    // Tau arrays for QR
    double *tauX = (double *)malloc(Kx * sizeof(double));
    double *tauY = (double *)malloc(Ky * sizeof(double));

    // QR of X
    LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, N, P, X, P, tauX);

    // Extract Rx
    memset(Rx, 0, Kx * P * sizeof(double));
    for (long i = 0; i < Kx; i++) {
        for (long j = i; j < P; j++) {
            Rx[i * P + j] = X[i * P + j];
        }
    }

    // Form Qx
    LAPACKE_dorgqr(LAPACK_ROW_MAJOR, N, Kx, Kx, X, P, tauX);

    // QR of Y
    LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, N, Q, Y, Q, tauY);

    // Extract Ry
    memset(Ry, 0, Ky * Q * sizeof(double));
    for (long i = 0; i < Ky; i++) {
        for (long j = i; j < Q; j++) {
            Ry[i * Q + j] = Y[i * Q + j];
        }
    }

    // Form Qy
    LAPACKE_dorgqr(LAPACK_ROW_MAJOR, N, Ky, Ky, Y, Q, tauY);

    // Compute M = Qx^T * Qy
    double *M = (double *)malloc(Kx * Ky * sizeof(double));

    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                Kx, Ky, N,
                1.0, X, P,
                Y, Q,
                0.0, M, Ky);

    // SVD of M -> U S V^T
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

    // Calculate Canonical Vectors
    *A_vec = (double *)malloc(nvec * P * sizeof(double));
    *B_vec = (double *)malloc(nvec * Q * sizeof(double));

    // For A: Solve Rx * A = U
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

    free(S_Rx);
    free(Rx_tmp);
    free(RHS_A);

    // Repeat for B: Solve Ry * B = V (columns of V are rows of Vt)
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

    free(S_Ry);
    free(Ry_tmp);
    free(RHS_B);

    free(Rx); free(Ry);
    free(tauX); free(tauY);
    free(M); free(S); free(U); free(Vt); free(superb);
}
