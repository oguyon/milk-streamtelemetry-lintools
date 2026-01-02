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

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <nvec> <A.fits> <B.fits>\n", argv[0]);
        return 1;
    }

    int nvec = atoi(argv[1]);
    const char *fileA = argv[2];
    const char *fileB = argv[3];

    printf("CCA Configuration:\n");
    printf("  nvec: %d\n", nvec);
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

    if (nvec > Pa) {
        fprintf(stderr, "Warning: nvec (%d) > Pa (%ld). Reducing nvec to Pa.\n", nvec, Pa);
        nvec = (int)Pa;
    }
    if (nvec > Qb) {
        fprintf(stderr, "Warning: nvec (%d) > Qb (%ld). Reducing nvec to Qb.\n", nvec, Qb);
        nvec = (int)Qb;
    }

    // Center Data
    printf("Centering data...\n");
    center_columns(X, N, Pa);
    center_columns(Y, N, Qb);

    // Perform CCA
    printf("Performing CCA...\n");
    double *A_vec = NULL, *B_vec = NULL;

    // We expect output vectors to be (xa*ya) x nvec and (xb*yb) x nvec.
    // In memory, we'll store them such that they can be written as (xa, ya, nvec) fits cube.
    // FITS 3D: (xa, ya, nvec). In C flat array: nvec * (xa*ya).
    // So perform_cca should allocate and fill these.

    perform_cca(X, N, Pa, Y, Qb, nvec, &A_vec, &B_vec);

    // Write Outputs
    printf("Writing ccaA.fits...\n");
    write_fits("ccaA.fits", A_vec, xa, ya, nvec);

    printf("Writing ccaB.fits...\n");
    write_fits("ccaB.fits", B_vec, xb, yb, nvec);

    printf("Done.\n");

    free(X);
    free(Y);
    free(A_vec);
    free(B_vec);

    return 0;
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

// Perform CCA
// Input: X (N x P), Y (N x Q)
// Output: A_vec (P x nvec), B_vec (Q x nvec)
// Stored as flattened arrays. A_vec layout:
// We want to write A_vec as (xa, ya, nvec).
// FITS writer expects: Vector 0 (xa*ya), Vector 1 (xa*ya)...
// So A_vec should be nvec x P in memory (row-major where rows are vectors).
// But for calculation, we compute columns of A.
// So we will compute P x nvec matrix, but we need to transpose it for writing if we want
// the first vector to be contiguous.
// Wait. write_fits takes (xa, ya, nvec).
// fits_write_pix writes linearly.
// The file should have:
// Slice 0 (naxis1=xa, naxis2=ya) -> Vector 0.
// Slice 1 -> Vector 1.
// So in memory, data should be:
// [Vector 0 pixels] [Vector 1 pixels] ...
// So `data` passed to write_fits must be `nvec` blocks of `P` pixels.
// Effectively an `nvec x P` matrix stored row-major.
//
// Our algebra gives `a` as columns.
// So we get `P x nvec`. We need to transpose to `nvec x P` before writing?
// Yes.

void perform_cca(double *X, long N, long P, double *Y, long Q, int nvec, double **A_vec, double **B_vec) {
    // 1. QR Decomposition of X -> Qx Rx
    // X is N x P.
    // Use LAPACKE_dgeqrf.
    // If N < P, R is N x P (trapezoidal).
    // If N >= P, R is P x P (triangular) + zeros.

    // We work with economy size QR if possible, but LAPACK stores Q as reflectors in A.
    // We need to extract R first, then form Q.

    // Allocate memory for R parts.
    // Rx needs to be stored.
    // If N >= P, Rx is P x P.
    // If N < P, Rx is N x P.
    // We will use Rx later to solve Rx * a = u.
    // Let Kx = min(N, P). Rx will be Kx x P (upper trapezoidal).
    // Actually, if N >= P, Rx is P x P upper triangular.
    // If N < P, Rx is N x P upper trapezoidal.

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
    // X is overwritten by reflectors and R diagonal/upper.
    LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, N, P, X, P, tauX);

    // Extract Rx
    // Rx is the upper trapezoidal part of X (first Kx rows).
    memset(Rx, 0, Kx * P * sizeof(double));
    for (long i = 0; i < Kx; i++) {
        for (long j = i; j < P; j++) {
            Rx[i * P + j] = X[i * P + j];
        }
    }

    // Form Qx
    // We only need the first Kx columns of Qx?
    // Wait. CCA involves Qx^T Qy.
    // If we use economy QR, Qx is N x Kx.
    // LAPACKE_dorgqr forms Q.
    // We overwrite X with Qx.
    // We want Qx to be N x Kx.
    // So we call dorgqr with m=N, n=Kx, k=Kx.
    // Note: X was N x P.
    // If we want Qx in place, we need N x Kx.
    // If P > Kx (i.e. P > N), then X is large enough to hold N x N Q? No.
    // X is N x P. If P >= N, Kx = N. Qx is N x N. X fits it.
    // If P < N, Kx = P. Qx is N x P. X fits it.
    // So X buffer is large enough for Qx.

    LAPACKE_dorgqr(LAPACK_ROW_MAJOR, N, Kx, Kx, X, P, tauX);
    // Now X contains Qx (N x Kx) stored in N x P buffer (strided).
    // Wait, dorgqr on N x Kx outputs N x Kx matrix.
    // If lda=P, it writes it into the buffer with stride P.
    // Correct.

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
    // Now Y contains Qy (N x Ky) in N x Q buffer.

    // Compute M = Qx^T * Qy
    // M is Kx x Ky.
    double *M = (double *)malloc(Kx * Ky * sizeof(double));

    // Matrix multiplication: M = X^T * Y
    // X is (N x Kx) with stride P.
    // Y is (N x Ky) with stride Q.
    // We use cblas_dgemm.
    // Op(A)=Trans, Op(B)=NoTrans.
    // M = alpha * A^T * B + beta * C
    // A=X, B=Y, C=M.
    // Dimensions: (Kx x N) * (N x Ky) = Kx x Ky.

    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                Kx, Ky, N,
                1.0, X, P,
                Y, Q,
                0.0, M, Ky);

    // SVD of M -> U S V^T
    // M is Kx x Ky.
    // We need U (Kx x Kx) and V^T (Ky x Ky).
    // Actually we only need first `nvec` columns of U and V.
    // We can use LAPACKE_dgesdd or dgesvd.

    double *S = (double *)malloc(((Kx < Ky) ? Kx : Ky) * sizeof(double));
    double *U = (double *)malloc(Kx * Kx * sizeof(double));
    double *Vt = (double *)malloc(Ky * Ky * sizeof(double));
    double *superb = (double *)malloc((((Kx < Ky) ? Kx : Ky) - 1) * sizeof(double)); // for dgesvd

    LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', Kx, Ky, M, Ky,
                   S, U, Kx, Vt, Ky, superb);

    printf("Canonical Correlations:\n");
    for(int i=0; i<nvec; i++) {
        // nvec is already checked to be <= N, but we should also check if i < min(Kx, Ky)
        // because S has size min(Kx, Ky).
        if (i < ((Kx < Ky) ? Kx : Ky)) {
            printf("  Mode %d: %g\n", i, S[i]);
        }
    }

    // Calculate Canonical Vectors
    // Solve Rx * A = U (first nvec columns)
    // Solve Ry * B = V (first nvec columns of V = rows of Vt^T... wait)
    // SVD: M = U * S * Vt.
    // Canonical vectors for X are related to U.
    // Canonical vectors for Y are related to V.
    // Columns of U are u_i.
    // Rows of Vt are v_i^T. So columns of V are rows of Vt.
    // V_matrix (Ky x Ky) = Vt^T.
    // We need first nvec columns of U and V.

    // We need to output nvec vectors.
    // Output A_out is (nvec x P) [vector 0; vector 1; ...]
    // But we solve for columns.

    // Let's process A vectors.
    // We have Rx (Kx x P).
    // We want a (P x 1) such that Rx * a = u (Kx x 1).
    // If Kx = P (N >= P), Rx is square upper triangular.
    // Use dtrtrs.
    // If Kx = N (N < P), Rx is N x P. Underdetermined.
    // Use min norm solution.
    // Min norm solution to R a = u is a = R^T (R R^T)^-1 u.
    // Or simpler: a = R^T * z, where (R R^T) z = u.
    // Since R comes from QR of X, R R^T = (Q^T X) (Q^T X)^T = Q^T X X^T Q.
    // Wait. R is upper trapezoidal.
    // Using `dgels` with transpose?
    // Solve R * a = u.
    // `dgels` solves min ||b - A x|| or min ||x|| s.t. A x = b.
    // We want min ||a|| s.t. Rx * a = u.
    // Routine `dgels` with 'T' (transpose) solves min ||y|| s.t. A^T y = b ?? No.
    // `dgels` for underdetermined system A x = b (A is m x n, m < n):
    // Finds min norm solution.
    // We pass A = Rx (Kx x P). B = u (Kx x 1).
    // Since Rx is row major? `dgels` modifies A.
    // We need to be careful.

    // We'll process nvec vectors.
    // A_out allocation: nvec * P.
    *A_vec = (double *)malloc(nvec * P * sizeof(double));
    *B_vec = (double *)malloc(nvec * Q * sizeof(double));

    // For A:
    // We solve for multiple RHS if possible, but nvec might be small.
    // Also Rx is modified by dgels? dgels overwrites A with QR factors.
    // We can't reuse Rx easily unless we copy it.
    // Or we solve for all nvec columns at once.
    // Matrix U_nvec (Kx x nvec) contains the first nvec columns of U.
    // We want A_sol (P x nvec) such that Rx * A_sol = U_nvec.

    // Construct RHS matrix for A.
    // U is Kx x Kx (row major).
    // We want first nvec columns.
    // Copy into a new buffer.
    // If we use dgels, B (on input) is max(M, N) x NRHS.
    // Rx is Kx x P.
    // If Kx < P (underdetermined): B should be size P x nvec.
    // Input B: first Kx rows are RHS.
    // Output B: P rows are solution.

    double *RHS_A = (double *)calloc(P * nvec, sizeof(double));
    // Copy first nvec columns of U into RHS_A.
    // RHS_A is P x nvec (stored column-wise? No, LAPACK expects column-major B usually?)
    // LAPACKE_dgels expects B in layout specified.
    // If LAPACK_ROW_MAJOR: B is LDB x NRHS.
    // Wait. If Row Major, B is dim_rhs_A x NRHS ??
    // No. B is usually max(M,N) x NRHS.
    // In Row Major: B has dimensions (LDB, NRHS) where LDB >= NRHS?
    // Let's check docs.
    // LAPACKE_dgels(int matrix_layout, char trans, lapack_int m, lapack_int n, lapack_int nrhs, double* a, lapack_int lda, double* b, lapack_int ldb);
    // Solves minimize ||b - A*x|| or A*x = b.
    // If A is m x n.
    // If m < n (underdetermined):
    // b must be size max(m,n) x nrhs. i.e. n x nrhs.
    // On entry, first m rows contain RHS.
    // On exit, n rows contain solution.
    //
    // Layout:
    // If Row Major:
    // A is m x n.
    // B is max(m,n) x nrhs.
    // LDB is stride of B. In Row Major, LDB is distance between rows?
    // No, LDB is usually column count?
    // In Row Major, B[i][j] is B[i * LDB + j].
    // So LDB must be >= NRHS.
    // Correct.

    // Fill RHS_A with U columns.
    // U is Kx x Kx.
    // RHS_A is P x nvec.
    // We want RHS_A[i, j] = U[i, j] for i < Kx, 0 for i >= Kx.
    for (int j = 0; j < nvec; j++) {
        for (int i = 0; i < Kx; i++) {
            RHS_A[i * nvec + j] = U[i * Kx + j];
        }
    }

    // Copy Rx to a temp buffer because dgelss destroys it.
    double *Rx_tmp = (double *)malloc(Kx * P * sizeof(double));
    memcpy(Rx_tmp, Rx, Kx * P * sizeof(double));

    // Solve Rx * A = U using dgelss (SVD based least squares) to handle rank deficiency
    // Rx is Kx x P.
    // dgelss args: m, n, nrhs, a, lda, b, ldb, s, rcond, rank

    double *S_Rx = (double *)malloc(((Kx < P) ? Kx : P) * sizeof(double));
    lapack_int rank_Rx;
    LAPACKE_dgelss(LAPACK_ROW_MAJOR, Kx, P, nvec, Rx_tmp, P, RHS_A, nvec, S_Rx, -1, &rank_Rx);

    // Solution is in RHS_A (P x nvec).
    // We need to transpose it to A_vec (nvec x P) for writing.
    for (int j = 0; j < nvec; j++) {
        for (long i = 0; i < P; i++) {
            (*A_vec)[j * P + i] = RHS_A[i * nvec + j];
        }
    }

    free(S_Rx);
    free(Rx_tmp);
    free(RHS_A);

    // Repeat for B
    // Solve Ry * B_sol = V
    // Ry is Ky x Q.
    // V's columns are rows of Vt.
    // V[i, j] = Vt[j, i].

    double *RHS_B = (double *)calloc(Q * nvec, sizeof(double));

    // Fill RHS_B with V columns.
    // Vt is Ky x Ky.
    // We want first nvec columns of V => first nvec rows of Vt.
    // Vt is row major.
    // V_col_j is Vt_row_j.
    // RHS_B[i, j] = V[i, j] = Vt[j, i].
    // We want j-th eigenvector (j < nvec).
    // The j-th singular vector for Y is the j-th row of Vt (as returned by SVD: M = U S Vt).
    // Wait. SVD of M = U S Vt.
    // U columns are left singular vectors.
    // Vt rows are right singular vectors.
    // So the j-th pair is (U_col_j, Vt_row_j).
    // Correct.
    // So for RHS_B, we want the j-th column to be Vt_row_j.
    // RHS_B[i, j] = Vt[j, i].
    // Wait. Vt is Ky x Ky.
    // Vt_row_j has Ky elements.
    // So i goes from 0 to Ky-1.

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

    // Transpose to B_vec
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
