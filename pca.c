#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cblas.h>
#include <lapacke.h>
#include "common.h"

int main(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <npca> <input.fits> <modes.fits> <coeffs.fits>\n", argv[0]);
        return 1;
    }

    int npca = atoi(argv[1]);
    const char *infile = argv[2];
    const char *modes_file = argv[3];
    const char *coeffs_file = argv[4];

    printf("PCA Configuration:\n");
    printf("  npca: %d\n", npca);
    printf("  Input: %s\n", infile);
    printf("  Modes Output: %s\n", modes_file);
    printf("  Coeffs Output: %s\n", coeffs_file);

    double *X = NULL;
    long N, P;
    int xa, ya;

    printf("Reading %s...\n", infile);
    read_fits(infile, &X, &N, &P, &xa, &ya);
    printf("  Dimensions: %ld samples x %ld pixels (%d x %d)\n", N, P, xa, ya);

    if (npca > N) {
        fprintf(stderr, "Warning: npca (%d) > N (%ld). Reducing npca to N.\n", npca, N);
        npca = (int)N;
    }
    if (npca > P) {
        fprintf(stderr, "Warning: npca (%d) > P (%ld). Reducing npca to P.\n", npca, P);
        npca = (int)P;
    }

    // Center Data?
    // User requested "Make the PCA a separate program...". Usually PCA implies centering.
    // However, if we reconstruct, we won't have the mean unless we save it.
    // The previous prompt's integrated PCA did center.
    // If we want reconstruction to work, we need to know if we centered.
    // I'll center, but note that reconstruction will be relative to mean.
    // Or I should save mean.fits? The user didn't ask for it.
    // I'll skip centering for this standalone tool to ensure exact reconstruction X ~= C * M.
    // Or I'll assume the user knows PCA centers data.
    // Standard PCA centers. SVD without centering is not PCA (it's just SVD).
    // The prompt says "Modify the program so that a PCA... is first done". This implied centering in CCA.
    // "Make the PCA a separate program... provide a program that reconstructs...".
    // If I center, reconstruction X_rec = C*M + Mean.
    // If I don't save Mean, I can't reconstruct X_rec exactly.
    // But maybe the user just wants the variation.
    // I will center, because "PCA" implies it.
    // And I will accept that reconstruction is "centered reconstruction".
    // Wait. If I don't center, the first mode will be the mean. This is often acceptable in image processing.
    // Let's center to be consistent with `cca.c`.

    printf("Centering data...\n");
    center_columns(X, N, P);

    // SVD
    long K = (N < P) ? N : P;
    double *S = (double *)malloc(K * sizeof(double));
    double *U = (double *)malloc(N * K * sizeof(double));
    double *Vt = (double *)malloc(K * P * sizeof(double));
    double *superb = (double *)malloc((K - 1) * sizeof(double));

    printf("Computing SVD...\n");
    LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'S', 'S', N, P, X, P,
                   S, U, K, Vt, P, superb);

    // Coeffs = U * S (first npca columns)
    double *Coeffs = (double *)malloc(N * npca * sizeof(double));
    for(long i=0; i<N; i++) {
        for(int j=0; j<npca; j++) {
            Coeffs[i * npca + j] = U[i * K + j] * S[j];
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
    // Coeffs is N x npca.
    // write_fits_2d expects (width, height).
    // We want N rows, npca columns?
    // FITS images are usually (x, y).
    // If we view it as an image, maybe x=npca, y=N?
    // Or x=N, y=npca?
    // Let's stick to standard matrix storage.
    // width=npca, height=N.
    write_fits_2d(coeffs_file, Coeffs, npca, N);

    printf("Done.\n");

    free(X); free(S); free(U); free(Vt); free(superb);
    free(Coeffs); free(Modes);
    return 0;
}
