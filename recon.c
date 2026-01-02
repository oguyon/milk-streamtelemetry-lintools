#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cblas.h>
#include "common.h"

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <modes.fits> <coeffs.fits> <output.fits>\n", argv[0]);
        return 1;
    }

    const char *modes_file = argv[1];
    const char *coeffs_file = argv[2];
    const char *outfile = argv[3];

    printf("Reconstruction Configuration:\n");
    printf("  Modes: %s\n", modes_file);
    printf("  Coeffs: %s\n", coeffs_file);
    printf("  Output: %s\n", outfile);

    double *Modes = NULL, *Coeffs = NULL;
    long N_modes, P_modes; // Modes is npca x P
    int xa, ya;
    long N_coeffs, P_coeffs; // Coeffs is N x npca

    printf("Reading Modes...\n");
    // Modes file is 3D (xa, ya, npca) or maybe 2D if npca=1?
    // read_fits handles both.
    // N will be naxis3 (or naxis2 if 2D), P will be xa*ya.
    // If modes is (xa, ya, npca), then N_modes = npca, P_modes = xa*ya.
    read_fits(modes_file, &Modes, &N_modes, &P_modes, &xa, &ya);
    printf("  Modes Dim: %ld modes x %ld pixels (%d x %d)\n", N_modes, P_modes, xa, ya);

    printf("Reading Coeffs...\n");
    // Coeffs file is 2D (npca, N).
    // read_fits for 2D: N=height, P=width.
    // We wrote it as width=npca, height=N.
    // So N_coeffs = N, P_coeffs = npca.
    int dummy_xa, dummy_ya;
    read_fits(coeffs_file, &Coeffs, &N_coeffs, &P_coeffs, &dummy_xa, &dummy_ya);
    printf("  Coeffs Dim: %ld samples x %ld coeffs\n", N_coeffs, P_coeffs);

    if (N_modes != P_coeffs) {
        fprintf(stderr, "Error: Dimension mismatch. Modes has %ld modes, Coeffs has %ld coeffs.\n", N_modes, P_coeffs);
        free(Modes); free(Coeffs);
        return 1;
    }

    long npca = N_modes;
    long N = N_coeffs;
    long P = P_modes;

    // Reconstruction: X = Coeffs * Modes
    // Coeffs: N x npca
    // Modes: npca x P
    // X: N x P

    double *X = (double *)malloc(N * P * sizeof(double));
    if (!X) {
        fprintf(stderr, "Error: Memory allocation failed for output.\n");
        return 1;
    }

    printf("Reconstructing (Matrix Multiplication)...\n");
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                N, P, npca,
                1.0, Coeffs, npca,
                Modes, P,
                0.0, X, P);

    printf("Writing Output to %s...\n", outfile);
    // X is N x P.
    // write_fits_3d expects (xa, ya, z_dim).
    // Here z_dim = N.
    write_fits_3d(outfile, X, xa, ya, N);

    printf("Done.\n");

    free(Modes); free(Coeffs); free(X);
    return 0;
}
