#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "common.h"

int main() {
    int N = 10;
    int K = 5;
    int width = 10;
    int height = 10;
    int P = width * height;

    double *frames = (double *)malloc(N * P * sizeof(double));
    double *modes = (double *)malloc(K * P * sizeof(double));

    // Create random frames
    for (int i = 0; i < N * P; i++) {
        frames[i] = ((double)rand() / RAND_MAX);
    }

    // Create random orthonormal modes (simplified, just random here)
    // For proper test, we should orthogonalize, but random is fine for functional test
    for (int i = 0; i < K * P; i++) {
        modes[i] = ((double)rand() / RAND_MAX);
    }

    write_fits_3d("frames.fits", frames, width, height, N);
    write_fits_3d("modes.fits", modes, width, height, K);

    free(frames);
    free(modes);
    return 0;
}
