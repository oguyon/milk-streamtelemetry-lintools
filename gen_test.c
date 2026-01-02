#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "common.h"

int main() {
    int N = 10;
    int width = 10;
    int height = 10;
    int P = width * height;

    double *data = (double *)malloc(N * P * sizeof(double));

    // Create some synthetic data
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < P; j++) {
            data[i * P + j] = (double)(i + j) + ((double)rand() / RAND_MAX);
        }
    }

    write_fits_3d("A.fits", data, width, height, N);

    // Create B slightly different
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < P; j++) {
            data[i * P + j] = (double)(i - j) + ((double)rand() / RAND_MAX);
        }
    }
    write_fits_3d("B.fits", data, width, height, N);

    free(data);
    return 0;
}
