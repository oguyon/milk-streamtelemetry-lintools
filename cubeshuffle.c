#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "common.h"

void print_help(const char *progname) {
    fprintf(stderr, "Usage: %s <input.fits> <output.fits>\n", progname);
    fprintf(stderr, "\n");
    fprintf(stderr, "Randomly shuffles the slices (time axis) of a 3D FITS cube.\n");
    fprintf(stderr, "Used for testing and destroying temporal correlations.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  <input.fits>    Input 3D FITS cube.\n");
    fprintf(stderr, "  <output.fits>   Output shuffled 3D FITS cube.\n");
    fprintf(stderr, "\n");
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Error: Incorrect number of arguments. Expected input and output filenames.\n");
        print_help(argv[0]);
        return 1;
    }

    if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        print_help(argv[0]);
        return 0;
    }

    const char *infile = argv[1];
    const char *outfile = argv[2];

    double *data = NULL;
    long N, P;
    int xa, ya;
    int naxis;

    printf("Reading %s...\n", infile);
    read_fits(infile, &data, &N, &P, &xa, &ya, &naxis);
    printf("  Dimensions: %ld samples x %ld pixels (%d x %d)\n", N, P, xa, ya);

    if (N < 2) {
        fprintf(stderr, "Warning: Only %ld sample(s). Shuffling meaningless.\n", N);
    }

    // Seed random number generator
    srand((unsigned int)time(NULL));

    // Fisher-Yates Shuffle
    // We want to shuffle the N samples. Each sample is P doubles.
    // data is N x P (row-major).
    // Slice i is at &data[i*P].

    double *temp_slice = (double *)malloc(P * sizeof(double));
    if (!temp_slice) {
        fprintf(stderr, "Error: Memory allocation failed for temp slice.\n");
        free(data);
        return 1;
    }

    printf("Shuffling %ld slices...\n", N);
    for (long i = N - 1; i > 0; i--) {
        // Pick a random index j such that 0 <= j <= i
        long j = (long)(rand() / (RAND_MAX + 1.0) * (i + 1));

        if (i != j) {
            // Swap slice i and slice j
            double *slice_i = &data[i * P];
            double *slice_j = &data[j * P];

            memcpy(temp_slice, slice_i, P * sizeof(double));
            memcpy(slice_i, slice_j, P * sizeof(double));
            memcpy(slice_j, temp_slice, P * sizeof(double));
        }
    }

    free(temp_slice);

    printf("Writing %s...\n", outfile);
    // Determine output format. If input was 3D, write 3D.
    // read_fits sets naxis.
    if (naxis == 3) {
        write_fits_3d(outfile, data, xa, ya, (int)N);
    } else {
        // Assume 2D (matrix form). write_fits_2d uses (width, height).
        // width = P, height = N.
        write_fits_2d(outfile, data, P, N);
    }

    printf("Done.\n");

    free(data);
    return 0;
}
