#include "common.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <input.fits>\n", argv[0]);
        return 1;
    }

    const char *input_file = argv[1];

    double *data = NULL;
    long rows, cols;
    int xa, ya, naxis;

    // read_fits returns flattened array.
    // Dimensions N (height/rows), P (width/cols).
    // For 2D: xa=width, N=height.
    read_fits(input_file, &data, &rows, &cols, &xa, &ya, &naxis);

    // We iterate row by row and print columns separated by spaces.
    for (long r = 0; r < rows; r++) {
        for (long c = 0; c < cols; c++) {
            // Use standard format
            printf("%.10g", data[r * cols + c]);
            if (c < cols - 1) {
                printf(" ");
            }
        }
        printf("\n");
    }

    free(data);
    return 0;
}
