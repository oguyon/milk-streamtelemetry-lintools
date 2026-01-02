#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void read_fits(const char *filename, double **data, long *N, long *P, int *xa, int *ya, int *naxis_out) {
    fitsfile *fptr;
    int status = 0;
    int naxis;
    long naxes[3] = {0, 0, 0};

    fits_open_file(&fptr, filename, READONLY, &status);
    CHECK_STATUS(status);

    fits_get_img_dim(fptr, &naxis, &status);
    CHECK_STATUS(status);

    if (naxis < 2 || naxis > 3) {
        fprintf(stderr, "Error: %s is not a 2D or 3D FITS file (NAXIS=%d).\n", filename, naxis);
        fits_close_file(fptr, &status);
        exit(1);
    }

    if (naxis_out) *naxis_out = naxis;

    fits_get_img_size(fptr, naxis, naxes, &status);
    CHECK_STATUS(status);

    // Check for integer overflow if dimensions are too large for int (LAPACK/BLAS limitation)
    for(int i=0; i<naxis; i++) {
        if (naxes[i] > 2147483647L) {
            fprintf(stderr, "Error: Dimension %d size %ld exceeds INT_MAX. LAPACK/BLAS operations may fail.\n", i, naxes[i]);
            fits_close_file(fptr, &status);
            exit(1);
        }
    }

    if (naxis == 3) {
        *xa = (int)naxes[0];
        *ya = (int)naxes[1];
        *N = naxes[2];
        *P = (*xa) * (*ya);
    } else { // naxis == 2
        *xa = (int)naxes[0]; // width (P)
        *ya = 1;
        *N = naxes[1]; // height (rows)
        *P = (*xa);
    }

    long total_elements = (*N) * (*P);
    *data = (double *)malloc(total_elements * sizeof(double));
    if (*data == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for %ld doubles.\n", total_elements);
        exit(1);
    }

    long fpixel[3] = {1, 1, 1};
    fits_read_pix(fptr, TDOUBLE, fpixel, total_elements, NULL, *data, NULL, &status);
    CHECK_STATUS(status);

    fits_close_file(fptr, &status);
    CHECK_STATUS(status);
}

void write_fits_3d(const char *filename, double *data, int xa, int ya, int z_dim) {
    if (xa < 1 || ya < 1 || z_dim < 1) {
        fprintf(stderr, "Error: Invalid dimensions for 3D FITS write: %d x %d x %d\n", xa, ya, z_dim);
        exit(1);
    }

    fitsfile *fptr;
    int status = 0;
    long naxes[3];
    naxes[0] = xa;
    naxes[1] = ya;
    naxes[2] = z_dim;

    remove(filename);

    fits_create_file(&fptr, filename, &status);
    CHECK_STATUS(status);

    fits_create_img(fptr, DOUBLE_IMG, 3, naxes, &status);
    CHECK_STATUS(status);

    long total_elements = (long)xa * ya * z_dim;
    long fpixel[3] = {1, 1, 1};
    fits_write_pix(fptr, TDOUBLE, fpixel, total_elements, data, &status);
    CHECK_STATUS(status);

    fits_close_file(fptr, &status);
    CHECK_STATUS(status);
}

void write_fits_2d(const char *filename, double *data, long width, long height) {
    if (width < 1 || height < 1) {
        fprintf(stderr, "Error: Invalid dimensions for 2D FITS write: %ld x %ld\n", width, height);
        exit(1);
    }

    fitsfile *fptr;
    int status = 0;
    long naxes[2];
    naxes[0] = width;
    naxes[1] = height;

    remove(filename);

    fits_create_file(&fptr, filename, &status);
    CHECK_STATUS(status);

    fits_create_img(fptr, DOUBLE_IMG, 2, naxes, &status);
    CHECK_STATUS(status);

    long total_elements = width * height;
    long fpixel[2] = {1, 1};
    fits_write_pix(fptr, TDOUBLE, fpixel, total_elements, data, &status);
    CHECK_STATUS(status);

    fits_close_file(fptr, &status);
    CHECK_STATUS(status);
}
