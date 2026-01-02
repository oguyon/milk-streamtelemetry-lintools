#ifndef COMMON_H
#define COMMON_H

#include <fitsio.h>

// Helper macro for error checking
#define CHECK_STATUS(status) if (status) { fits_report_error(stderr, status); exit(status); }

// Reads a FITS file (2D or 3D).
// For 3D: N=naxis3, P=naxis1*naxis2. xa=naxis1, ya=naxis2.
// For 2D: N=naxis2, P=naxis1. xa=naxis1, ya=1.
// data is allocated and filled with flattened array (row-major: N x P).
void read_fits(const char *filename, double **data, long *N, long *P, int *xa, int *ya);

// Writes a 3D FITS cube (xa, ya, z_dim).
// data is expected to be z_dim * (xa*ya) elements.
void write_fits_3d(const char *filename, double *data, int xa, int ya, int z_dim);

// Writes a 2D FITS image (width, height).
// data is expected to be height * width elements.
// Note: FITS stores (width, height). C stores (height, width).
// width = P (columns/pixels), height = N (samples).
void write_fits_2d(const char *filename, double *data, long width, long height);

// Centers the columns of the N x P data matrix.
void center_columns(double *data, long N, long P);

#endif
