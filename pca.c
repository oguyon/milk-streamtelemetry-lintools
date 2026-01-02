#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cblas.h>
#include <lapacke.h>
#include "common.h"

// Comparator for sorting variance (descending)
int compare_double_desc(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    return (db > da) - (da > db);
}

// Comparator for sorting variance (descending) - float
int compare_float_desc(const void *a, const void *b) {
    float da = *(const float *)a;
    float db = *(const float *)b;
    return (db > da) - (da > db);
}

int main(int argc, char *argv[]) {
    int ncpu = 0;
    int use_float = 0;
    int arg_offset = 0;
    char *mask_file = NULL;
    char *automask_arg = NULL;

    // Argument parsing
    int i = 1;
    while (i < argc) {
        if (strcmp(argv[i], "-ncpu") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: -ncpu requires an argument.\n");
                return 1;
            }
            ncpu = atoi(argv[i+1]);
            if (ncpu < 1) {
                fprintf(stderr, "Error: ncpu must be >= 1. Got %d.\n", ncpu);
                return 1;
            }
            i += 2;
            arg_offset += 2;
        } else if (strcmp(argv[i], "-float") == 0) {
            use_float = 1;
            i += 1;
            arg_offset += 1;
        } else if (strcmp(argv[i], "-mask") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: -mask requires an argument.\n");
                return 1;
            }
            mask_file = argv[i+1];
            i += 2;
            arg_offset += 2;
        } else if (strcmp(argv[i], "-automask") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: -automask requires an argument.\n");
                return 1;
            }
            automask_arg = argv[i+1];
            i += 2;
            arg_offset += 2;
        } else {
            break; // Positional argument found
        }
    }

    if (argc - arg_offset != 5) {
        fprintf(stderr, "Usage: %s [-ncpu <n>] [-float] [-mask <file>] [-automask <arg>] <npca> <input.fits> <modes.fits> <coeffs.fits>\n", argv[0]);
        return 1;
    }

    if (ncpu > 0) {
        openblas_set_num_threads(ncpu);
        printf("Using %d CPU cores.\n", ncpu);
    }
    if (use_float) {
        printf("Using Single Precision (float).\n");
    } else {
        printf("Using Double Precision (double).\n");
    }

    int npca = atoi(argv[1 + arg_offset]);
    const char *infile = argv[2 + arg_offset];
    const char *modes_file = argv[3 + arg_offset];
    const char *coeffs_file = argv[4 + arg_offset];

    if (npca < 1) {
        fprintf(stderr, "Error: npca must be >= 1. Got %d.\n", npca);
        return 1;
    }

    printf("PCA Configuration:\n");
    printf("  npca: %d\n", npca);
    printf("  Input: %s\n", infile);
    printf("  Modes Output: %s\n", modes_file);
    printf("  Coeffs Output: %s\n", coeffs_file);
    if (mask_file) printf("  Mask Input: %s\n", mask_file);
    if (automask_arg) printf("  Auto Mask: %s\n", automask_arg);

    if (mask_file && automask_arg) {
        fprintf(stderr, "Error: Cannot specify both -mask and -automask.\n");
        return 1;
    }

    if (use_float) {
        float *X = NULL;
        long N, P;
        int xa, ya;
        int naxis;

        printf("Reading %s...\n", infile);
        read_fits_float(infile, &X, &N, &P, &xa, &ya, &naxis);
        printf("  Dimensions: %ld samples x %ld pixels (%d x %d)\n", N, P, xa, ya);

        // Mask generation / loading
        unsigned char *mask = NULL;
        if (automask_arg) {
            printf("Computing Auto Mask...\n");
            float *variances = (float *)malloc(P * sizeof(float));
            float *sorted_vars = (float *)malloc(P * sizeof(float));

            // Compute variance for each pixel
            for(long j=0; j<P; j++) {
                double mean = 0.0;
                double var = 0.0;
                for(long k=0; k<N; k++) {
                    mean += X[k*P + j];
                }
                mean /= N;
                for(long k=0; k<N; k++) {
                    double diff = X[k*P + j] - mean;
                    var += diff * diff;
                }
                variances[j] = (float)(var / N); // Or N-1
                sorted_vars[j] = variances[j];
            }

            qsort(sorted_vars, P, sizeof(float), compare_float_desc);

            long n_active = 0;
            if (automask_arg[0] == 'n') {
                n_active = atoi(automask_arg + 1);
            } else if (automask_arg[0] == 'f') {
                float frac = atof(automask_arg + 1);
                n_active = (long)(frac * P);
            } else {
                 fprintf(stderr, "Error: Invalid automask argument format. Use n%%d or f%%f.\n");
                 return 1;
            }

            if (n_active > P) n_active = P;
            if (n_active < 1) n_active = 1;

            float threshold = sorted_vars[n_active - 1];
            printf("  Auto Mask Threshold: %g (keeping %ld pixels)\n", threshold, n_active);

            mask = (unsigned char *)calloc(P, sizeof(unsigned char));
            float *mask_float = (float *)malloc(P * sizeof(float)); // For writing to FITS

            // Re-scan to handle ties correctly? Simple thresholding works if unique or we don't care about exact count ties.
            // Strict check: if multiple pixels have exactly threshold value, we might overshoot.
            // Better to iterate and select top n_active indices, but threshold is simpler.
            // Let's use threshold.

            long count = 0;
            for(long j=0; j<P; j++) {
                if (variances[j] >= threshold && count < n_active) {
                    mask[j] = 1;
                    mask_float[j] = 1.0f;
                    count++;
                } else {
                    mask[j] = 0;
                    mask_float[j] = 0.0f;
                }
            }
            // If ties prevented reaching n_active, we might add more?
            // If ties caused us to pick too few (e.g. threshold is 0 because many are 0), we might need logic.
            // But usually variance is continuous enough.

            printf("Writing pca.automask.fits...\n");
            write_fits_2d_float("pca.automask.fits", mask_float, xa, ya); // FITS is 2D image

            free(variances);
            free(sorted_vars);
            free(mask_float);

        } else if (mask_file) {
            printf("Loading mask from %s...\n", mask_file);
            float *mask_data = NULL;
            long M_h, M_w;
            int m_xa, m_ya;
            int m_naxis;
            read_fits_float(mask_file, &mask_data, &M_h, &M_w, &m_xa, &m_ya, &m_naxis);

            if (M_w * M_h != P) { // Dimensions check
                 fprintf(stderr, "Error: Mask dimensions (%ldx%ld) do not match input spatial dimensions (%dx%d=%ld).\n",
                         (m_naxis==3?M_w:m_xa), (m_naxis==3?M_h:1), xa, ya, P);
                 return 1;
            }

            mask = (unsigned char *)malloc(P * sizeof(unsigned char));
            for(long j=0; j<P; j++) {
                mask[j] = (mask_data[j] > 0.5f) ? 1 : 0;
            }
            free(mask_data);
        }

        float *X_reduced = NULL;
        long P_reduced = P;
        if (mask) {
            long count = 0;
            for(long j=0; j<P; j++) if(mask[j]) count++;
            P_reduced = count;
            printf("  Mask active pixels: %ld / %ld\n", P_reduced, P);

            if (P_reduced == 0) {
                fprintf(stderr, "Error: Mask has 0 active pixels.\n");
                return 1;
            }

            X_reduced = (float *)malloc(N * P_reduced * sizeof(float));
            for(long i=0; i<N; i++) {
                long col_idx = 0;
                for(long j=0; j<P; j++) {
                    if (mask[j]) {
                        X_reduced[i * P_reduced + col_idx] = X[i * P + j];
                        col_idx++;
                    }
                }
            }
            free(X); // Free original large X
            X = X_reduced;
        } else {
             // No mask, P_reduced = P. X stays as is.
        }

        if (npca > N) {
            fprintf(stderr, "Warning: npca (%d) > N (%ld). Reducing npca to N.\n", npca, N);
            npca = (int)N;
        }
        if (npca > P_reduced) {
            fprintf(stderr, "Warning: npca (%d) > P_active (%ld). Reducing npca to P_active.\n", npca, P_reduced);
            npca = (int)P_reduced;
        }

        // SVD
        long K = (N < P_reduced) ? N : P_reduced;
        float *S = (float *)malloc(K * sizeof(float));
        float *U = (float *)malloc(N * K * sizeof(float));
        float *Vt = (float *)malloc(K * P_reduced * sizeof(float));

        printf("Computing SVD (using sgesdd)...\n");
        lapack_int info = LAPACKE_sgesdd(LAPACK_ROW_MAJOR, 'S', N, P_reduced, X, P_reduced,
                                         S, U, K, Vt, P_reduced);
        if (info != 0) {
            fprintf(stderr, "LAPACKE_sgesdd failed with error code %d\n", info);
            return 1;
        }

        // Write Eigenvalues
        FILE *feig = fopen("pca.eigenvalues.txt", "w");
        if (feig) {
            for(long k=0; k<K; k++) {
                fprintf(feig, "%.6g\n", S[k]*S[k]);
            }
            fclose(feig);
            printf("Wrote pca.eigenvalues.txt\n");
        } else {
            fprintf(stderr, "Warning: Could not write pca.eigenvalues.txt\n");
        }

        // Coeffs = U * S (first npca columns)
        float *Coeffs = (float *)malloc(N * npca * sizeof(float));
        for(long i=0; i<N; i++) {
            for(int j=0; j<npca; j++) {
                Coeffs[i * npca + j] = U[i * K + j] * S[j];
            }
        }

        // Check sign of the first mode
        double coeff_mean = 0.0;
        for(long i=0; i<N; i++) {
            coeff_mean += Coeffs[i * npca + 0];
        }
        coeff_mean /= N;

        if (coeff_mean < 0) {
            printf("Flipping sign of Mode 0 to ensure positive average coefficient.\n");
            // Flip Coeffs column 0
            for(long i=0; i<N; i++) {
                Coeffs[i * npca + 0] *= -1.0f;
            }
            // Flip Vt row 0 (which becomes Mode 0)
            for(long j=0; j<P_reduced; j++) {
                Vt[0 * P_reduced + j] *= -1.0f;
            }
        }

        // Modes reconstruction
        // Modes should be size npca x P (full size). Inactive pixels are 0.
        float *Modes = (float *)calloc(npca * P, sizeof(float)); // calloc inits to 0
        for(int i=0; i<npca; i++) {
            long reduced_idx = 0;
            for(long j=0; j<P; j++) {
                if (!mask || mask[j]) { // If no mask, mask is NULL, so effectively all active
                     Modes[i * P + j] = Vt[i * P_reduced + reduced_idx];
                     reduced_idx++;
                } else {
                     Modes[i * P + j] = 0.0f;
                }
            }
        }

        printf("Writing Modes to %s...\n", modes_file);
        write_fits_3d_float(modes_file, Modes, xa, ya, npca);

        printf("Writing Coeffs to %s...\n", coeffs_file);
        write_fits_2d_float(coeffs_file, Coeffs, npca, N);

        printf("Done.\n");

        free(X); free(S); free(U); free(Vt);
        free(Coeffs); free(Modes);
        if (mask) free(mask);

    } else {
        // Double precision (original code + mask logic)
        double *X = NULL;
        long N, P;
        int xa, ya;
        int naxis;

        printf("Reading %s...\n", infile);
        read_fits(infile, &X, &N, &P, &xa, &ya, &naxis);
        printf("  Dimensions: %ld samples x %ld pixels (%d x %d)\n", N, P, xa, ya);

        // Mask generation / loading (Double precision path)
        unsigned char *mask = NULL;
        if (automask_arg) {
            printf("Computing Auto Mask...\n");
            double *variances = (double *)malloc(P * sizeof(double));
            double *sorted_vars = (double *)malloc(P * sizeof(double));

            for(long j=0; j<P; j++) {
                double mean = 0.0;
                double var = 0.0;
                for(long k=0; k<N; k++) {
                    mean += X[k*P + j];
                }
                mean /= N;
                for(long k=0; k<N; k++) {
                    double diff = X[k*P + j] - mean;
                    var += diff * diff;
                }
                variances[j] = var / N;
                sorted_vars[j] = variances[j];
            }

            qsort(sorted_vars, P, sizeof(double), compare_double_desc);

            long n_active = 0;
            if (automask_arg[0] == 'n') {
                n_active = atoi(automask_arg + 1);
            } else if (automask_arg[0] == 'f') {
                float frac = atof(automask_arg + 1);
                n_active = (long)(frac * P);
            } else {
                 fprintf(stderr, "Error: Invalid automask argument format. Use n%%d or f%%f.\n");
                 return 1;
            }

            if (n_active > P) n_active = P;
            if (n_active < 1) n_active = 1;

            double threshold = sorted_vars[n_active - 1];
            printf("  Auto Mask Threshold: %g (keeping %ld pixels)\n", threshold, n_active);

            mask = (unsigned char *)calloc(P, sizeof(unsigned char));
            double *mask_double = (double *)malloc(P * sizeof(double));

            long count = 0;
            for(long j=0; j<P; j++) {
                if (variances[j] >= threshold && count < n_active) {
                    mask[j] = 1;
                    mask_double[j] = 1.0;
                    count++;
                } else {
                    mask[j] = 0;
                    mask_double[j] = 0.0;
                }
            }

            printf("Writing pca.automask.fits...\n");
            // FITS images are usually width=xa, height=ya.
            // write_fits_2d uses width, height.
            // P = xa*ya.
            write_fits_2d("pca.automask.fits", mask_double, xa, ya);

            free(variances);
            free(sorted_vars);
            free(mask_double);

        } else if (mask_file) {
            printf("Loading mask from %s...\n", mask_file);
            double *mask_data = NULL;
            long M_h, M_w;
            int m_xa, m_ya;
            int m_naxis;
            read_fits(mask_file, &mask_data, &M_h, &M_w, &m_xa, &m_ya, &m_naxis);

            if (M_w * M_h != P) {
                 fprintf(stderr, "Error: Mask dimensions do not match input spatial dimensions.\n");
                 return 1;
            }

            mask = (unsigned char *)malloc(P * sizeof(unsigned char));
            for(long j=0; j<P; j++) {
                mask[j] = (mask_data[j] > 0.5) ? 1 : 0;
            }
            free(mask_data);
        }

        double *X_reduced = NULL;
        long P_reduced = P;
        if (mask) {
            long count = 0;
            for(long j=0; j<P; j++) if(mask[j]) count++;
            P_reduced = count;
            printf("  Mask active pixels: %ld / %ld\n", P_reduced, P);

            if (P_reduced == 0) {
                fprintf(stderr, "Error: Mask has 0 active pixels.\n");
                return 1;
            }

            X_reduced = (double *)malloc(N * P_reduced * sizeof(double));
            for(long i=0; i<N; i++) {
                long col_idx = 0;
                for(long j=0; j<P; j++) {
                    if (mask[j]) {
                        X_reduced[i * P_reduced + col_idx] = X[i * P + j];
                        col_idx++;
                    }
                }
            }
            free(X);
            X = X_reduced;
        }

        if (npca > N) {
            fprintf(stderr, "Warning: npca (%d) > N (%ld). Reducing npca to N.\n", npca, N);
            npca = (int)N;
        }
        if (npca > P_reduced) {
            fprintf(stderr, "Warning: npca (%d) > P (%ld). Reducing npca to P.\n", npca, P_reduced);
            npca = (int)P_reduced;
        }

        // SVD
        long K = (N < P_reduced) ? N : P_reduced;
        double *S = (double *)malloc(K * sizeof(double));
        double *U = (double *)malloc(N * K * sizeof(double));
        double *Vt = (double *)malloc(K * P_reduced * sizeof(double));

        printf("Computing SVD (using dgesdd)...\n");
        // Use dgesdd for speed (Divide and Conquer)
        lapack_int info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'S', N, P_reduced, X, P_reduced,
                                         S, U, K, Vt, P_reduced);
        if (info != 0) {
            fprintf(stderr, "LAPACKE_dgesdd failed with error code %d\n", info);
            return 1;
        }

        // Write Eigenvalues
        FILE *feig = fopen("pca.eigenvalues.txt", "w");
        if (feig) {
            for(long k=0; k<K; k++) {
                fprintf(feig, "%.15g\n", S[k]*S[k]);
            }
            fclose(feig);
            printf("Wrote pca.eigenvalues.txt\n");
        } else {
            fprintf(stderr, "Warning: Could not write pca.eigenvalues.txt\n");
        }

        // Coeffs = U * S (first npca columns)
        double *Coeffs = (double *)malloc(N * npca * sizeof(double));
        for(long i=0; i<N; i++) {
            for(int j=0; j<npca; j++) {
                Coeffs[i * npca + j] = U[i * K + j] * S[j];
            }
        }

        // Check sign of the first mode
        double coeff_mean = 0.0;
        for(long i=0; i<N; i++) {
            coeff_mean += Coeffs[i * npca + 0];
        }
        coeff_mean /= N;

        if (coeff_mean < 0) {
            printf("Flipping sign of Mode 0 to ensure positive average coefficient.\n");
            // Flip Coeffs column 0
            for(long i=0; i<N; i++) {
                Coeffs[i * npca + 0] *= -1.0;
            }
            // Flip Vt row 0 (which becomes Mode 0)
            for(long j=0; j<P_reduced; j++) {
                Vt[0 * P_reduced + j] *= -1.0;
            }
        }

        // Modes = Vt (first npca rows) expanded
        double *Modes = (double *)calloc(npca * P, sizeof(double));
        for(int i=0; i<npca; i++) {
             long reduced_idx = 0;
            for(long j=0; j<P; j++) {
                if (!mask || mask[j]) {
                    Modes[i * P + j] = Vt[i * P_reduced + reduced_idx];
                    reduced_idx++;
                } else {
                    Modes[i * P + j] = 0.0;
                }
            }
        }

        printf("Writing Modes to %s...\n", modes_file);
        write_fits_3d(modes_file, Modes, xa, ya, npca);

        printf("Writing Coeffs to %s...\n", coeffs_file);
        write_fits_2d(coeffs_file, Coeffs, npca, N);

        printf("Done.\n");

        free(X); free(S); free(U); free(Vt);
        free(Coeffs); free(Modes);
        if (mask) free(mask);
    }

    return 0;
}
