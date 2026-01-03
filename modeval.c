#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cblas.h>
#include "common.h"

void print_help(const char *progname) {
    fprintf(stderr, "Usage: %s [options] <frames.fits> <modes.fits> <output.fits>\n", progname);
    fprintf(stderr, "\n");
    fprintf(stderr, "Performs modal decomposition of a series of frames onto a set of modes.\n");
    fprintf(stderr, "Assumes modes are orthonormal.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  <frames.fits>   Input frames (3D FITS cube: x, y, N_frames).\n");
    fprintf(stderr, "  <modes.fits>    Input modes (3D FITS cube: x, y, N_modes).\n");
    fprintf(stderr, "  <output.fits>   Output coefficients (2D FITS: N_frames x N_modes).\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -ncpu <n>       Set number of CPU threads (OpenBLAS).\n");
    fprintf(stderr, "  -float          Use single precision (float).\n");
    fprintf(stderr, "  -xp <file>      Compute and save mode cross-products (Modes * Modes^T) to check orthogonality.\n");
    fprintf(stderr, "  -ascii <file>   Save coefficients to a text file in addition to FITS.\n");
    fprintf(stderr, "\n");
}

int main(int argc, char *argv[]) {
    int ncpu = 0;
    int use_float = 0;
    char *xp_file = NULL;
    char *ascii_file = NULL;
    int arg_offset = 0;

    // Argument parsing
    int i = 1;
    while (i < argc) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_help(argv[0]);
            return 0;
        } else if (strcmp(argv[i], "-ncpu") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: -ncpu requires an argument.\n");
                return 1;
            }
            ncpu = atoi(argv[i+1]);
            i += 2;
            arg_offset += 2;
        } else if (strcmp(argv[i], "-float") == 0) {
            use_float = 1;
            i += 1;
            arg_offset += 1;
        } else if (strcmp(argv[i], "-xp") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: -xp requires an argument.\n");
                return 1;
            }
            xp_file = argv[i+1];
            i += 2;
            arg_offset += 2;
        } else if (strcmp(argv[i], "-ascii") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: -ascii requires an argument.\n");
                return 1;
            }
            ascii_file = argv[i+1];
            i += 2;
            arg_offset += 2;
        } else {
            break; // Positional argument
        }
    }

    if (argc - arg_offset != 3) {
        if (argc - arg_offset < 3) {
            fprintf(stderr, "Error: Missing required arguments. Expected 3 positional arguments, found %d.\n", argc - arg_offset - 1);
        } else {
            fprintf(stderr, "Error: Too many arguments.\n");
        }
        print_help(argv[0]);
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

    const char *frames_file = argv[1 + arg_offset];
    const char *modes_file = argv[2 + arg_offset];
    const char *out_file = argv[3 + arg_offset];

    if (use_float) {
        float *Frames = NULL;
        long N, P_frames;
        int xa, ya, naxis;

        printf("Reading Frames %s...\n", frames_file);
        read_fits_float(frames_file, &Frames, &N, &P_frames, &xa, &ya, &naxis);
        printf("  Frames: %ld samples x %ld pixels (%d x %d)\n", N, P_frames, xa, ya);

        float *Modes = NULL;
        long K, P_modes;
        int mxa, mya, mnaxis;

        printf("Reading Modes %s...\n", modes_file);
        read_fits_float(modes_file, &Modes, &K, &P_modes, &mxa, &mya, &mnaxis);
        printf("  Modes: %ld modes x %ld pixels (%d x %d)\n", K, P_modes, mxa, mya);

        if (P_frames != P_modes) {
            fprintf(stderr, "Error: Pixel dimension mismatch (Frames: %ld, Modes: %ld).\n", P_frames, P_modes);
            free(Frames); free(Modes);
            return 1;
        }
        long P = P_frames;

        if (xp_file) {
            printf("Computing Mode Cross-Products (-xp)...\n");
            // XP = M * M^T. M is K x P. XP is K x K.
            float *XP = (float *)malloc(K * K * sizeof(float));
            // syrk: C = alpha*A*A' + beta*C.
            // 'U'pper or 'L'ower. dgemm is safer for full matrix output if not using syrk-specific pack.
            // But we can use dgemm: M (KxP) * M^T (PxK) -> KxK.
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                        K, K, P,
                        1.0f, Modes, P,
                        Modes, P,
                        0.0f, XP, K);

            printf("Writing XP to %s...\n", xp_file);
            write_fits_2d_float(xp_file, XP, K, K);
            free(XP);
        }

        // Coefficients = Frames * Modes^T
        // Frames: N x P
        // Modes: K x P
        // Coeffs: N x K

        printf("Computing Coefficients...\n");
        float *Coeffs = (float *)malloc(N * K * sizeof(float));
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                    N, K, P,
                    1.0f, Frames, P,
                    Modes, P,
                    0.0f, Coeffs, K);

        printf("Writing Coefficients to %s...\n", out_file);
        write_fits_2d_float(out_file, Coeffs, K, N); // Width=K, Height=N

        if (ascii_file) {
            printf("Writing ASCII output to %s...\n", ascii_file);
            FILE *fp = fopen(ascii_file, "w");
            if (fp) {
                for(long i=0; i<N; i++) {
                    for(long j=0; j<K; j++) {
                        fprintf(fp, "%14.7e ", Coeffs[i*K + j]);
                    }
                    fprintf(fp, "\n");
                }
                fclose(fp);
            } else {
                fprintf(stderr, "Error: Could not write to %s\n", ascii_file);
            }
        }

        free(Frames);
        free(Modes);
        free(Coeffs);

    } else {
        // Double precision
        double *Frames = NULL;
        long N, P_frames;
        int xa, ya, naxis;

        printf("Reading Frames %s...\n", frames_file);
        read_fits(frames_file, &Frames, &N, &P_frames, &xa, &ya, &naxis);
        printf("  Frames: %ld samples x %ld pixels (%d x %d)\n", N, P_frames, xa, ya);

        double *Modes = NULL;
        long K, P_modes;
        int mxa, mya, mnaxis;

        printf("Reading Modes %s...\n", modes_file);
        read_fits(modes_file, &Modes, &K, &P_modes, &mxa, &mya, &mnaxis);
        printf("  Modes: %ld modes x %ld pixels (%d x %d)\n", K, P_modes, mxa, mya);

        if (P_frames != P_modes) {
            fprintf(stderr, "Error: Pixel dimension mismatch (Frames: %ld, Modes: %ld).\n", P_frames, P_modes);
            free(Frames); free(Modes);
            return 1;
        }
        long P = P_frames;

        if (xp_file) {
            printf("Computing Mode Cross-Products (-xp)...\n");
            // XP = M * M^T. M is K x P. XP is K x K.
            double *XP = (double *)malloc(K * K * sizeof(double));
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                        K, K, P,
                        1.0, Modes, P,
                        Modes, P,
                        0.0, XP, K);

            printf("Writing XP to %s...\n", xp_file);
            write_fits_2d(xp_file, XP, K, K);
            free(XP);
        }

        // Coefficients = Frames * Modes^T
        // Frames: N x P
        // Modes: K x P
        // Coeffs: N x K

        printf("Computing Coefficients...\n");
        double *Coeffs = (double *)malloc(N * K * sizeof(double));
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                    N, K, P,
                    1.0, Frames, P,
                    Modes, P,
                    0.0, Coeffs, K);

        printf("Writing Coefficients to %s...\n", out_file);
        write_fits_2d(out_file, Coeffs, K, N); // Width=K, Height=N

        if (ascii_file) {
            printf("Writing ASCII output to %s...\n", ascii_file);
            FILE *fp = fopen(ascii_file, "w");
            if (fp) {
                for(long i=0; i<N; i++) {
                    for(long j=0; j<K; j++) {
                        fprintf(fp, "%18.10e ", Coeffs[i*K + j]);
                    }
                    fprintf(fp, "\n");
                }
                fclose(fp);
            } else {
                fprintf(stderr, "Error: Could not write to %s\n", ascii_file);
            }
        }

        free(Frames);
        free(Modes);
        free(Coeffs);
    }

    printf("Done.\n");
    return 0;
}
