#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cblas.h>
#include <lapacke.h>
#include "common.h"

// Forward declaration
void perform_merge(char **mode_files, char **coeff_files, int n_files, const char *out_modes, const char *out_coeffs, int n_modes_keep, int ncpu);

void print_help(const char *progname) {
    fprintf(stderr, "Usage: %s [options] <input_list.txt> <output_modes.fits> <output_coeffs.fits>\n", progname);
    fprintf(stderr, "\n");
    fprintf(stderr, "Merges multiple PCA results (partial spatial modes) into a global PCA basis.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  <input_list.txt>     ASCII file listing pairs of 'modes.fits coeffs.fits'.\n");
    fprintf(stderr, "  <output_modes.fits>  Output global spatial modes.\n");
    fprintf(stderr, "  <output_coeffs.fits> Output global temporal coefficients.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -nmodes <n>     Number of global modes to retain (default: rank of merged basis).\n");
    fprintf(stderr, "  -ncpu <n>       Set number of CPU threads (OpenBLAS).\n");
    fprintf(stderr, "\n");
}

int main(int argc, char *argv[]) {
    int ncpu = 0;
    int n_modes_keep = 0; // 0 = all
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
            if (ncpu < 1) {
                fprintf(stderr, "Error: ncpu must be >= 1. Got %d.\n", ncpu);
                return 1;
            }
            i += 2;
            arg_offset += 2;
        } else if (strcmp(argv[i], "-nmodes") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "Error: -nmodes requires an argument.\n");
                return 1;
            }
            n_modes_keep = atoi(argv[i+1]);
            if (n_modes_keep < 1) {
                fprintf(stderr, "Error: nmodes must be >= 1. Got %d.\n", n_modes_keep);
                return 1;
            }
            i += 2;
            arg_offset += 2;
        } else {
            break; // Positional argument
        }
    }

    if (argc - arg_offset != 3) {
        print_help(argv[0]);
        return 1;
    }

    if (ncpu > 0) {
        openblas_set_num_threads(ncpu);
        printf("Using %d CPU cores.\n", ncpu);
    }

    const char *list_file = argv[1 + arg_offset];
    const char *out_modes = argv[2 + arg_offset];
    const char *out_coeffs = argv[3 + arg_offset];

    // Read list file
    FILE *fp = fopen(list_file, "r");
    if (!fp) {
        fprintf(stderr, "Error: Could not open list file %s\n", list_file);
        return 1;
    }

    // Count lines first
    int n_files = 0;
    char line[1024];
    while (fgets(line, sizeof(line), fp)) {
        // Skip empty lines or comments
        if (strlen(line) < 2 || line[0] == '#') continue;
        n_files++;
    }
    rewind(fp);

    if (n_files < 1) {
        fprintf(stderr, "Error: List file is empty or invalid.\n");
        fclose(fp);
        return 1;
    }

    char **mode_files = (char **)malloc(n_files * sizeof(char *));
    char **coeff_files = (char **)malloc(n_files * sizeof(char *));

    int idx = 0;
    while (fgets(line, sizeof(line), fp)) {
        if (strlen(line) < 2 || line[0] == '#') continue;

        char m_file[512];
        char c_file[512];
        if (sscanf(line, "%s %s", m_file, c_file) != 2) {
            fprintf(stderr, "Error: Invalid line format in list file: %s\n", line);
            // Cleanup
            return 1;
        }
        mode_files[idx] = strdup(m_file);
        coeff_files[idx] = strdup(c_file);
        idx++;
        if (idx >= n_files) break;
    }
    fclose(fp);

    printf("Merging PCA from %d file pairs...\n", n_files);

    perform_merge(mode_files, coeff_files, n_files, out_modes, out_coeffs, n_modes_keep, ncpu);

    // Cleanup filenames
    for(int j=0; j<n_files; j++) {
        free(mode_files[j]);
        free(coeff_files[j]);
    }
    free(mode_files);
    free(coeff_files);

    return 0;
}

// Helper to free matrix
void free_matrix(double *mat) {
    if (mat) free(mat);
}

void perform_merge(char **mode_files, char **coeff_files, int n_files, const char *out_modes, const char *out_coeffs, int n_modes_keep, int ncpu) {
    long P = 0; // Spatial dimension
    int xa = 0, ya = 0; // Spatial dims for FITS

    // Pass 1: Determine Dimensions and Total K
    long total_K = 0;
    long *Ki = (long *)malloc(n_files * sizeof(long));
    long *Ni = (long *)malloc(n_files * sizeof(long));
    long total_N = 0;

    for (int i=0; i<n_files; i++) {
        fitsfile *fptr;
        int status = 0;
        int naxis;
        long naxes[3] = {0,0,0};

        // Check Modes Dimensions
        if (fits_open_file(&fptr, mode_files[i], READONLY, &status)) {
            fits_report_error(stderr, status);
            exit(1);
        }
        fits_get_img_dim(fptr, &naxis, &status);
        fits_get_img_size(fptr, naxis, naxes, &status);
        fits_close_file(fptr, &status);

        if (naxis != 3) {
            fprintf(stderr, "Error: Mode file %s is not 3D.\n", mode_files[i]);
            exit(1);
        }
        long curr_P = naxes[0] * naxes[1];
        long curr_K = naxes[2];

        if (i == 0) {
            P = curr_P;
            xa = naxes[0];
            ya = naxes[1];
        } else {
            if (curr_P != P) {
                fprintf(stderr, "Error: Mismatched spatial dimensions in %s.\n", mode_files[i]);
                exit(1);
            }
        }
        Ki[i] = curr_K;
        total_K += curr_K;

        // Check Coeffs Dimensions
        status = 0;
        if (fits_open_file(&fptr, coeff_files[i], READONLY, &status)) {
            fits_report_error(stderr, status);
            exit(1);
        }
        fits_get_img_dim(fptr, &naxis, &status);
        fits_get_img_size(fptr, naxis, naxes, &status);
        fits_close_file(fptr, &status);

        if (naxis != 2) {
            fprintf(stderr, "Error: Coeff file %s is not 2D.\n", coeff_files[i]);
            exit(1);
        }
        // FITS 2D: width=naxes[0], height=naxes[1]
        // Coeffs logic in pca.c writes width=npca (K), height=N.
        // So width should match K.
        if (naxes[0] != curr_K) {
             fprintf(stderr, "Error: Coeff file %s width (%ld) does not match Modes depth (%ld).\n", coeff_files[i], naxes[0], curr_K);
             exit(1);
        }
        Ni[i] = naxes[1];
        total_N += Ni[i];
    }

    printf("Total Modes Input: %ld\n", total_K);
    printf("Total Samples Input: %ld\n", total_N);
    printf("Spatial Dimension: %ld (%dx%d)\n", P, xa, ya);

    // Allocate 'Weighted Stacked Modes' Matrix L (Total_K x P)
    // Beware of memory usage. total_K * P * 8 bytes.
    // 1000 * 100000 = 100MB. 10000 * 1000000 = 80GB.
    double *L = (double *)malloc(total_K * P * sizeof(double));
    if (!L) {
        fprintf(stderr, "Error: Could not allocate memory for merged modes (%ld x %ld).\n", total_K, P);
        exit(1);
    }

    long current_row_offset = 0;

    // Pass 2: Load and Scale Modes
    for (int i=0; i<n_files; i++) {
        double *Mi = NULL;
        double *Ci = NULL;
        long N_curr = Ni[i];
        long K_curr = Ki[i];

        long dummyN, dummyP;
        int dxa, dya, dnaxis;

        // Read Modes (K x P)
        read_fits(mode_files[i], &Mi, &dummyN, &dummyP, &dxa, &dya, &dnaxis);
        // Mi is flat K*P. Row-major?
        // write_fits_3d writes: xa, ya, z_dim. FITS order.
        // C order: z_dim * (xa*ya).
        // So Mi is K rows of P pixels. Correct.

        // Read Coeffs (N x K)
        read_fits(coeff_files[i], &Ci, &dummyN, &dummyP, &dxa, &dya, &dnaxis);
        // read_fits for 2D: width=P, height=N.
        // Here P=K, N=N.
        // Ci is N rows of K columns.

        // SVD of Ci: Ci = Uc * Sc * Vc^T
        // We need Vc^T and Sc.
        // Dimension of Ci is N x K. Typically N >> K.
        // K is dim(Sc) = K (or min(N,K)). Assuming N >= K.
        long min_NK = (N_curr < K_curr) ? N_curr : K_curr;

        double *Sc = (double *)malloc(min_NK * sizeof(double));
        double *Uc = (double *)malloc(N_curr * min_NK * sizeof(double)); // Not needed really
        double *Vt_c = (double *)malloc(min_NK * K_curr * sizeof(double));

        // LAPACKE_dgesdd
        // Ci is row-major N x K.
        // We want SVD.
        // Using dgesdd 'S' option to get min(N,K) singular vectors.
        lapack_int info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'S', N_curr, K_curr, Ci, K_curr,
                                         Sc, Uc, min_NK, Vt_c, K_curr);
        if (info != 0) {
            fprintf(stderr, "Error: SVD failed on coeffs file %s\n", coeff_files[i]);
            exit(1);
        }

        // Compute Li = Sc * (Vt_c * Mi)
        // Vt_c is (min_NK x K). Mi is (K x P). Result Temp is (min_NK x P).
        // Then scale by Sc (vector).

        double *Li_part = (double *)malloc(min_NK * P * sizeof(double));

        // Matrix Multiply: Vt_c * Mi
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    min_NK, P, K_curr,
                    1.0, Vt_c, K_curr,
                    Mi, P,
                    0.0, Li_part, P);

        // Scale rows by Sc
        for(long r=0; r<min_NK; r++) {
            double s_val = Sc[r];
            for(long c=0; c<P; c++) {
                Li_part[r*P + c] *= s_val;
            }
        }

        // Copy Li_part to Global L
        // L is (total_K x P). We assume total_K includes all K_curr.
        // Wait, did we sum K_curr or min_NK?
        // We allocated L based on total_K = sum(Ki).
        // But we only generated min_NK rows.
        // If N < K, we have fewer rows.
        // This is fine, we just advance offset by min_NK.
        // But we allocated more. Re-allocation or just track used rows?
        // Let's copy min_NK rows.

        memcpy(&L[current_row_offset * P], Li_part, min_NK * P * sizeof(double));
        current_row_offset += min_NK;

        free(Mi);
        free(Ci);
        free(Sc);
        free(Uc);
        free(Vt_c);
        free(Li_part);
    }

    // Update total_K to actual used rows (if N < K anywhere)
    long used_rows_L = current_row_offset;
    printf("Constructed Weighted Mode Matrix L: %ld x %ld\n", used_rows_L, P);

    // Pass 3: SVD of L -> Global Modes
    // L = U_L * S_L * V_L^T
    // We want V_L^T (or top rows of it). These are the global spatial modes.
    // L is (used_rows_L x P).
    // If P is very large, computing full V is expensive.
    // If used_rows_L << P, we only need first used_rows_L V vectors.
    // dgesdd 'S' returns min(M,N) singular vectors.
    // min(used_rows_L, P). Usually used_rows_L is smaller.

    long min_dim = (used_rows_L < P) ? used_rows_L : P;

    double *S_glob = (double *)malloc(min_dim * sizeof(double));
    double *U_glob = (double *)malloc(used_rows_L * min_dim * sizeof(double));
    double *Vt_glob = (double *)malloc(min_dim * P * sizeof(double));

    printf("Computing Global SVD...\n");
    lapack_int info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'S', used_rows_L, P, L, P,
                                     S_glob, U_glob, min_dim, Vt_glob, P);
    if (info != 0) {
        fprintf(stderr, "Error: Global SVD failed.\n");
        exit(1);
    }

    free(L); // Done with L
    free(U_glob); // Don't need U of L
    // S_glob are global singular values (sqrt of variance).
    // We could output eigenvalues here too?
    // pca.eigenvalues.txt from merged result?
    // Let's add it for consistency.
    FILE *feig = fopen("pca.eigenvalues.txt", "w");
    if (feig) {
        for(long k=0; k<min_dim; k++) {
            fprintf(feig, "%.15g\n", S_glob[k]*S_glob[k]);
        }
        fclose(feig);
        printf("Wrote pca.eigenvalues.txt (global)\n");
    }

    // Select number of output modes
    int K_out = n_modes_keep;
    if (K_out <= 0) K_out = (int)min_dim; // Default to rank of L
    if (K_out > min_dim) K_out = (int)min_dim;

    printf("Keeping %d global modes.\n", K_out);

    // M_global is top K_out rows of Vt_glob
    // Vt_glob is (min_dim x P).
    // We can just use Vt_glob pointer as M_global source (first K_out*P elements).

    printf("Writing Global Modes...\n");
    write_fits_3d(out_modes, Vt_glob, xa, ya, K_out);

    // Pass 4: Compute Global Coefficients
    // We need to re-read files to compute projection.
    // C_i_global = C_i * (M_i * M_global^T).
    // M_global^T is (P x K_out).
    // Let T_glob = M_global^T.
    // Actually we have M_global (rows). Vt_glob.
    // M_i * M_global^T: (K_i x P) * (P x K_out) -> (K_i x K_out).

    // We can't keep all C_i_global in memory easily if N is massive.
    // But we need to concatenate them.
    // allocate C_glob (total_N x K_out).

    double *C_glob = (double *)malloc(total_N * K_out * sizeof(double));
    if (!C_glob) {
        fprintf(stderr, "Error: Memory allocation for global coefficients failed.\n");
        exit(1);
    }

    long current_sample_offset = 0;

    // We need Vt_glob (Global Modes) for projection.
    // Let's keep it in memory. It is (min_dim x P).
    // We only need top K_out rows.

    for (int i=0; i<n_files; i++) {
        double *Mi = NULL;
        double *Ci = NULL;
        long N_curr = Ni[i];
        long K_curr = Ki[i];
        long dummyN, dummyP;
        int dxa, dya, dnaxis;

        read_fits(mode_files[i], &Mi, &dummyN, &dummyP, &dxa, &dya, &dnaxis);
        read_fits(coeff_files[i], &Ci, &dummyN, &dummyP, &dxa, &dya, &dnaxis);

        // Projection Matrix Pi = Mi * M_global^T
        // Mi: K_curr x P.
        // M_global: K_out x P (rows of Vt_glob).
        // We want Mi * M_global^T.
        // Using cblas_dgemm with TransB.

        double *Pi = (double *)malloc(K_curr * K_out * sizeof(double));

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                    K_curr, K_out, P,
                    1.0, Mi, P,
                    Vt_glob, P, // Vt_glob effectively M_global
                    0.0, Pi, K_out);

        // Ci_global = Ci * Pi
        // Ci: N_curr x K_curr
        // Pi: K_curr x K_out
        // Result: N_curr x K_out

        // Pointer to slice of C_glob
        double *C_dest = &C_glob[current_sample_offset * K_out];

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    N_curr, K_out, K_curr,
                    1.0, Ci, K_curr,
                    Pi, K_out,
                    0.0, C_dest, K_out);

        current_sample_offset += N_curr;

        free(Mi);
        free(Ci);
        free(Pi);
    }

    printf("Writing Global Coeffs...\n");
    write_fits_2d(out_coeffs, C_glob, K_out, total_N);

    free(C_glob);
    free(Vt_glob);
    free(S_glob);
    free(Ki);
    free(Ni);
}
