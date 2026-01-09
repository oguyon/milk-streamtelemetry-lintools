#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <png.h>

// Helper to parse comma-separated list of integers
int* parse_int_list(const char *str, int *count) {
    int capacity = 10;
    int *list = malloc(capacity * sizeof(int));
    *count = 0;

    const char *p = str;
    while (*p) {
        char *end;
        long val = strtol(p, &end, 10);
        if (p == end) break; // No number found

        if (*count >= capacity) {
            capacity *= 2;
            list = realloc(list, capacity * sizeof(int));
        }
        list[(*count)++] = (int)val;

        p = end;
        if (*p == ',') p++;
    }
    return list;
}

void write_png(const char *filename, unsigned char *image, int width, int height) {
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        fprintf(stderr, "Error: Could not open file %s for writing\n", filename);
        exit(1);
    }

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) {
        fprintf(stderr, "Error: Could not allocate write struct\n");
        exit(1);
    }

    png_infop info = png_create_info_struct(png);
    if (!info) {
        fprintf(stderr, "Error: Could not allocate info struct\n");
        exit(1);
    }

    if (setjmp(png_jmpbuf(png))) {
        fprintf(stderr, "Error during png creation\n");
        exit(1);
    }

    png_init_io(png, fp);

    // 8-bit depth, grayscale
    png_set_IHDR(png, info, width, height,
                 8, PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

    png_write_info(png, info);

    // Allocate row pointers
    png_bytep *row_pointers = (png_bytep *) malloc(sizeof(png_bytep) * height);
    for (int y = 0; y < height; y++) {
        row_pointers[y] = &image[y * width];
    }

    png_write_image(png, row_pointers);
    png_write_end(png, NULL);

    free(row_pointers);
    png_destroy_write_struct(&png, &info);
    fclose(fp);
}

void print_usage(const char *progname) {
    fprintf(stderr, "Usage: %s <input.fits> <output.png> [options]\n", progname);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -n <N>            Number of automatically selected slices (default: 5)\n");
    fprintf(stderr, "  -slices <list>    Comma-separated list of slice indices (e.g. 0,10,20)\n");
    fprintf(stderr, "  -geom <cols>x<rows> Tile geometry (default: N x 1)\n");
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        print_usage(argv[0]);
        return 1;
    }

    const char *input_file = argv[1];
    const char *output_file = argv[2];

    int n_auto = 5;
    int *explicit_indices = NULL;
    int explicit_count = 0;
    int geom_cols = 0;
    int geom_rows = 0;

    // Parse arguments
    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "-n") == 0) {
            if (i + 1 < argc) {
                n_auto = atoi(argv[i+1]);
                i++;
            }
        } else if (strcmp(argv[i], "-slices") == 0) {
            if (i + 1 < argc) {
                explicit_indices = parse_int_list(argv[i+1], &explicit_count);
                i++;
            }
        } else if (strcmp(argv[i], "-geom") == 0) {
            if (i + 1 < argc) {
                if (sscanf(argv[i+1], "%dx%d", &geom_cols, &geom_rows) != 2) {
                    fprintf(stderr, "Error: Invalid geometry format. Use COLSxROWS (e.g. 4x3)\n");
                    return 1;
                }
                i++;
            }
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            print_usage(argv[0]);
            return 1;
        }
    }

    // Read FITS
    float *data = NULL;
    long n_frames_long = 0, pixels_per_frame = 0;
    int width, height, nframes;
    int xa, ya, naxis;

    read_fits_float(input_file, &data, &n_frames_long, &pixels_per_frame, &xa, &ya, &naxis);

    if (naxis != 3) {
        fprintf(stderr, "Error: Input FITS file must be 3D cube (found NAXIS=%d).\n", naxis);
        free(data);
        if (explicit_indices) free(explicit_indices);
        return 1;
    }

    nframes = (int)n_frames_long;
    width = xa;
    height = ya;
    long pix_per_frame = (long)width * height;

    // Determine slices to extract
    int N_slices = 0;
    int *indices = NULL;

    if (explicit_indices) {
        N_slices = explicit_count;
        indices = explicit_indices;
        // Validate indices
        for (int i=0; i<N_slices; i++) {
            if (indices[i] < 0 || indices[i] >= nframes) {
                fprintf(stderr, "Error: Slice index %d out of bounds [0, %d]\n", indices[i], nframes-1);
                free(data);
                free(indices);
                return 1;
            }
        }
    } else {
        N_slices = n_auto;
        if (N_slices <= 0) N_slices = 1;
        if (N_slices > nframes) {
             fprintf(stderr, "Warning: N (%d) > nframes (%d). Setting N = nframes.\n", N_slices, nframes);
             N_slices = nframes;
        }
        indices = (int*) malloc(N_slices * sizeof(int));
        if (N_slices == 1) {
            indices[0] = 0;
        } else {
            double step = (double)(nframes - 1) / (double)(N_slices - 1);
            for (int i = 0; i < N_slices; i++) {
                indices[i] = (int)(i * step + 0.5);
                if (indices[i] >= nframes) indices[i] = nframes - 1;
            }
        }
    }

    // Determine Geometry
    if (geom_cols > 0 && geom_rows > 0) {
        // User specified geometry
        // Check if fits
        if (geom_cols * geom_rows < N_slices) {
            fprintf(stderr, "Warning: Geometry %dx%d (%d slots) is smaller than number of slices (%d).\n",
                    geom_cols, geom_rows, geom_cols*geom_rows, N_slices);
            // We'll just display the first slots
        }
    } else {
        // Default: Side by side
        geom_cols = N_slices;
        geom_rows = 1;
    }

    printf("Input: %s (%dx%dx%d)\n", input_file, width, height, nframes);
    printf("Extracting %d slices -> Output Grid %dx%d\n", N_slices, geom_cols, geom_rows);

    // Find Global Min/Max across selected slices
    float min_val = 1e30f;
    float max_val = -1e30f;
    int first = 1;

    for (int i = 0; i < N_slices; i++) {
        int t = indices[i];
        long offset = t * pix_per_frame;
        for (long p = 0; p < pix_per_frame; p++) {
            float v = data[offset + p];
            if (first || v < min_val) { min_val = v; first = 0; }
            if (v > max_val) max_val = v;
        }
    }

    printf("Global range: [%g, %g]\n", min_val, max_val);
    if (max_val == min_val) max_val = min_val + 1.0f;

    // Construct Output Image
    int out_width = geom_cols * width;
    int out_height = geom_rows * height;

    unsigned char *png_data = (unsigned char*) calloc(out_width * out_height, sizeof(unsigned char));

    for (int i = 0; i < N_slices; i++) {
        // Check if i exceeds geometry capacity
        if (i >= geom_cols * geom_rows) break;

        int grid_x = i % geom_cols;
        int grid_y = i / geom_cols;

        int t = indices[i];
        long offset = t * pix_per_frame;

        int start_x = grid_x * width;
        int start_y = grid_y * height;

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                float v = data[offset + y * width + x];
                float norm = (v - min_val) / (max_val - min_val);
                if (norm < 0) norm = 0;
                if (norm > 1) norm = 1;
                unsigned char byte_val = (unsigned char)(norm * 255.0f);

                png_data[(start_y + y) * out_width + (start_x + x)] = byte_val;
            }
        }
    }

    write_png(output_file, png_data, out_width, out_height);
    printf("Wrote %s\n", output_file);

    free(png_data);
    free(indices);
    free(data);

    return 0;
}
