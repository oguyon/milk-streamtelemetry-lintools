#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <png.h>

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

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input.fits> <N> <output.png>\n", argv[0]);
        return 1;
    }

    const char *input_file = argv[1];
    int N_slices = atoi(argv[2]);
    const char *output_file = argv[3];

    if (N_slices <= 0) {
        fprintf(stderr, "Error: N must be greater than 0.\n");
        return 1;
    }

    // Read FITS
    float *data = NULL;
    long n_frames_long = 0, pixels_per_frame = 0;
    int width, height, nframes;
    int xa, ya, naxis;

    // read_fits_float(const char *filename, float **data, long *N, long *P, int *xa, int *ya, int *naxis_out);
    // For 3D: N=naxis3 (frames), P=naxis1*naxis2 (pixels/frame). xa=naxis1, ya=naxis2.
    // read_fits_float allocates data.

    read_fits_float(input_file, &data, &n_frames_long, &pixels_per_frame, &xa, &ya, &naxis);

    if (naxis != 3) {
        fprintf(stderr, "Error: Input FITS file must be 3D cube (found NAXIS=%d).\n", naxis);
        free(data);
        return 1;
    }

    nframes = (int)n_frames_long;
    width = xa;
    height = ya;
    long pix_per_frame = (long)width * height;

    if (N_slices > nframes) {
        fprintf(stderr, "Warning: N (%d) > nframes (%d). Setting N = nframes.\n", N_slices, nframes);
        N_slices = nframes;
    }

    printf("Input: %s (%dx%dx%d)\n", input_file, width, height, nframes);
    printf("Extracting %d slices...\n", N_slices);

    // Determine indices
    int *indices = (int*) malloc(N_slices * sizeof(int));
    // Use stride
    // if N=1, take index 0.
    // if N=nframes, take all.
    // if N < nframes, take linearly spaced.
    // stride = (double)(nframes - 1) / (N_slices - 1) ?
    // Let's use simple integer arithmetic for stability or float stepping.

    if (N_slices == 1) {
        indices[0] = 0;
    } else {
        double step = (double)(nframes - 1) / (double)(N_slices - 1);
        for (int i = 0; i < N_slices; i++) {
            indices[i] = (int)(i * step + 0.5);
            if (indices[i] >= nframes) indices[i] = nframes - 1;
        }
    }

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

    if (max_val == min_val) max_val = min_val + 1.0f; // Avoid div/0

    // Construct Output Image
    // Side by side: Width = N * width, Height = height
    int out_width = N_slices * width;
    int out_height = height;

    unsigned char *png_data = (unsigned char*) calloc(out_width * out_height, sizeof(unsigned char));

    for (int i = 0; i < N_slices; i++) {
        int t = indices[i];
        long offset = t * pix_per_frame;

        // Copy slice i to png_data at x-offset i*width
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                float v = data[offset + y * width + x];
                // Normalize
                float norm = (v - min_val) / (max_val - min_val);
                if (norm < 0) norm = 0;
                if (norm > 1) norm = 1;

                unsigned char byte_val = (unsigned char)(norm * 255.0f);

                // PNG row pointers usually start from top.
                // FITS usually stores bottom-to-top or top-to-bottom depending on convention.
                // common.c likely reads it into array order.
                // Let's assume standard C array order matching write_fits_2d logic.

                int out_x = i * width + x;
                int out_y = y; // Keep same orientation

                png_data[out_y * out_width + out_x] = byte_val;
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
