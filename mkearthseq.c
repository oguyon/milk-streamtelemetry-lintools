#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Options
int width = 200;
int height = 200;
float rotation_speed = 1.0; // degrees per frame

// Texture function
float get_earth_intensity(float lat, float lon) {
    // lat: [-PI/2, PI/2], lon: [-PI, PI]
    // Simple procedural earth-like pattern

    // Base level
    float val = 0.0;

    // Continents (low freq)
    val += 0.5 * sin(2.0 * lat) * cos(lon);
    val += 0.4 * cos(3.0 * lat + 1.0) * sin(2.0 * lon + 0.5);

    // Details (high freq)
    val += 0.1 * sin(10.0 * lat) * cos(10.0 * lon);
    val += 0.05 * sin(20.0 * lat + lon);

    // Normalize roughly to 0-1
    float norm = (val + 1.05) / 2.1;
    if (norm < 0) norm = 0;
    if (norm > 1) norm = 1;

    // Threshold for land/sea
    // Land is brighter
    if (norm > 0.55) {
        return 0.6 + 0.4 * norm; // Land [0.82, 1.0]
    } else {
        return 0.1 + 0.2 * norm; // Sea [0.1, 0.21]
    }
}

void parse_filename_arg(const char *arg, char *basename, size_t basename_size, float *lat, float *lon) {
    const char *bracket = strchr(arg, '[');
    if (bracket) {
        size_t len = bracket - arg;
        if (len >= basename_size) {
            fprintf(stderr, "Error: Filename too long.\n");
            exit(1);
        }
        strncpy(basename, arg, len);
        basename[len] = '\0';

        // Parse lat,lon
        // Expected format: [lat,lon]
        if (sscanf(bracket, "[%f,%f]", lat, lon) != 2) {
             fprintf(stderr, "Error parsing coordinate in filename: %s\n", arg);
             exit(1);
        }
    } else {
        if (strlen(arg) >= basename_size) {
            fprintf(stderr, "Error: Filename too long.\n");
            exit(1);
        }
        strcpy(basename, arg);
        *lat = 0.0;
        *lon = 0.0;
    }
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <nframes> <basename[lat,lon]> [options]\n", argv[0]);
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  -size <w> <h>    Set image size (default 200 200)\n");
        fprintf(stderr, "  -speed <deg>     Set rotation speed in deg/frame (default 1.0)\n");
        return 1;
    }

    print_args(argc, argv);

    int nframes = atoi(argv[1]);
    char basename[256];
    float lat0 = 0.0, lon0 = 0.0;

    parse_filename_arg(argv[2], basename, sizeof(basename), &lat0, &lon0);

    // Parse options
    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "-size") == 0) {
            if (i + 2 < argc) {
                width = atoi(argv[i+1]);
                height = atoi(argv[i+2]);
                i += 2;
            }
        } else if (strcmp(argv[i], "-speed") == 0) {
            if (i + 1 < argc) {
                rotation_speed = atof(argv[i+1]);
                i++;
            }
        }
    }

    printf("Generating %d frames of Earth sequence.\n", nframes);
    printf("Basename: %s\n", basename);
    printf("Start Center: Lat %.2f, Lon %.2f\n", lat0, lon0);
    printf("Size: %dx%d\n", width, height);
    printf("Speed: %.2f deg/frame\n", rotation_speed);

    // Allocate 3D buffer
    // Layout: nframes * width * height
    size_t npix = (size_t)width * height * nframes;
    float *cube = (float*) malloc(npix * sizeof(float));
    if (!cube) {
        fprintf(stderr, "Failed to allocate memory for cube.\n");
        return 1;
    }

    // Precompute basis matrix M for t=0
    // M = [Right, Up, W]
    // W points to (lat0, lon0)
    double phi0 = lat0 * M_PI / 180.0;
    double lam0 = lon0 * M_PI / 180.0;

    double Wx = cos(phi0) * cos(lam0);
    double Wy = cos(phi0) * sin(lam0);
    double Wz = sin(phi0);

    // North pole
    double Nx = 0, Ny = 0, Nz = 1;

    // Up = N - (N.W)W
    double dotNW = Nz * Wz; // Nx=0, Ny=0
    double Ux = Nx - dotNW * Wx;
    double Uy = Ny - dotNW * Wy;
    double Uz = Nz - dotNW * Wz;

    // Normalize Up
    double lenU = sqrt(Ux*Ux + Uy*Uy + Uz*Uz);
    if (lenU < 1e-6) {
        // Special case: looking at pole
        if (phi0 > 0) { // North Pole
             Ux = -1; Uy = 0; Uz = 0;
        } else { // South Pole
             Ux = 1; Uy = 0; Uz = 0;
        }
        lenU = 1.0;
    }
    Ux /= lenU; Uy /= lenU; Uz /= lenU;

    // Right = Up x W
    double Rx = Uy * Wz - Uz * Wy;
    double Ry = Uz * Wx - Ux * Wz;
    double Rz = Ux * Wy - Uy * Wx;

    double rad_per_frame = rotation_speed * M_PI / 180.0;
    int min_dim = (width < height) ? width : height;

    // Generate frames
    for (int t = 0; t < nframes; t++) {
        double alpha = t * rad_per_frame;
        double cos_a = cos(alpha);
        double sin_a = sin(alpha);

        // Offset in buffer
        size_t frame_offset = (size_t)t * width * height;

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {

                double v_scale_x = 2.0 * (x - width/2.0 + 0.5) / min_dim;
                double v_scale_y = 2.0 * (y - height/2.0 + 0.5) / min_dim;

                float intensity = 0.0;

                double r2 = v_scale_x*v_scale_x + v_scale_y*v_scale_y;
                if (r2 <= 1.0) {
                    double vz = sqrt(1.0 - r2);

                    // P_I = (v_scale_x, v_scale_y, vz)
                    // P_E0 = v_scale_x * R + v_scale_y * U + vz * W
                    double px = v_scale_x * Rx + v_scale_y * Ux + vz * Wx;
                    double py = v_scale_x * Ry + v_scale_y * Uy + vz * Wy;
                    double pz = v_scale_x * Rz + v_scale_y * Uz + vz * Wz;

                    // Rotate Earth Frame (reverse rotation of Earth)
                    // P_E = R_z(-alpha) * P_E0

                    double final_x = px * cos_a + py * sin_a;
                    double final_y = -px * sin_a + py * cos_a;
                    double final_z = pz;

                    // Convert to lat, lon
                    float lat = asin(final_z);
                    float lon = atan2(final_y, final_x);

                    intensity = get_earth_intensity(lat, lon);
                }

                cube[frame_offset + y * width + x] = intensity;
            }
        }
    }

    char outname[300];
    snprintf(outname, sizeof(outname), "%s.fits", basename);
    write_fits_3d_float(outname, cube, width, height, nframes);

    free(cube);
    return 0;
}
