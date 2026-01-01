# milk-streamtelemetry-cca
Canonical Correlation Analysis between two synchronized telemetry sets

## Description

Performs Canonical Correlation Analysis (CCA) between two sets of images.
The input consists of two 3D FITS cubes, noted here A and B, which are each a series of N images.
Each image is a slice of the 3D cube, so if images A have size (xa;ya) and images B have size (xb;yb), then A.fits and B.fits are 3D cubes of size (xa;ya;N) and (xb;yb;N) respectively.

The program uses CFITSIO to read the image cubes, and OpenBLAS for the linear algebra.


## Compilation

Coded in C, using gcc and cmake for compilation.
```
mkdir build
cd build
cmake ..
make
```

## Usage

```
milk-streamtelemetry-cca <nvec> A.fits B.fits
```
Where `nvec` is the number of canonical vector pairs computed.

The program supports CFITSIO's extended filename syntax, allowing you to read a subset of the image (crop) or a specific extension.

Example:
```
milk-streamtelemetry-cca 5 "A.fits[1:50,1:50,1:1000]" "B.fits[1:50,1:50,1:1000]"
```
This will perform CCA on a 50x50 pixel region of the original images. Note that for 3D cubes (x, y, t), you should specify the range for all 3 dimensions.


## Output

The program outputs the canonical vectors as images, so the output will be ccaA.fits and ccaB.fits, of sizes respecively (xa;ya;nvec) and (xb;yb;nvec).

