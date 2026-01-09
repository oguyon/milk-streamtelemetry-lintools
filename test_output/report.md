# Milk Stream Telemetry Test Report

This report documents the testing of the Earth sequence generation and analysis tools.

## 1. Earth Datasets Generation

Two simulated datasets were generated representing a rotating Earth from different viewpoints.

**Commands:**
```bash
# Generate Earth A (View: Lat 0, Lon 0)
../build/milk-streamtelemetry-mkearthseq 200 "earthA[0,0]" -speed 1.0 -size 100 100

# Generate Earth B (View: Lat 30, Lon 45)
../build/milk-streamtelemetry-mkearthseq 200 "earthB[30,45]" -speed 1.0 -size 100 100
```

### Earth A
![Earth A Frames](earthA_frames.png)

### Earth B
![Earth B Frames](earthB_frames.png)

## 2. PCA Analysis

Principal Component Analysis was performed on both datasets (10 modes).

**Commands:**
```bash
# Run PCA
../build/milk-streamtelemetry-pca 10 earthA.fits modesA.fits coeffsA.fits
../build/milk-streamtelemetry-pca 10 earthB.fits modesB.fits coeffsB.fits
```

**Note on Mode 0:** The PCA implementation does NOT re-center (de-average) the input data. Therefore, the first mode (Mode 0) represents the average intensity of the dataset.

### Earth A Modes (First 5)
![Earth A Modes](modesA.png)

### Earth B Modes (First 5)
![Earth B Modes](modesB.png)

## 3. PCA Reconstruction

The datasets were reconstructed using the computed modes and coefficients.

**Commands:**
```bash
# Reconstruct
../build/milk-streamtelemetry-recon modesA.fits coeffsA.fits reconA.fits
../build/milk-streamtelemetry-recon modesB.fits coeffsB.fits reconB.fits
```

### Earth A Reconstruction
![Earth A Reconstruction](reconA_frames.png)

### Earth B Reconstruction
![Earth B Reconstruction](reconB_frames.png)

## 4. Visualization Commands

The PNG images in this report were generated using the `3DFITS-to-png` tool:

```bash
# Visualize Frames
../build/milk-streamtelemetry-3DFITS-to-png earthA.fits earthA_frames.png -n 5 -geom 5x1
../build/milk-streamtelemetry-3DFITS-to-png earthB.fits earthB_frames.png -n 5 -geom 5x1

# Visualize Modes
../build/milk-streamtelemetry-3DFITS-to-png modesA.fits modesA.png -slices 0,1,2,3,4 -geom 5x1
../build/milk-streamtelemetry-3DFITS-to-png modesB.fits modesB.png -slices 0,1,2,3,4 -geom 5x1

# Visualize Reconstruction
../build/milk-streamtelemetry-3DFITS-to-png reconA.fits reconA_frames.png -n 5 -geom 5x1
../build/milk-streamtelemetry-3DFITS-to-png reconB.fits reconB_frames.png -n 5 -geom 5x1
```
