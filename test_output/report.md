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

## 4. CCA Analysis

Canonical Correlation Analysis was performed between Earth A and Earth B, using the top 10 PCA modes for data reduction.

**Commands:**
```bash
# Run CCA (5 vectors, using top 10 PCA modes)
../build/milk-streamtelemetry-cca -npca 10 5 earthA.fits earthB.fits
```

### Earth A Canonical Vectors (First 5)
![Earth A CCA](ccaA_vectors.png)

### Earth B Canonical Vectors (First 5)
![Earth B CCA](ccaB_vectors.png)

## 5. Visualization Commands

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

# Visualize CCA
../build/milk-streamtelemetry-3DFITS-to-png ccaAvec.fits ccaA_vectors.png -slices 0,1,2,3,4 -geom 5x1
../build/milk-streamtelemetry-3DFITS-to-png ccaBvec.fits ccaB_vectors.png -slices 0,1,2,3,4 -geom 5x1
```

# 8. Generate Corner Plot
echo "Generating Corner Plot of Canonical Variables..."

# Convert FITS variables to ASCII
"../build/milk-streamtelemetry-fits2ascii" "./ccaAvar.fits" > "./ccaAvar.txt"
"../build/milk-streamtelemetry-fits2ascii" "./ccaBvar.fits" > "./ccaBvar.txt"

# Combine side-by-side
paste "./ccaAvar.txt" "./ccaBvar.txt" > "./cca_vars.dat"

# Gnuplot script
GP_SCRIPT="./plot_cca.gp"
cat <<ENDGP > ""
set terminal pngcairo size 1000,1000 enhanced font "Arial,10"
set output "./cca_corner.png"
set multiplot layout 5,5 rowsfirst title "CCA Canonical Variables Corner Plot"

# Data file structure:
# Cols 1-5: A0..A4
# Cols 6-10: B0..B4

do for [j=0:4] {
    do for [i=0:4] {
        colA = i + 1
        colB = j + 6

        # Only show labels on edges to reduce clutter
        if (j == 4) { set xlabel sprintf("A%d", i) } else { unset xlabel }
        if (i == 0) { set ylabel sprintf("B%d", j) } else { unset ylabel }

        unset key
        # Make plots square-ish
        set size square

        plot "./cca_vars.dat" using colA:colB pt 7 ps 0.5 lc rgb "black"
    }
}
unset multiplot
ENDGP

# Run gnuplot
gnuplot ""

# Append to report
cat <<EOF >> "./report.md"

## 6. CCA Corner Plot

Scatter plots of canonical variables $ vs $. The diagonal (=j$) should show the strongest correlations.

![CCA Corner Plot](cca_corner.png)

**Gnuplot Script:**
```gnuplot

```
