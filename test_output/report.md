# Milk Stream Telemetry Test Report

This report documents the testing of the Earth sequence generation and analysis tools.

## 1. Earth Datasets Generation

Two simulated datasets were generated representing a rotating Earth from different viewpoints.

**Commands:**
```bash
# Generate Earth A (View: Lat 0, Lon 0)
../build/milk-streamtelemetry-mkearthseq 1000 "earthA[0,0]" -speed 1.0 -size 100 100

# Generate Earth B (View: Lat 10, Lon 20)
../build/milk-streamtelemetry-mkearthseq 1000 "earthB[10,20]" -speed 1.0 -size 100 100
```

### Earth A
![Earth A Frames](earthA_frames.png)

### Earth B
![Earth B Frames](earthB_frames.png)

## 2. PCA Analysis

Principal Component Analysis was performed on both datasets (20 modes).

**Commands:**
```bash
# Run PCA
../build/milk-streamtelemetry-pca 20 earthA.fits modesA.fits coeffsA.fits
../build/milk-streamtelemetry-pca 20 earthB.fits modesB.fits coeffsB.fits
```

**Note on Mode 0:** The PCA implementation does NOT re-center (de-average) the input data. Therefore, the first mode (Mode 0) represents the average intensity of the dataset.

### Earth A Modes (First 20)
![Earth A Modes](modesA.png)

### Earth B Modes (First 20)
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

Canonical Correlation Analysis was performed between Earth A and Earth B, using the top 40 PCA modes for data reduction and computing 20 canonical vectors.

**Commands:**
```bash
# Run CCA (20 vectors, using top 40 PCA modes)
../build/milk-streamtelemetry-cca -npca 40 20 earthA.fits earthB.fits
```

### Earth A Canonical Vectors (First 20)
![Earth A CCA](ccaA_vectors.png)

### Earth B Canonical Vectors (First 20)
![Earth B CCA](ccaB_vectors.png)

## 5. Visualization Commands

The PNG images in this report were generated using the `3DFITS-to-png` tool:

```bash
# Visualize Frames
../build/milk-streamtelemetry-3DFITS-to-png earthA.fits earthA_frames.png -n 20 -geom 10x2
../build/milk-streamtelemetry-3DFITS-to-png earthB.fits earthB_frames.png -n 20 -geom 10x2

# Visualize Modes
../build/milk-streamtelemetry-3DFITS-to-png modesA.fits modesA.png -slices 0..19 -geom 10x2
../build/milk-streamtelemetry-3DFITS-to-png modesB.fits modesB.png -slices 0..19 -geom 10x2

# Visualize Reconstruction
../build/milk-streamtelemetry-3DFITS-to-png reconA.fits reconA_frames.png -n 5 -geom 5x1
../build/milk-streamtelemetry-3DFITS-to-png reconB.fits reconB_frames.png -n 5 -geom 5x1

# Visualize CCA
../build/milk-streamtelemetry-3DFITS-to-png ccaAvec.fits ccaA_vectors.png -slices 0..19 -geom 10x2
../build/milk-streamtelemetry-3DFITS-to-png ccaBvec.fits ccaB_vectors.png -slices 0..19 -geom 10x2
```


## 8. Covariance Matrix of Canonical Variables

The following table shows the covariance between canonical variables $ and $.

| A \ B | B0 | B1 | B2 | B3 | B4 | B5 | B6 | B7 | B8 | B9 | B10 | B11 | B12 | B13 | B14 | B15 | B16 | B17 | B18 | B19 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **A0** | **0.0010** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| **A1** | 0.0000 | **0.0010** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| **A2** | 0.0000 | 0.0000 | **0.0010** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| **A3** | 0.0000 | 0.0000 | 0.0000 | **0.0010** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| **A4** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | **0.0010** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| **A5** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | **0.0010** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| **A6** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | **0.0010** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| **A7** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | **0.0010** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| **A8** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | **0.0010** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| **A9** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | **0.0010** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| **A10** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | **0.0010** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| **A11** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | **0.0010** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| **A12** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | **0.0010** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| **A13** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | **0.0010** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| **A14** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | **0.0010** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| **A15** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | **0.0010** | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| **A16** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | **0.0010** | 0.0000 | 0.0000 | 0.0000 |
| **A17** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | **0.0010** | 0.0000 | 0.0000 |
| **A18** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | **0.0010** | 0.0000 |
| **A19** | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | **0.0010** |

## 6. CCA variables

Scatter plots of canonical variables $ vs $. The diagonal (=j$) should show the strongest correlations.
Correlation $ is shown on each plot. $\sigma$ values indicate standard deviation of the variables.

![CCA Variables](cca_corner.png)

### CCA Gnuplot Script
```gnuplot
set terminal pngcairo size 2000,2000 enhanced font "Arial,8"
set output "./cca_corner.png"

# Define style for boxed labels
set style textbox opaque noborder
set style fill solid 1.0 noborder

# Use margins to leave space for outer labels, spacing 0,0 to touch plots
set multiplot layout 20,20 rowsfirst title "CCA Variables (A vs B)" margins 0.04, 0.96, 0.04, 0.94 spacing 0,0

# Pre-calculate standard deviations
array stdA[20]
array stdB[20]
do for [k=1:20] {
    stats "./cca_vars.dat" u k nooutput
    stdA[k] = STATS_stddev
    stats "./cca_vars.dat" u (k+20) nooutput
    stdB[k] = STATS_stddev
}

do for [j=0:19] {
    do for [i=0:19] {
        colA = i + 1
        colB = j + 21

        # Calculate Correlation
        stats "./cca_vars.dat" u colA:colB nooutput
        corr = STATS_correlation

        # Reset labels
        unset label

        # Correlation Label (Top Right): Bold Red, White Box
        # Using a simpler label for density
        if (abs(corr) > 0.1) {
            set label 1 sprintf("%.1f", corr) at graph 0.90, 0.90 right font "Arial-Bold,7" textcolor rgb "red" front boxed
        }

        # Standard Deviation Labels
        # Top of Column
        if (j == 0) {
            set label 2 sprintf("σA%d\n%.1g", i, stdA[i+1]) at graph 0.5, 1.05 center font "Arial,7"
        }
        # Right of Row
        if (i == 19) {
            set label 3 sprintf("σB%d\n%.1g", j, stdB[j+1]) at graph 1.05, 0.5 left font "Arial,7"
        }

        # Axis Labels (Outer Only)
        set format x ""
        set format y ""
        unset xlabel
        unset ylabel
        unset xtics
        unset ytics

        # Restore tics structure but no values (or remove ticks completely?)
        # User said "remove the values on the x and y axes ticks", imply ticks might remain?
        # But for corner plot usually inner ticks are removed.
        # Let's keep ticks but remove format as requested.
        set tics scale 0.5

        if (j == 19) { set xlabel sprintf("A%d", i) font "Arial,8" }
        if (i == 0) { set ylabel sprintf("B%d", j) offset 1 font "Arial,8" }

        unset key
        # Make plots square-ish? 'set size square' can mess with layout filling.
        # With spacing 0,0, layout determines size. 'set size ratio 1' might create gaps.
        # Let's rely on layout.

        # Zero Cross
        set xzeroaxis lt 1 lc rgb "black"
        set yzeroaxis lt 1 lc rgb "black"

        plot "./cca_vars.dat" using colA:colB pt 7 ps 0.2 lc rgb "black"
    }
}
unset multiplot
```

## 7. PCA variables

Scatter plots of PCA coefficients $ vs $.

![PCA Variables](pca_corner.png)

### PCA Gnuplot Script
```gnuplot
set terminal pngcairo size 2000,2000 enhanced font "Arial,8"
set output "./pca_corner.png"

# Define style for boxed labels
set style textbox opaque noborder
set style fill solid 1.0 noborder

# Use margins to leave space for outer labels, spacing 0,0 to touch plots
set multiplot layout 20,20 rowsfirst title "PCA Variables (Coeffs A vs B)" margins 0.04, 0.96, 0.04, 0.94 spacing 0,0

# Pre-calculate standard deviations
array stdA[20]
array stdB[20]
do for [k=1:20] {
    stats "./pca_vars.dat" u k nooutput
    stdA[k] = STATS_stddev
    stats "./pca_vars.dat" u (k+20) nooutput
    stdB[k] = STATS_stddev
}

do for [j=0:19] {
    do for [i=0:19] {
        colA = i + 1
        colB = j + 21

        # Calculate Correlation
        stats "./pca_vars.dat" u colA:colB nooutput
        corr = STATS_correlation

        # Reset labels
        unset label

        # Correlation Label (Top Right): Bold Red, White Box
        # Using a simpler label for density
        if (abs(corr) > 0.1) {
            set label 1 sprintf("%.1f", corr) at graph 0.90, 0.90 right font "Arial-Bold,7" textcolor rgb "red" front boxed
        }

        # Standard Deviation Labels
        # Top of Column
        if (j == 0) {
            set label 2 sprintf("σA%d\n%.1g", i, stdA[i+1]) at graph 0.5, 1.05 center font "Arial,7"
        }
        # Right of Row
        if (i == 19) {
            set label 3 sprintf("σB%d\n%.1g", j, stdB[j+1]) at graph 1.05, 0.5 left font "Arial,7"
        }

        # Axis Labels (Outer Only)
        set format x ""
        set format y ""
        unset xlabel
        unset ylabel
        unset xtics
        unset ytics

        # Restore tics structure but no values (or remove ticks completely?)
        # User said "remove the values on the x and y axes ticks", imply ticks might remain?
        # But for corner plot usually inner ticks are removed.
        # Let's keep ticks but remove format as requested.
        set tics scale 0.5

        if (j == 19) { set xlabel sprintf("A%d", i) font "Arial,8" }
        if (i == 0) { set ylabel sprintf("B%d", j) offset 1 font "Arial,8" }

        unset key
        # Make plots square-ish? 'set size square' can mess with layout filling.
        # With spacing 0,0, layout determines size. 'set size ratio 1' might create gaps.
        # Let's rely on layout.

        # Zero Cross
        set xzeroaxis lt 1 lc rgb "black"
        set yzeroaxis lt 1 lc rgb "black"

        plot "./pca_vars.dat" using colA:colB pt 7 ps 0.2 lc rgb "black"
    }
}
unset multiplot
```

