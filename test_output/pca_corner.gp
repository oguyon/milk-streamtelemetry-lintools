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
