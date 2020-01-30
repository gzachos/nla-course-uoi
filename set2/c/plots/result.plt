#!/usr/bin/env gnuplot --persist
  
set terminal pngcairo enhanced font "arial,11" fontscale 1.0
set output 'results.png'
set bar 1.000000 front
set title 'Method Comparison'
set xlabel "Array Size"
set ylabel "Number of Iterations"
set xrange [90:11000]
set logscale x 10
set logscale y 10
set grid
set colorbox vertical origin screen 0.9, 0.2, 0 size screen 0.05, 0.6, 0 front bdefault
set style line 1 lc rgb '#d6b447' lt 1 pt 7 lw 2 ps 1.5
set style line 2 lc rgb '#02b8b8' lt 1 pt 5 lw 2 ps 1.5
set style line 3 lc rgb '#ed4d4a' lt 5 pt 9 lw 2 ps 2
set style line 4 lc rgb '#53a0ed' lt 5 pt 13 lw 2 ps 2
set key on outside bottom
plot 'plot.dat' using 1:2 title 'SD-S1' with linespoints ls 1, \
     'plot.dat' using 1:3 title 'CG-S1' with linespoints ls 2, \
     'plot.dat' using 1:4 title 'SD-S2' with linespoints ls 3, \
     'plot.dat' using 1:5 title 'CG-S2' with linespoints ls 4

