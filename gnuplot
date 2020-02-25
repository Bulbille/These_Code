#!/usr/bin/gnuplot
set style func linespoints
set xlabel 'T/TC'
set ylabel 'CPU-time per MC-step'
set logscale y


show style line
plot 'Interface/echelon/thermo_kaw-0.tri' u 1:($4/$3) lt 3 title 'Continous Time', \
 'DiscreteMC/local80x50/thermo_kaw-0.tri' u 1:($4/$3) lt 4 pt 1 title 'Discrete Time'

set term png
set output 'bench.png'
replot

