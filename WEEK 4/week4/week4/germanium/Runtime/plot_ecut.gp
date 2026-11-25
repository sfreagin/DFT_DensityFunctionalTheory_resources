set terminal pngcairo size 900,600
set output 'ecut_convergence.png'

set xlabel "ecutwfc (Ry)"
set ylabel "Total energy (Ry)"
set title "Ge total energy vs ecutwfc"
set grid

set arrow from 56, graph 0 to 56, graph 1 nohead dt 2 lc rgb 'red'
set label "56 Ry" at 56, graph 0.08 right tc rgb 'red'

plot 'ecut.dat' using 1:2 with linespoints pt 7 lc rgb 'blue' title "Total energy"