set terminal png size 600,400
set output "demo.png"
plot 'ecut.dat' using 1:2 with lines
