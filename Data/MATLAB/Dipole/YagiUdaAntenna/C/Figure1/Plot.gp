set terminal wxt
set xrange [0:14]
set yrange [0:1]
set ylabel("Element number")
set ylabel("I [mA]")
plot "DataPeakCurrent.dat" with lines lw 2 lc 16 