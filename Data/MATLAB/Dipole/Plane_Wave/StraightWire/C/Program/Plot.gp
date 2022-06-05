set terminal wxt
set xrange [0:180]
set yrange [-20:10]
plot "DataRCS.dat" using 1:2 with lines lw 2 lc 16 