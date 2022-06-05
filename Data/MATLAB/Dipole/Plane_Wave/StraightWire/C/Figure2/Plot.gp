set terminal wxt
set xrange [0:4]
set yrange [-30:10]
set xlabel("L/lambda")
set ylabel("RCS [dB]")
plot "DataRCS.dat" using 1:2 with lines lw 2 lc 16 