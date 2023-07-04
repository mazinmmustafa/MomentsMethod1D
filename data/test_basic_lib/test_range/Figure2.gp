reset 
set terminal cairolatex standalone color background rgbcolor 'white' \
size 12 cm, 7.5 cm

set output 'Figure2.tex'

# set key out vert
set key reverse Left
set key top left
set key spacing 2

set xrange [1E-3:1E1]
set yrange [-0.1:+0.4]

# set ytics 200

# set mxtics 0.1
# set mytics 0.1

set grid xtics ytics mxtics mytics lt 1 lc 'grey' dt 1 lw 1

set format x '$10^{%T}$'

# set arrow from 1, -1 to 1, +1 nohead linestyle 1 lc 'black' lw 2 dt 4

set logscale x

set xlabel '$x$'
set ylabel '$f(x)$'
set title ''

plot 'data2.dat' using 1:2 with linespoints pt 7 ps 0.4 lw 2 dt 1 lc 'blue' notitle'


