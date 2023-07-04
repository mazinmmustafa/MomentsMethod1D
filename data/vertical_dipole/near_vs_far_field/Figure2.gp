reset 
set terminal cairolatex standalone color background rgbcolor 'white' \
size 12 cm, 7.5 cm

set output 'Figure2.tex'

set key out vert
set key reverse Left
set key top right
set key spacing 2

# set xrange [0:4]
set yrange [-100:0]

# set xtics 10
set ytics 20

# set mxtics 1
# set mytics 1

# set grid xtics ytics mxtics mytics lt 1 lc 'grey' dt 1 lw 1

set format x '$10^{%T}$'

# set arrow from 1, -1 to 1, +1 nohead linestyle 1 lc 'black' lw 2 dt 4

set logscale x

set xlabel '$x/\lambda$'
set ylabel '$|H|$ [dB]'
set title ''

plot 'data2.dat' using 1:2 with lines lw 2 dt 1 lc 'blue' title 'Near Field', \
     'data2.dat' using 1:3 with lines lw 2 dt 2 lc 'red' title 'Far FIeld'
