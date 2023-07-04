reset 
set terminal cairolatex standalone color background rgbcolor 'white' \
size 12 cm, 7.5 cm

set output 'Figure1.tex'

set key out vert
set key reverse Left
set key top right
set key spacing 2
# set key box opaque
# set border back

set xrange [1E5:1E8]
set yrange [-400:+1000]

set ytics 200

# set mxtics 5
# set mytics 5

# set grid xtics ytics mxtics mytics lt 1 lc 'grey' dt 1 lw 1

set format x '$10^{%T}$'

# set arrow from 1, -1 to 1, +1 nohead linestyle 1 lc 'black' lw 2 dt 4

set logscale x

set xlabel 'Freq [Hz]'
set ylabel '$Z=R+jX$ [ohm]'
set title ''

plot 'data1.dat' using 1:2 with lines lw 2 dt 1 lc 'blue' title '$R$ (MM)', \
     'data1.dat' using 1:3 with lines lw 2 dt 1 lc 'red' title '$X$ (MM)', \
     'data1.dat' using 1:4 with lines lw 2 dt 4 lc 'blue' title '$R$ (TL)', \
     'data1.dat' using 1:5 with lines lw 2 dt 4 lc 'red' title '$X$ (TL)'


