reset 
set terminal cairolatex standalone color background rgbcolor 'white' \
size 12 cm, 7.5 cm

set output 'Figure1.tex'

set key out vert
set key reverse Left
set key top right 
set key spacing 2

set polar
set angle radian
set size square

set tmargin 3
set bmargin 3

set grid polar ls 1 lc 'black'

unset border
unset xtics
unset ytics
set border 0
                    
tickstep = 10                                   

r_min = -30
r_max = 0
set trange [-pi:3*pi/2]
set rrange [r_min:r_max]
set rtics tickstep

r_add = (r_max-r_min)/5
set size ratio -1

set label '$+90^{\circ}$' at (r_add+r_max-r_min)*cos(0*pi/180), (r_add+r_max-r_min)*sin(0*pi/180) center
set label '$+60^{\circ}$' at (r_add+r_max-r_min)*cos(30*pi/180), (r_add+r_max-r_min)*sin(30*pi/180) center 
set label '$+30^{\circ}$' at (r_add+r_max-r_min)*cos(60*pi/180), (r_add+r_max-r_min)*sin(60*pi/180) center
set label '$0^{\circ}$' at (r_add+r_max-r_min)*cos(90*pi/180), (r_add+r_max-r_min)*sin(90*pi/180) center
set label '$-30^{\circ}$' at (r_add+r_max-r_min)*cos(120*pi/180), (r_add+r_max-r_min)*sin(120*pi/180) center 
set label '$-60^{\circ}$' at (r_add+r_max-r_min)*cos(150*pi/180), (r_add+r_max-r_min)*sin(150*pi/180) center
set label '$-90^{\circ}$' at (r_add+r_max-r_min)*cos(180*pi/180), (r_add+r_max-r_min)*sin(180*pi/180) center
set label '$+120^{\circ}$' at (r_add+r_max-r_min)*cos(-30*pi/180), (r_add+r_max-r_min)*sin(-30*pi/180) center 
set label '$+150^{\circ}$' at (r_add+r_max-r_min)*cos(-60*pi/180), (r_add+r_max-r_min)*sin(-60*pi/180) center
set label '$\pm 180^{\circ}$' at (r_add+r_max-r_min)*cos(-90*pi/180), (r_add+r_max-r_min)*sin(-90*pi/180) center
set label '$-150^{\circ}$' at (r_add+r_max-r_min)*cos(-120*pi/180), (r_add+r_max-r_min)*sin(-120*pi/180) center 
set label '$-120^{\circ}$' at (r_add+r_max-r_min)*cos(-150*pi/180), (r_add+r_max-r_min)*sin(-150*pi/180) center

plot 'data1.dat' using (pi/2-$1):($3 > r_min ? $3 : r_min) with lines lw 2 dt 1 lc 'blue' title 'Azimuth'
