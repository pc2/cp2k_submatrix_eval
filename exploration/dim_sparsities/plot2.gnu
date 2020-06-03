set terminal epslatex standalone size 4.5in,3.2in color

set key autotitle column nobox samplen 5 noenhanced
unset title
set output "plot.tex"
set multiplot

lw=5

col1= 57*65536 + 106*256 + 177 #(57,106,177)
col2=218*65536 + 124*256 +  48 #(218,124,48)
col3= 62*65536 + 150*256 +  81 #(62,150,81)
col4=204*65536 +  37*256 +  41 #(204,37,41)
col5= 83*65536 +  81*256 +  84 #(83,81,84)
col6=107*65536 +  76*256 + 154 #(107,76,154)
col7=146*65536 +  36*256 +  40 #(146,36,40)
col8=148*65536 + 139*256 +  61 #(148,139,61)

f="eps_1e-5.dat"

set logscale x
set logscale y
set style fill solid noborder
set tics scale 2
set key left top
#set xtics ("1" 1,"2" 2,"3" 1e-4,"$10^{-3}$" 1e-3,"$10^{-2}$" 1e-2,"$10^{-1}$" 1e-1,"$10^{-0}$" 1e-0)
#set ytics ("$10^{0}$" 1, "" 2 1,"" 3 1,"" 4 1,"" 5 1,"" 6 1,"" 7 1,"" 8 1,"" 9 1, "$10^{1}$" 10, "" 20 1, "" 30 1, "" 40 1, "" 50 1, "" 60 1, "" 70 1, "" 80 1, "" 90 1,"$10^{2}$" 100,"" 200 1,"" 300 1,"" 400 1,"" 500 1,"" 600 1,"" 700 1,"" 800 1,"" 900 1,"$10^{3}$" 1000)
set xlabel "Number of water molecules"
set ylabel "" 
set key font ",12"

set border 31 linewidth 2

#set key tmargin 2.0
#set key left bottom
set key at 100,0.07
set key spacing 2.0
set ylabel "Fraction of non-zero blocks/elements"

plot [30:15000][0.008:1] \
f u (32*($1)*$1*$1):($2/(32*($1)*$1*$1)**2) w lp lc rgb col1 lw lw pt 4 t "$\\tilde K$ SZV, block-wise", \
f u (32*($1)*$1*$1):($3/(32*($1)*$1*$1)**2) w lp lc rgb col2 lw lw pt 4 t "$\\tilde K$ DZVP, block-wise", \
f u (32*($1)*$1*$1):($6/($4)**2) w lp lc rgb col1 lw lw pt 6 dt 2 t "$SM$ SZV, block-wise", \
f u (32*($1)*$1*$1):($7/($5)**2) w lp lc rgb col2 lw lw pt 6 dt 2 t "$SM$ DZVP, block-wise", \
f u (32*($1)*$1*$1):($8/($4*6)**2) w lp lc rgb col1 lw lw pt 8 dt "." t "$SM$ SZV, element-wise", \
f u (32*($1)*$1*$1):($9/($5*23)**2) w lp lc rgb col2 lw lw pt 8 dt "." t "$SM$ DZVP, element-wise"

unset multiplot

set output "bla"
!latex plot.tex
!dvips plot.dvi
!mv plot.ps sparse.ps
!rm sparse.eps
!ps2eps sparse.ps
!rm *.dvi *.aux 

