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
set format y "$10^{%T}$"
set yrange [1E2:2E6]
set style fill solid noborder
set tics scale 2
set key left top
#set xtics ("1" 1,"2" 2,"3" 1e-4,"$10^{-3}$" 1e-3,"$10^{-2}$" 1e-2,"$10^{-1}$" 1e-1,"$10^{-0}$" 1e-0)
#set ytics ("$10^{0}$" 1, "" 2 1,"" 3 1,"" 4 1,"" 5 1,"" 6 1,"" 7 1,"" 8 1,"" 9 1, "$10^{1}$" 10, "" 20 1, "" 30 1, "" 40 1, "" 50 1, "" 60 1, "" 70 1, "" 80 1, "" 90 1,"$10^{2}$" 100,"" 200 1,"" 300 1,"" 400 1,"" 500 1,"" 600 1,"" 700 1,"" 800 1,"" 900 1,"$10^{3}$" 1000)
set xlabel "Number of water molecules"
set ylabel "Dimension" 
set key font ",12"

set border 31 linewidth 2

#set key tmargin 2.0
#set key right top
set key at 240,1300000
set key spacing 2.0

#EREF=-74.249i194445722438
EREF=0
conv=1.0 #27.211/(32*3)*1000
plot [30:15000][:] \
f u (32*($1)*$1*$1):(6*32*($1)*$1*$1) w lp lc rgb col1 lw lw pt 6 t "$\\mathrm{dim}(\\tilde K)$ SZV-MOLOPT-SR-GTH", \
f u (32*($1)*$1*$1):(23*32*($1)*$1*$1) w lp lc rgb col2 lw lw pt 4 t "$\\mathrm{dim}(\\tilde K)$ DZVP-MOLOPT-SR-GTH", \
f u (32*($1)*$1*$1):($4*6) w lp lc rgb col1 lw lw pt 6 dt 2 t "$\\mathrm{dim}(SM)$ SZV-MOLOPT-SR-GTH", \
f u (32*($1)*$1*$1):($5*23) w lp lc rgb col2 lw lw pt 4 dt 2 t "$\\mathrm{dim}(SM)$ DZVP-MOLOPT-SR-GTH"


unset multiplot


set output "bla"
!latex plot.tex
!dvips plot.dvi
!mv plot.ps dim.ps
!rm dim.eps
!ps2eps dim.ps
!rm *.dvi *.aux 

